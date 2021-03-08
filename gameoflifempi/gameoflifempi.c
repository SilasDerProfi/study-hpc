#include <openmpi/mpi.h>
#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <setjmp.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

void writeParallelVTK(long timestep, int w, int h, int px, int py)
{
	int nx = w / px; 
	int ny = h / py;
	
  char filename[2048];

  snprintf(filename, sizeof(filename), "%s-%05ld%s", "parallel", timestep, ".pvti");
  FILE *fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<PImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", 0, w, 0, h, 0, 0, 1.0, 1.0, 0.0);
  fprintf(fp, "<PCellData Scalars=\"%s\">\n", "gol");
  fprintf(fp, "<PDataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", "gol");
  fprintf(fp, "</PCellData>\n");
  for (int i = 0; i < nx * ny; i++)
  {
    int xOffset = (i % nx) * px;
	int yOffset = (i / nx) * py;

    fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s-%d-%05ld%s\"/>", xOffset, xOffset + px, yOffset, yOffset + py, 0, 0, "gol", i, timestep, ".vti");
  }

  fprintf(fp, "</PImageData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void writeVTK(long timestep, double *data, char prefix[1024], int w, int h, int Gw, int offsetX, int offsetY, int num) {
  char filename[2048];  
  int x,y; 
  
  //int offsetX=0;
  //int offsetY=0;
  float deltax=1.0;
  long  nxy = w * h * sizeof(float);  

  snprintf(filename, sizeof(filename), "%s-%d-%05ld%s", prefix, num, timestep, ".vti");
  FILE* fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX, offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax, 0.0);
  fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fwrite((unsigned char*)&nxy, sizeof(long), 1, fp);

  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
        float value = data[calcIndex(Gw, x,y)];
        fwrite((unsigned char*)&value, sizeof(float), 1, fp);
    }
  }
  
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void filling(double* currentfield, int h, int w, int seed){
    srand(seed * 59 + 984);
    for (size_t i = 0; i < h * w; i++) {
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0;
    }
    return;
}

void evolve(double* currentField, double* nextField, int h, int w) {
    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        int neighbourCount = 0;
        for (int y1 = y - 1; y1 <= y + 1; y1++)
          for (int x1 = x - 1; x1 <= x + 1; x1++)
            if (x1 >= 0 && x1 < w && y1 >= -1 && y1 <= h && (x1 != x || y1 != y) && currentField[calcIndex(w, x1, y1) + w])
              neighbourCount++;

        nextField[calcIndex(w, x,y) + w] = neighbourCount == 3 || currentField[calcIndex(w, x,y) + w] && neighbourCount == 2;
      }
    }
}

void game_step(double* currentField, double* nextField, int partialHeight, int w, int size, int rank) {
    
    MPI_Request request;
    MPI_Isend(currentField + w, w, MPI_DOUBLE, ((rank - 1) + size) % size, 123, MPI_COMM_WORLD, &request);
    MPI_Isend(currentField + partialHeight * w, w, MPI_DOUBLE, ((rank + 1) + size) % size, 123, MPI_COMM_WORLD, &request);

    MPI_Status status;
    MPI_Recv(currentField, w, MPI_DOUBLE, ((rank + 1) + size) % size, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(currentField + (partialHeight + 1) * w, w, MPI_DOUBLE, ((rank - 1) + size) % size, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    evolve(currentField, nextField, partialHeight, w);
}

bool check_identical(double* currentField, double* nextField, int partialHeight, int w, int size) {
    int identical = 1;
    for (size_t y = 0; y < partialHeight; y++) {
        for (size_t x = 0; x < w; x++) {
            if(currentField[calcIndex(w, x, y + 1)] != nextField[calcIndex(w, x, y + 1)]) {
                identical = 0;
            }
        }
    }

    int allIdentical;
    MPI_Allreduce(&identical, &allIdentical, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return allIdentical == size;
}

void game(int h, int w, int size, int rank){

    int partialHeight = h / size;

    double *currentfield = calloc(w*(partialHeight + 2), sizeof(double));
    double *newfield     = calloc(w*(partialHeight + 2), sizeof(double));

    filling(currentfield + w, partialHeight, h, rank);

    for (size_t i = 0; i < 100; i++)
    {
        if(rank == 0) {
            writeParallelVTK(i, w, h, w, partialHeight);
        }

        game_step(currentfield, newfield, partialHeight, w, size, rank);

        writeVTK(i, currentfield + w, "gol", w, partialHeight, w, 0, rank * partialHeight, rank);

        if (check_identical(currentfield, newfield, partialHeight, w, size)) {
            printf("STOP\n");
            break;
        }

        double* tmpField = currentfield;
        currentfield = newfield;
        newfield = tmpField;
    }
    
}

int main(int argc, char* argv[]) {
    int rank, size, i;    

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("size %d, rank %d\n", size, rank);

    int h = 16, w = 16;

    game(h, w, size, rank);

    MPI_Finalize();
    return 0;
}