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

void filling(double* currentfield, int partialHeight, int partialWidth, int seed){
    srand(seed * 59 + 984);
    for (size_t y = 0; y < partialHeight; y++) {
      for (size_t x = 0; x < partialWidth; x++) {
        currentfield[calcIndex(partialWidth + 2, x + 1, y + 1)] = (rand() < RAND_MAX / 10) ? 1 : 0;
      }
    }
    return;
}

void evolve(double* currentField, double* nextField, int h, int w) {
    for (int y = 1; y < h + 1; y++) {
      for (int x = 1; x < w + 1; x++) {
        int neighbourCount = 0;
        for (int y1 = y - 1; y1 <= y + 1; y1++)
          for (int x1 = x - 1; x1 <= x + 1; x1++)
            if (x1 >= 0 && x1 < w && y1 >= -1 && y1 <= h && (x1 != x || y1 != y) && currentField[calcIndex(w + 2, x1, y1)])
              neighbourCount++;
        nextField[calcIndex(w + 2, x,y)] = neighbourCount == 3 || currentField[calcIndex(w + 2, x,y)] && neighbourCount == 2;
      }
    }
}

void game_step(double* currentField, double* nextField, int partialHeight, int partialWidth, int rank,  MPI_Comm comm_cart) {
    
    int upperNeighbor, lowerNeighbor, leftNeighbor, rightNeighbor;

    MPI_Cart_shift(comm_cart, 1, 1, &upperNeighbor, &lowerNeighbor);
    MPI_Cart_shift(comm_cart, 0, 1, &leftNeighbor, &rightNeighbor);

    int w = partialWidth + 2;
    int h = partialHeight + 2;

    int gSizes[2] = {h, w};
    int lSizesX[2] = {h, 1};
    int lSizesY[2] = {1, w};

    MPI_Datatype gl_left_outer, gl_right_outer, gl_upper_outer, gl_lower_outer;
    MPI_Datatype gl_left_inner, gl_right_inner, gl_upper_inner, gl_lower_inner;

    int startIndices[2] = {0,0};
    MPI_Type_create_subarray(2, gSizes, lSizesX, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_left_outer);
    MPI_Type_commit(&gl_left_outer);

    startIndices[1] = 1;
    MPI_Type_create_subarray(2, gSizes, lSizesX, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_left_inner);
    MPI_Type_commit(&gl_left_inner);

    startIndices[1] = w - 1;
    MPI_Type_create_subarray(2, gSizes, lSizesX, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_right_outer);
    MPI_Type_commit(&gl_right_outer);

    startIndices[1] = w - 2;
    MPI_Type_create_subarray(2, gSizes, lSizesX, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_right_inner);
    MPI_Type_commit(&gl_right_inner);

    startIndices[0] = 0;
    startIndices[1] = 0;
    MPI_Type_create_subarray(2, gSizes, lSizesY, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_upper_outer);
    MPI_Type_commit(&gl_upper_outer);

    startIndices[0] = 1;
    MPI_Type_create_subarray(2, gSizes, lSizesY, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_upper_inner);
    MPI_Type_commit(&gl_upper_inner);

    startIndices[0] = h - 1;
    MPI_Type_create_subarray(2, gSizes, lSizesY, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_lower_outer);
    MPI_Type_commit(&gl_lower_outer);

    startIndices[0] = h - 2;
    MPI_Type_create_subarray(2, gSizes, lSizesY, startIndices, MPI_ORDER_C, MPI_DOUBLE, &gl_lower_inner);
    MPI_Type_commit(&gl_lower_inner);

    MPI_Request request;
    MPI_Isend(currentField, 1, gl_left_inner, leftNeighbor, 123, comm_cart, &request);
    MPI_Isend(currentField, 1, gl_right_inner, rightNeighbor, 123, comm_cart, &request);
    MPI_Isend(currentField, 1, gl_upper_inner, upperNeighbor, 123, comm_cart, &request);
    MPI_Isend(currentField, 1, gl_lower_inner, lowerNeighbor, 123, comm_cart, &request);

    MPI_Status status;
    MPI_Recv(currentField, 1, gl_left_outer, leftNeighbor, MPI_ANY_TAG, comm_cart, &status);
    MPI_Recv(currentField, 1, gl_right_outer, rightNeighbor, MPI_ANY_TAG, comm_cart, &status);
    MPI_Recv(currentField, 1, gl_upper_outer, upperNeighbor, MPI_ANY_TAG, comm_cart, &status);
    MPI_Recv(currentField, 1, gl_lower_outer, lowerNeighbor, MPI_ANY_TAG, comm_cart, &status);

    evolve(currentField, nextField, partialHeight, partialWidth);
}

bool check_identical(double* currentField, double* nextField, int partialHeight, int partialWidth, int size) {
    int identical = 1;
    for (size_t y = 0; y < partialHeight; y++) {
        for (size_t x = 0; x < partialWidth; x++) {
            if(currentField[calcIndex(partialWidth, x + 1, y + 1)] != nextField[calcIndex(partialWidth, x + 1, y + 1)]) {
                identical = 0;
            }
        }
    }

    int allIdentical;
    MPI_Allreduce(&identical, &allIdentical, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return allIdentical == size;
}

void game(int h, int w, int* dims, int rank, int size, MPI_Comm comm_cart){

    int partialHeight = h / dims[1];
    int partialWidth = w / dims[0];

    double *currentfield = calloc((partialWidth + 2) * (partialHeight + 2), sizeof(double));
    double *newfield     = calloc((partialWidth + 2) * (partialHeight + 2), sizeof(double));

    filling(currentfield, partialHeight, partialWidth, rank);

    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);

    for (size_t i = 0; i < 100; i++)
    {
        if(rank == 0) {
            writeParallelVTK(i, w, h, partialWidth, partialHeight);
        }

        game_step(currentfield, newfield, partialHeight, partialWidth, rank, comm_cart);

        int num = coords[0] + coords[1] * dims[0];
        writeVTK(i, currentfield + w, "gol", partialWidth, partialHeight, w, coords[0] * partialWidth, coords[1] * partialHeight, num);

        // if (check_identical(currentfield, newfield, partialHeight, partialWidth, size)) {
        //     printf("STOP\n");
        //     break;
        // }

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

    int dims[2];
    dims[0] = 0;
    dims[1] = 0;
    MPI_Dims_create(size, 2, dims);

    int periods[2] = {true, true};

    MPI_Comm comm_cart;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &comm_cart);

    int h = 16, w = 16;

    game(h, w, dims, rank, size, comm_cart);

    MPI_Comm_free(&comm_cart);

    MPI_Finalize();
    return 0;
}