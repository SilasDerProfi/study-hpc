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

void writeVTK_MPI(size_t timestep, char prefix[1024], double* currentfield, int threadX, int threadY, int partialWidth, int partialHeight, int width, int height, int rank) {

  long nxy = width * height * sizeof(double); 

  MPI_File file;

  char filename[2048];
  snprintf(filename, sizeof(filename), "%s-%ld.vti", prefix, timestep);

  char header[2048];
  snprintf(header, sizeof(header), 
    "<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
    "<ImageData WholeExtent=\"%d %d %d %d 0 0\" Origin=\"0 0 0\" Spacing=\"1.0 1.0 0.0\">\n"
    "<CellData Scalars=\"%s\">\n"
    "<DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n"
    "</CellData>\n"
    "</ImageData>\n"
    "<AppendedData encoding=\"raw\">\n"
    "_", 0, width, 0, height, prefix, prefix);

  char footer[128];
  snprintf(footer, sizeof(footer), "\n</AppendedData>\n</VTKFile>");

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);

  int headerSize = strlen(header);
  int footerSize = strlen(footer);

  if(rank == 0) {
    int contentLength = width * height * sizeof(double);
    MPI_File_write(file, header, headerSize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_write(file, &nxy, 1, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_seek(file, headerSize + sizeof(nxy) + contentLength, MPI_SEEK_SET);
    MPI_File_write(file, footer, footerSize, MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_Datatype fileType, fieldType;
  int fieldSizes[2] = {height, width};
  int partialSizes[2] = {partialHeight, partialWidth};
  int memorySizes[2] = {partialHeight + 2, partialWidth + 2};
  int startIndices[2];

  startIndices[0] = threadY * partialHeight;
  startIndices[1] = threadX * partialWidth;
  MPI_Type_create_subarray(2, fieldSizes, partialSizes, startIndices, MPI_ORDER_C, MPI_DOUBLE, &fileType);
  MPI_Type_commit(&fileType);

  startIndices[0] = 1;
  startIndices[1] = 1;
  MPI_Type_create_subarray(2, memorySizes, partialSizes, startIndices, MPI_ORDER_C, MPI_DOUBLE, &fieldType);
  MPI_Type_commit(&fieldType);

  MPI_File_set_view(file, headerSize + sizeof(nxy), MPI_DOUBLE, fileType, "native", MPI_INFO_NULL);
  MPI_File_write_all(file, currentfield, 1, fieldType, MPI_STATUS_IGNORE);

  MPI_File_close(&file);
}

void show(double* currentfield, int w, int h) {
  printf("\033[H");
  int x,y;
  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
    //printf("\033[E");
    printf("\n");
  }
  fflush(stdout);
}

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

  for (y = 1; y <= h; y++) {
    for (x = 1; x <= w; x++) {
      float value = data[calcIndex(w + 2, x,y)];
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
        currentfield[calcIndex(partialWidth + 2, x + 1, y + 1)] = (rand() < RAND_MAX / 10)? 1: 0;

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
            if ((x1 != x || y1 != y) && currentField[calcIndex(w + 2, x1, y1)])
              neighbourCount++;
        nextField[calcIndex(w + 2, x,y)] = neighbourCount == 3 || currentField[calcIndex(w + 2, x,y)] && neighbourCount == 2;
        // nextField[calcIndex(w + 2, x,y)] = neighbourCount;
      }
    }
}

void game_step(double* currentField, double* nextField, int partialHeight, int partialWidth, int rank,  MPI_Comm comm_cart) {
    
    int upperNeighbor, lowerNeighbor, leftNeighbor, rightNeighbor;

    MPI_Cart_shift(comm_cart, 0, 1, &leftNeighbor, &rightNeighbor);
    MPI_Cart_shift(comm_cart, 1, 1, &upperNeighbor, &lowerNeighbor);

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
    MPI_Recv(currentField, 1, gl_right_outer, leftNeighbor, MPI_ANY_TAG, comm_cart, &status);
    MPI_Recv(currentField, 1, gl_left_outer, rightNeighbor, MPI_ANY_TAG, comm_cart, &status);
    MPI_Recv(currentField, 1, gl_lower_outer, upperNeighbor, MPI_ANY_TAG, comm_cart, &status);
    MPI_Recv(currentField, 1, gl_upper_outer, lowerNeighbor, MPI_ANY_TAG, comm_cart, &status);

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

void game(int time_steps, int h, int w, int* dims, int rank, int size, MPI_Comm comm_cart){

    int partialHeight = h / dims[1];
    int partialWidth = w / dims[0];

    double *currentfield = calloc((partialWidth + 2) * (partialHeight + 2), sizeof(double));
    double *newfield     = calloc((partialWidth + 2) * (partialHeight + 2), sizeof(double));

    filling(currentfield, partialHeight, partialWidth, rank);

    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);

    for (size_t i = 0; i < time_steps; i++)
    {
        game_step(currentfield, newfield, partialHeight, partialWidth, rank, comm_cart);

        // writeVTK_MPI(i, "gol", currentfield, coords[0], coords[1], partialWidth, partialHeight, w, h, rank);

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

    int tS = 32, nX = 8, nY = 8, pX = 2, pY = 2;

    if (argc > 1) tS = atoi(argv[1]); ///< read nX

    if (argc > 2) nX = atoi(argv[2]); ///< read nX
    if (argc > 3) nY = atoi(argv[3]); ///< read nY

    if (argc > 4) pX = atoi(argv[4]); ///< read pX
    if (argc > 5) pY = atoi(argv[5]); ///< read pY

    int h = nY * pY;
    int w = nX * pX;

    // printf("size %d, rank %d\n", size, rank);

    int dims[2];
    if(argc > 5) {
      dims[0] = pY;
      dims[1] = pX;
    }else {
      dims[0] = 0;
      dims[1] = 0;
      MPI_Dims_create(size, 2, dims);
    }

    int periods[2] = {true, true};

    MPI_Comm comm_cart;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &comm_cart);

    game(tS, h, w, dims, rank, size, comm_cart);

    MPI_Comm_free(&comm_cart);

    MPI_Finalize();
    return 0;
}