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

void evolve(){
    
}

void filling(int32_t* currentfield, int h, int w){
    int32_t i;
    for (i = 0; i < h * w; i++) {
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0;
    }
    return;
}

void game_step(int32_t* currentField, int partialHeight, int w, int size, int rank) {
    
    MPI_Request request;
    MPI_Isend(currentField + w, w, MPI_INT32_T, ((rank - 1) + size) % size, 123, MPI_COMM_WORLD, &request);
    MPI_Isend(currentField + partialHeight * w, w, MPI_INT32_T, ((rank + 1) + size) % size, 123, MPI_COMM_WORLD, &request);

    MPI_Status status;
    MPI_Recv(currentField, w, MPI_INT32_T, ((rank + 1) + size) % size, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(currentField + (partialHeight + 1) * w, w, MPI_INT32_T, ((rank - 1) + size) % size, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
}

void game(int h, int w, int size, int rank){

    int partialHeight = h / size;

    int32_t *currentfield = calloc(w*(partialHeight + 2), sizeof(int32_t));
    int32_t *newfield     = calloc(w*(partialHeight + 2), sizeof(int32_t));
  
    filling(currentfield + w, w, partialHeight);

    for (size_t i = 0; i < 10; i++)
    {
        game_step(currentfield, partialHeight, w, size, rank);

        int32_t* tmpField = currentfield;
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

    int h = 32, w = 32;

    game(h, w, size, rank);

    MPI_Finalize();
    return 0;
}