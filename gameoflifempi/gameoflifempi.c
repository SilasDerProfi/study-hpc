#include <openmpi/mpi.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
    int rank, size, i;    

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("size %d, rank %d\n", size, rank);

    MPI_Finalize();
    return 0;
}