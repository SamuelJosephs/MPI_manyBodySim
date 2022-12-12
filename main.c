#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>


#ifndef LEAPFROG
#define LEAPFROG
#include "leapFrog.c"
#endif


#define NUM_OBJECTS 1000





int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int num_procs;
	int proc_rank;

	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&proc_rank);

	printf("Proc %d says Hello!\n",proc_rank);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;

}
