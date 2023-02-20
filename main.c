#include "common.h"

#ifndef MPI
#define MPI
#include <mpi.h>
#endif

#ifndef LEAPFROG
#define LEAPFROG
#include "leapFrog.c"
#endif

#include "systems.c"


const double T = 4000*(365.25*24.*60.*60.);
const double dt = 10*(24.*60.*60.);
const double nsteps = T/dt;
const double epsilon = 1e3;

// Program goes as follows
// Each Core starts with an array of Earth, Sun, and Jupiter.
// Cores calculate the updated positions and velocities then Allgather them to eveyone 
// The last line is repeated for nsteps
// Each time the output positions, kinetic energy and potential energy is written to a file





int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);


	int num_procs;
	int proc_rank;
	int NUM_OBJECTS;

	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&proc_rank);

	// Set up system
	// object *inputArray = EarthSunJupiter(&NUM_OBJECTS,epsilon,dt);
	object* inputArray = sphericalDistribution(&NUM_OBJECTS,epsilon,dt);
	object* outputArray = malloc(NUM_OBJECTS * nsteps * sizeof(object));
	leapFrogSetup(inputArray,NUM_OBJECTS,epsilon,dt);
	FILE* file;
	if (proc_rank == 0){
		file = fopen("output.csv","w+");
	}
	if (proc_rank == 0){
		writeOutputToFile(file,inputArray,NUM_OBJECTS);
	}
	// deal with the edge case where num_procs > NUM_OBJECTS, use adjusted comm for all further MPI calls.
	MPI_Comm adjustedComm;
	MPI_Group worldGroup;
	MPI_Group adjustedGroup;
	
	printf("There are %d procs\n",num_procs);
	if (num_procs > NUM_OBJECTS){
		
		// Calculate ranks in worldGroup to exclude from adjustedGroup
		const int diff = num_procs - NUM_OBJECTS; 
		int ranks[diff]; // might need to malloc
		int counter = 0;
		for (int i = 0; i < num_procs; i++){ // i is proc_rank in MPI_COMM_WORLD
			if (i >= NUM_OBJECTS){
				ranks[counter] = i;
				counter++;
			}
		}

		// Create new group (adjustedGroup) excluding those ranks
		MPI_Comm_group(MPI_COMM_WORLD,&worldGroup);
		MPI_Group_excl(worldGroup,diff,ranks,&adjustedGroup);
		const int result = MPI_Comm_create(MPI_COMM_WORLD,adjustedGroup,&adjustedComm);
		if (result == MPI_SUCCESS)
			fprintf(stderr, "Result was a succsess!\n");
		else {fprintf(stderr, "Result was an error\n");};
		
		// recalculate the number of procs and ranks (Only for procs in adjustedComm)
		if (proc_rank < NUM_OBJECTS ){
			MPI_Comm_size(adjustedComm,&num_procs);
			MPI_Comm_rank(adjustedComm,&proc_rank);
		} else {
			// Else just exit the program.
			MPI_Finalize();
			fprintf(stderr,"Core %i is uneeded and will terminate normally",proc_rank);
			exit(0);
		}

	} else {
		MPI_Comm_dup(MPI_COMM_WORLD,&adjustedComm); // If none of this is applicable then set adjustedComm = MPI_COMM_WORLD
	}


	unsigned int batchsize = NUM_OBJECTS / num_procs;
	const unsigned int usualBatchSize = batchsize; // Used for calculating displacements
	unsigned int startIndex = batchsize * proc_rank;
	
	// If NUM_OBJECTS cannot be divided up nicely between cores then let the final core have the remainder
	// The possiblity of num_cores > NUM_OBJECTS was eliminated above
	if ((NUM_OBJECTS % num_procs != 0) && (proc_rank == num_procs - 1)){
		batchsize = NUM_OBJECTS % num_procs;
		startIndex = (NUM_OBJECTS - 1) - batchsize; // final index - batchsize
	}

	

	int* countsRecv = malloc(num_procs *sizeof(int));
	int* displs = malloc(num_procs * sizeof(int));
	// Compute the number of counts recieved from each processor and where they should go;
	for (int i = 0; i < num_procs; i++){
		if (i == num_procs - 1){
			countsRecv[i] = batchsize*sizeof(object);
			displs[i] =  usualBatchSize*i*sizeof(object); // Get the displacement in bytes
			continue;
		}
		countsRecv[i] = usualBatchSize*sizeof(object);
		displs[i] = usualBatchSize * i * sizeof(object); // Get the displacement in bytes
	}


		
	
	printf("StartIndex from proc_rank %d is %d with batchSize %d\n ",proc_rank,startIndex,batchsize);
	int n = 0;
	int counter = 0;
	for (int i = 0; i < nsteps; i++){
		// leapFrogStep(inputArray,NUM_OBJECTS,epsilon,dt,batchsize,startIndex);
		leapFrogStep(inputArray,NUM_OBJECTS,epsilon,dt,batchsize,startIndex);
		
		// TODO: MPI_AllGatherv
		MPI_Allgatherv(inputArray + startIndex,batchsize*sizeof(object),MPI_CHAR,inputArray,countsRecv,displs,MPI_CHAR,adjustedComm);
		// Each core copy batchsize objects to outputArray + n + startIndex

		// memcpy(outputArray + n + startIndex, inputArray,batchsize * sizeof(object));
		// Write contents to outputArray.
		MPI_Allgatherv(inputArray + startIndex,batchsize*sizeof(object),MPI_CHAR,outputArray + n,countsRecv,displs,MPI_CHAR,adjustedComm);
		n += NUM_OBJECTS;

		//TODO: Impliment dynamic formatable strings




	}


		if (proc_rank == 0)
			writeOutputArrayToFile(outputArray,file,NUM_OBJECTS,nsteps);

	if (proc_rank == 0)
		fclose(file);
	

	
	free(countsRecv); free(displs); free(inputArray); free(outputArray);
	
	MPI_Finalize();
	return 0;

}
