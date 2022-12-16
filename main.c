#include<stdio.h>
#include<stdlib.h>

#ifndef MPI
#define MPI
#include <mpi.h>
#endif

#ifndef LEAPFROG
#define LEAPFROG
#include "leapFrog.c"
#endif


#define NUM_OBJECTS 3
#define MASS_Earth 5.9722e24
#define MASS_SUN 1.989e30
#define MASS_JUPITER 1.898e27

#define RAD_EARTH 1.4959787e11
#define RAD_JUPITER 778e9
#define G 6.67e-11
#define V_RAD(R) sqrt((G*MASS_SUN)/R)

const double T = (365.25*24.*60.*60.);
const double dt = (24.*60.*60.);
const double nsteps = T/dt;
const double epsilon = 1e-4;

// Program goes as follows
// Each Core starts with an array of Earth, Sun, and Jupiter.
// Cores calculate the updated positions and velocities then Allgather them to eveyone 
// The last line is repeated for nsteps
// Each time the output positions, kinetic energy and potential energy is written to a file





int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);


	int num_procs;
	int proc_rank;

	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&proc_rank);

	// Set up system
	vec3 earthPos = vec3From(RAD_EARTH,0.0,0.0);
	vec3 earthVel = vec3From(0.0,V_RAD(RAD_EARTH),0.0);
	vec3 earthAcc = vec3From(0.0,0.0,0.0);
	object earth = obFrom(earthPos,earthVel,earthAcc,MASS_Earth,0.0,0.0);

	vec3 jupiterPos = vec3From(-RAD_JUPITER,0.0,0.0);
	vec3 jupiterVel = vec3From(0.0,-V_RAD(RAD_JUPITER),0.0);
	vec3 jupiterAcc = vec3From(0.0,0.0,0.0);
	object jupiter = obFrom(jupiterPos,jupiterVel,jupiterAcc,MASS_JUPITER,0.0,0.0);

	vec3 sunPos = vec3From(0.0,0.0,0.0);
	vec3 sunVel = vec3From(0.0,0.0,0.0);
	vec3 sunAcc = vec3From(0.0,0.0,0.0);
	object sun = obFrom(sunPos,sunVel,sunAcc,MASS_SUN,0.0,0.0);

	object inputArray[NUM_OBJECTS] = {earth,sun,jupiter};
	leapFrogSetup(inputArray,NUM_OBJECTS,epsilon,dt);
	// Core 0 will need to print the outputs to a file;
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
		MPI_Comm_dup(MPI_COMM_WORLD,&adjustedComm); // If none of this is applicable then set adjustedComm = M{I_COMM_WORLD
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
			countsRecv[i] = batchsize;
			displs[i] =  ((NUM_OBJECTS-1) - batchsize)*sizeof(object); // Get the displacement in bytes
			continue;
		}
		countsRecv[i] = NUM_OBJECTS / num_procs;
		displs[i] = usualBatchSize * i * sizeof(object); // Get the displacement in bytes
	}


		
	

	for (int i = 0; i < nsteps; i++){
		leapFrogStep(inputArray,NUM_OBJECTS,epsilon,dt,batchsize,startIndex);
		// TODO: MPI_AllGatherv
		MPI_Allgatherv(inputArray + (proc_rank*usualBatchSize*sizeof(object)),batchsize,MPI_CHAR,inputArray,countsRecv,displs,MPI_CHAR,adjustedComm);
	




		if (proc_rank == 0)
			writeOutputToFile(file,inputArray,NUM_OBJECTS);
	}

	if (proc_rank == 0)
		fclose(file);
	

	
	free(countsRecv); free(displs);
	
	MPI_Finalize();
	return 0;

}
