#ifndef SIMPLE_VEC
#define SIMPLE_VEC
#include "simpleVec.c"
#endif

#ifndef OBJECT
#define OBJECT
#include "Object.c"
#endif

#include <math.h>

//#ifndef MPI
//#define MPI
//#include <mpi.h>
//#endif

vec3 acceleration_and_GPE(object* const inputArray,const unsigned int inputArraySize, const unsigned int ithObject, const double epsilon){
    const double g = 6.67e-11;
    vec3 outputAcc = vec3From(0.0,0.0,0.0);
    const vec3 iPos = inputArray[ithObject].pos;
    const double iMass = inputArray[ithObject].mass;
    double GPE = 0;
    inputArray[ithObject].GPE = 0.; // Prevent build up of energy in system.
    for (unsigned int j = 0; j < inputArraySize; j++){
        if (j == ithObject){
            continue;
        }
        const vec3 jPos = inputArray[j].pos;
        const double jMass = inputArray[j].mass;

        const vec3 seperation = sub_vec3(&jPos,&iPos); // Should be the other way round but more computationally efficient to encode the minus sign here
        const double acc_mul = (g * jMass) / sqrt(pow(vec3_mag_squared(seperation) + epsilon*epsilon,3.0));
        const vec3 acc_mul_sep = scalar_mul_vec3(acc_mul,&seperation);
        outputAcc = add_vec3(&outputAcc,&acc_mul_sep);

        const double ith_GPE = (iMass * jMass) / sqrt(vec3_mag_squared(seperation));
        inputArray[ithObject].GPE += ith_GPE;

    }

    return outputAcc;
}

// updates pos, vel, KE, and GPE of every object in inputArray, then calls AllGather to update each core with the results of the operation
void leapFrogStep( object* const restrict inputArray, const unsigned int inputArraySize, const double epsilon, const double dt,const unsigned int batchSize,const unsigned int startIndex){
    
    for (int i = startIndex; i < (startIndex + batchSize); i++){ // possible segfault if startIndex + batchSize > inputArraySize
        vec3 acc = acceleration_and_GPE(inputArray,inputArraySize,i,epsilon);
        const vec3 acc_dt = scalar_mul_vec3(dt,&acc);
        const vec3 acc_half_dt = scalar_mul_vec3(0.5*dt,&acc);
        const vec3 velSyncedToPos = add_vec3(&inputArray[i].vel,&acc_half_dt); // to calculate KE
        inputArray[i].vel = add_vec3(&inputArray[i].vel,&acc_dt);

        inputArray[i].KE = 0.5 * inputArray[i].mass * vec3_mag_squared(velSyncedToPos);
        // GPE calculated in acceleration function
       

        
    }
    // Now update positions

    for (int i = startIndex; i < (startIndex + batchSize); i++){
        const vec3 v_dt = scalar_mul_vec3(dt,&inputArray[i].vel);
        inputArray[i].pos = add_vec3(&inputArray[i].pos,&v_dt);
        
    }

    

}


// updates vel by half a time step
void leapFrogSetup( object* const restrict inputArray, const unsigned int inputArraySize, const double epsilon, const double dt){
    
    for (int i = 0; i < inputArraySize; i++){ // possible segfault if startIndex + batchSize > inputArraySize
        vec3 acc = acceleration_and_GPE(inputArray,inputArraySize,i,epsilon);
        const vec3 acc_dt = scalar_mul_vec3(-0.5*dt,&acc);
        
        
        inputArray[i].vel = add_vec3(&inputArray[i].vel,&acc_dt);

        inputArray[i].KE = 0.5 * inputArray[i].mass * vec3_mag_squared(inputArray[i].vel);
        // GPE calculated in acceleration function
       

        
    }

}


void writeOutputToFile(FILE* file, const object* const restrict BUFFER, const unsigned int BUFF_LEN){
    //Calculate the total kinetic and potential energies,
    double KE = 0.;
    double GPE = 0.;
    printf("Printing %d objects: \n",BUFF_LEN);
    fflush(stdout);
    for (int i = 0; i < BUFF_LEN; i++){
        fprintf(file,"%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                BUFFER[i].mass,
                BUFFER[i].pos.x,
                BUFFER[i].pos.y,
                BUFFER[i].pos.z,
                BUFFER[i].vel.x,
                BUFFER[i].vel.y,
                BUFFER[i].vel.z,
                BUFFER[i].KE,
                BUFFER[i].GPE);
        KE += BUFFER[i].KE;
        GPE += BUFFER[i].GPE;
        fflush(file);
        

    }
    fprintf(file,"%f,%f,%f",KE,GPE,KE + GPE);
    fputc('\n',file);

}

void writeOutputArrayToFile(const object* const restrict outputArray, FILE* file, const int NUM_OBJECTS, const int nsteps){
    int n = 0;
    for (int step = 0; step < nsteps; step++){
        writeOutputToFile(file,outputArray + n,NUM_OBJECTS);
        n += NUM_OBJECTS;
    }
}

void writeOutpuDensityProfileToFile(mesh* inputMesh,FILE* file){
    double a = inputMesh->potentialCellWidth;
    for (int iLoop = 0; iLoop < inputMesh->numPotentialMeshCells; iLoop++){
       

        // Traverse linked list of cells
        double containedMass = 0.0;
        
        for (int jLoop = inputMesh->potentialMeshCells[iLoop].head; jLoop != -1; jLoop = inputMesh->objects[jLoop].nextPotentialMesh){
            

            containedMass += inputMesh->objects[jLoop].mass;
        }

        // containedMass /= vol;

        fprintf(file, "%f,%f,%f,%f,",
            containedMass,
            inputMesh->potentialMeshCells[iLoop].pos.x,
            inputMesh->potentialMeshCells[iLoop].pos.y,
            inputMesh->potentialMeshCells[iLoop].pos.z
        );
    }
    fprintf(file,"0.0\n");
}

void writeObjectsToFile(FILE* file, mesh* inputMesh){
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    object* objects = inputMesh->objects;
    int N = inputMesh->numObjects;
    for (int i = 0; i < N; i++){

        px += objects[i].mass * objects[i].vel.x;
        py += objects[i].mass * objects[i].vel.y;
        pz += objects[i].mass * objects[i].vel.z;
        const double pxTemp =  objects[i].mass * objects[i].vel.x;
        const double pyTemp =  objects[i].mass * objects[i].vel.y;
        const double pzTemp =  objects[i].mass * objects[i].vel.z;

        fprintf(file,"%f,%f,%f,%f,%f,%f,",
        objects[i].pos.x,
        objects[i].pos.y,
        objects[i].pos.z,
        pxTemp,
        pyTemp,
        pzTemp
        );

    }
    fprintf(file,"%f,%f,%f\n",px,py,pz);
}


void writeRhoArrayToFile(mesh* inputMesh, FILE* file){
   inputMesh->numPotentialMeshCellsPerSideLength; 
}