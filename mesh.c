#include <math.h>

#ifndef SIMPLE_VEC
#define SIMPLE_VEC
#include "simpleVec.c"
#endif

#ifndef OBJECT
#define OBJECT
#include "Object.c"
#endif

#include "common.h"
#include <complex.h>
#include "fftw3.h"
// TODO: Check if object moves to a different cell and impose periodic boundary's 
// I.E.: IF an object leaves the domain of the "universe" it appears in the opposite side
typedef struct meshCell {
    int head; // head of chain
    vec3 pos; // position of center of cell;
    double rho;
    double volume;

} meshCell;

typedef struct mesh {
    unsigned int numMeshCells;
    unsigned int numMeshCellsPerSideLength;
    unsigned int numPotentialMeshCellsPerSideLength;
    unsigned int numPotentialMeshCells;
    meshCell* meshCells;
    meshCell* potentialMeshCells;
    double cellWidth;
    double potentialCellWidth;
    double universeWidth;
    object* objects;
    unsigned int numObjects;
    double epsilon;
    int numIterationsPast;
    fftw_complex* rho_array;
    fftw_plan plan;
    fftw_complex* kxArray;
    fftw_plan kxplan;
    fftw_complex* kyArray;
    fftw_plan kyplan;
    fftw_complex* kzArray;
    fftw_plan kzplan;
} mesh;


int indexArray(int i, int j, int k, int N){
    return ((N*N*i) + N*j + k);
}
meshCell* indexMesh(const mesh* inputMesh, int i, int j, int k){
    const int N = inputMesh->numMeshCellsPerSideLength;
    return &(inputMesh->meshCells[((i * N*N) + (j*N) + k)]); // Row major order
}
meshCell* indexPotentialMesh(const mesh* inputMesh, int i, int j, int k){
    const int N = inputMesh->numPotentialMeshCellsPerSideLength;
    return &(inputMesh->potentialMeshCells[((i * N*N) + (j*N) + k)]);
}

fftw_complex* indexRhoArray(const mesh* inputMesh, int i, int j, int k){
    const int N = inputMesh->numPotentialMeshCellsPerSideLength;
    return &(inputMesh->rho_array[((i*N*N) + (j*N) + k)]);
    
}
void assignObjectToMeshCell(int obNum,mesh* inputMesh,int* iOut, int* jOut, int* kOut){
    // Work out fraction across x,y, and z axis
    object ob = inputMesh->objects[obNum];
    double width = inputMesh->universeWidth;
    int numCells = inputMesh->numMeshCellsPerSideLength;
    double fracX = ob.pos.x / width;
    double fracY = ob.pos.y / width;
    double fracZ = ob.pos.z / width;

    // Use this to work out where it should be in the Mesh
    int i = (int)(fracX * (double)numCells -0.5);
    int j = (int)(fracY * (double)numCells - 0.5);
    int k = (int)(fracZ * (double)numCells - 0.5);

    if (iOut != NULL || jOut != NULL || kOut != NULL){
        *iOut = i;
        *jOut = j;
        *kOut = k;
        
        return;
    }

    meshCell* cell = indexMesh(inputMesh,i,j,k);
    // Prepend to linked list in Cell
    inputMesh->objects[obNum].next = cell->head;
    cell->head = obNum;
    if (cell->head == inputMesh->objects[cell->head].next){
        fprintf(stderr,"cell->head = objects[cell->head].next, this will result in an infite loop when traversing the list\n");
    }
}

void assignObjectToPotentialMeshCell(int obNum,mesh* inputMesh,int* iOut, int* jOut, int* kOut){
    // Work out fraction across x,y, and z axis
    object ob = inputMesh->objects[obNum];
    double width = inputMesh->universeWidth;
    int numCells = inputMesh->numPotentialMeshCellsPerSideLength;
    double fracX = ob.pos.x / width;
    double fracY = ob.pos.y / width;
    double fracZ = ob.pos.z / width;

    // Use this to work out where it should be in the Mesh
    int i = (int)(fracX * (double)numCells -0.5);
    int j = (int)(fracY * (double)numCells - 0.5);
    int k = (int)(fracZ * (double)numCells - 0.5);

    if (iOut != NULL || jOut != NULL || kOut != NULL){
        *iOut = i;
        *jOut = j;
        *kOut = k;
        
        return;
    }
    meshCell* cell = indexPotentialMesh(inputMesh,i,j,k);
    // Prepend to linked list in Cell
    fprintf(stderr, "Potential: Cell Head = %d, next: %d\n",cell->head,inputMesh->objects[cell->head].nextPotentialMesh);
    inputMesh->objects[obNum].nextPotentialMesh = cell->head;
    cell->head = obNum;
    if (cell->head == inputMesh->objects[obNum].nextPotentialMesh){
        fprintf(stderr,"Potential: Cell head == obnum.nextPotentialMesh, %d,%d This will cause an infinite loop\n",cell->head,inputMesh->objects[obNum].nextPotentialMesh);
        // exit(1);
    }
}

void printMeshCellObjects(meshCell* inputCell,mesh* inputMesh){
    
    for (int i = inputCell->head; i != -1; i = inputMesh->objects[i].next){
        if (i > inputMesh->numObjects){
            break;
        }
        printf("i = %d, next = %d\n",i,inputMesh->objects[i].next);
        printf("%i Index:(%d,%d,%d), ",i,inputMesh->objects[i].i,inputMesh->objects[i].j,inputMesh->objects[i].k);
        fflush(stdout);
    }
    printf("\nEnd of objects\n");
}

void printPotentialMeshCellObjects(const meshCell* inputCell,const mesh* inputMesh){
    
    for (int itemp = inputCell->head; itemp != -1; itemp = inputMesh->objects[itemp].nextPotentialMesh){
        printf("i = %d, next = %d\n",itemp,inputMesh->objects[itemp].next);
        printf("%i Index:(%d,%d,%d), ",itemp,inputMesh->objects[itemp].i,inputMesh->objects[itemp].j,inputMesh->objects[itemp].k);
    }
    printf("\nEnd of objects\n");
}

void findObjectInMesh(mesh* inputMesh, int objectNum,int* iOut, int* jOut, int*kOut){
    const int numCells = inputMesh->numMeshCellsPerSideLength;
    for (int i = 0; i < numCells; i++){
        for (int j = 0; j < numCells; j++){
            for (int k = 0; k < numCells; k++){
                meshCell* cell = indexMesh(inputMesh,i,j,k);
                for (int l = cell->head ; ;l = inputMesh->objects[l].next){
                    if (l == objectNum){
                        *iOut = i;
                        *jOut = j;
                        *kOut = k;
                        return;
                    }
                    if (l == -1){
                        break;
                    }
                }
            }
        }
    }
    return;
}

void findObjectInPotentialMesh(mesh* inputMesh, int objectNum,int* iOut, int* jOut, int*kOut){
    const int numCells = inputMesh->numPotentialMeshCellsPerSideLength;
    for (int i = 0; i < numCells; i++){
        for (int j = 0; j < numCells; j++){
            for (int k = 0; k < numCells; k++){
                meshCell* cell = indexPotentialMesh(inputMesh,i,j,k);
                for (int l = cell->head ; ;l = inputMesh->objects[l].nextPotentialMesh){
                    if (l == objectNum){
                        *iOut = i;
                        *jOut = j;
                        *kOut = k;
                        return;
                    }
                    if (l == -1){
                        break;
                    }
                }
            }
        }
    }
    return;
}

void assignObjectsToMesh(object* objects,const unsigned int numObjects, mesh* inputMesh){
    for (int i = 0; i < numObjects; i++){
        int iPos;
        int kPos;
        int jPos;
        assignObjectToMeshCell(i,inputMesh,&iPos,&jPos,&kPos);
        objects[i].i = iPos;
        objects[i].j = jPos;
        objects[i].k = kPos;
        assignObjectToMeshCell(i,inputMesh,NULL,NULL,NULL);

        int iPosPotential;
        int jPosPotential;
        int kPosPotential;
        assignObjectToPotentialMeshCell(i,inputMesh,&iPosPotential,&jPosPotential,&kPosPotential);
        objects[i].iPotential = iPosPotential;
        objects[i].jPotential = jPosPotential;
        objects[i].kPotential = kPosPotential;
        assignObjectToPotentialMeshCell(i,inputMesh,NULL,NULL,NULL);

    }
    
}

mesh meshFrom(const double meshCellWidth, const double universeWidth, int numPotentialMeshCellsPerMeshCell, object* objects, int numObjects,double epsilon){
    mesh output;
    output.numIterationsPast = 0;
    output.objects = objects;
    output.numObjects = numObjects;
    output.universeWidth = universeWidth;
    output.numMeshCellsPerSideLength = (unsigned int) (universeWidth / meshCellWidth);
    output.numPotentialMeshCellsPerSideLength = (unsigned int) output.numMeshCellsPerSideLength * numPotentialMeshCellsPerMeshCell;
    output.numPotentialMeshCells = output.numPotentialMeshCellsPerSideLength * output.numPotentialMeshCellsPerSideLength * output.numPotentialMeshCellsPerSideLength; // cube it 
     
    output.numMeshCells = output.numMeshCellsPerSideLength * output.numMeshCellsPerSideLength * output.numMeshCellsPerSideLength; // cube it 
    output.cellWidth = meshCellWidth;
    output.potentialCellWidth = universeWidth / output.numPotentialMeshCellsPerSideLength;
    output.epsilon = epsilon;

    output.meshCells = malloc(output.numMeshCells * sizeof(meshCell));
    output.potentialMeshCells = malloc(output.numPotentialMeshCells * sizeof(meshCell));
    output.rho_array = fftw_malloc(output.numPotentialMeshCells * sizeof(fftw_complex));
    output.kxArray = fftw_malloc(output.numPotentialMeshCells * sizeof(fftw_complex));
    output.kyArray = fftw_malloc(output.numPotentialMeshCells * sizeof(fftw_complex));
    output.kzArray = fftw_malloc(output.numPotentialMeshCells * sizeof(fftw_complex));



    int N = output.numPotentialMeshCellsPerSideLength;
    printf("Setting up fftw plans\n");
    output.plan = fftw_plan_dft_3d(N,N,N,output.rho_array,output.rho_array,FFTW_FORWARD,FFTW_MEASURE);
    output.kxplan = fftw_plan_dft_3d(N,N,N,output.kxArray,output.kxArray,FFTW_BACKWARD,FFTW_MEASURE);
    output.kyplan = fftw_plan_dft_3d(N,N,N,output.kyArray,output.kyArray,FFTW_BACKWARD,FFTW_MEASURE);
    output.kzplan = fftw_plan_dft_3d(N,N,N,output.kzArray,output.kzArray,FFTW_BACKWARD,FFTW_MEASURE);



    if (output.meshCells == NULL){
        fprintf(stderr,"Failed to allocate meshCells\n");
        exit(1);
    }
    meshCell temp;
    for (int i = 0; i < output.numMeshCells; i++){ // initialise next pointers to each meshcell to be -1 (empty)
        
        output.meshCells[i] = temp; // Maybe undefined behaviour. Initialise heads to be empty (-1).
        output.meshCells[i].head = -1;
        
    }
    double vol = output.potentialCellWidth * output.potentialCellWidth * output.potentialCellWidth; 
    for (int i = 0; i < output.numPotentialMeshCells; i++){
        output.potentialMeshCells[i] = temp;
        output.potentialMeshCells[i].head = -1;
        output.potentialMeshCells[i].volume = vol;
    }
    printf("Got Past Initialisation of meshcells\n");

    for (int i = 0; i < output.numMeshCellsPerSideLength; i ++) {
        for (int j = 0; j < output.numMeshCellsPerSideLength; j++){
            for (int k = 0; k < output.numMeshCellsPerSideLength; k ++){
                vec3 pos;
                pos.x = (((double) i) / ((double) output.numMeshCellsPerSideLength)) * output.universeWidth - (0.5 * meshCellWidth); 
                pos.y = (((double) j) / ((double) output.numMeshCellsPerSideLength)) * output.universeWidth - (0.5 * meshCellWidth);
                pos.z = (((double) k) / ((double) output.numMeshCellsPerSideLength)) * output.universeWidth - (0.5 * meshCellWidth);
                indexMesh(&output,i,j,k)->pos = pos;
            }
        }
    }
    for (int i = 0; i < output.numPotentialMeshCellsPerSideLength; i ++) {
        for (int j = 0; j < output.numPotentialMeshCellsPerSideLength; j++){
            for (int k = 0; k < output.numPotentialMeshCellsPerSideLength; k ++){
                vec3 pos;
                pos.x = (((double) i) / ((double) output.numPotentialMeshCellsPerSideLength)) * output.universeWidth - (0.5 * output.potentialCellWidth); 
                pos.y = (((double) j) / ((double) output.numPotentialMeshCellsPerSideLength)) * output.universeWidth - (0.5 * output.potentialCellWidth);
                pos.z = (((double) k) / ((double) output.numPotentialMeshCellsPerSideLength)) * output.universeWidth - (0.5 * output.potentialCellWidth);
                indexPotentialMesh(&output,i,j,k)->pos = pos;
            }
        }
    }
    printf("Got Past Initialisation of posititions of meshcells and potential mesh cells\n");

    assignObjectsToMesh(objects,numObjects,&output);

    return output;




}

vec3 accelerationBetweenObjects(mesh* inputMesh, int A, int B){
    const double g = 6.67e-11;
    const double epsilon = inputMesh->epsilon;
    // vec3 outputAcc = inputMesh->objects[A].acc;
    vec3 outputAcc = vec3From(0.0,0.0,0.0);
    const vec3 iPos = inputMesh->objects[A].pos;
   

    const vec3 jPos = inputMesh->objects[B].pos;
    const double jMass = inputMesh->objects[B].mass;

    const vec3 seperation = sub_vec3(&jPos,&iPos); // Should be the other way round but more computationally efficient to encode the minus sign here
    const double acc_mul = (g * jMass) / sqrt(pow(vec3_mag_squared(seperation) + epsilon*epsilon,3.0));
    const vec3 acc_mul_sep = scalar_mul_vec3(acc_mul,&seperation);
    outputAcc = add_vec3(&outputAcc,&acc_mul_sep);
    // return outputAcc;
    return acc_mul_sep;
}

void getNeighbhors(int* inputArray, const int i, const int j, const int k, mesh* inputMesh){
    int counter = 0;
    const int numCells = inputMesh->numMeshCellsPerSideLength;
    
    for (int iLoop = -1; iLoop < 2; iLoop++){
         for (int jLoop = -1; jLoop < 2; jLoop++){
            for (int kLoop = -1; kLoop < 2; kLoop++){ // This gives the correct number of 
                
                // Logic goes here
                
                const int newI = i + iLoop;
                const int newJ = j + jLoop;
                const int newK = k + kLoop;
                // // Check that the new i,j,k's are valid
                
                const unsigned char iNotValid = (newI >= numCells) || (newI < 0);
                const unsigned char jNotValid = (newJ >= numCells) || (newJ < 0);
                const unsigned char kNotValid = (newK >= numCells) || (newK < 0);
                if (iNotValid || jNotValid || kNotValid){
                    inputArray[counter] = -1;
                    counter++;
                    continue;
                }

                
                meshCell* temp = indexMesh(inputMesh,newI,newJ,newK);
                const int head = temp->head;
 
                fflush(stdout);
                inputArray[counter] = head;
                
                // and ends before here
                counter++;

                
        
            }
        }   
    }
}

void accelerationFromCell(mesh* inputMesh,const unsigned int i, const unsigned int j,const unsigned int k ){
    meshCell* cell = indexMesh(inputMesh,i,j,k);
    int neighbhors[27]; // max neighbhors = 27 = 3**3 
    getNeighbhors(neighbhors,i,j,k,inputMesh);  
    const double a = 0.7*inputMesh->potentialCellWidth;
    const double aSquared = a*a;
    // Loop through every i'th object in cell
    for (int i = cell->head; i != -1; i = inputMesh->objects[i].next){
        vec3 accel = vec3From(0.0,0.0,0.0);
        vec3 posi = inputMesh->objects[i].pos;
        for (int j = 0; j < 27; j++){
            if (neighbhors[j] == -1){
                continue;
            }           
            // Calculate forces between i'th object and every neighbhoring cell
            for (int k = neighbhors[j]; k != -1; k = inputMesh->objects[k].next){
                vec3 posj = inputMesh->objects[k].pos;
                vec3 sep = sub_vec3(&posj,&posi);
                double magSquared = vec3_mag_squared(sep);
                if (magSquared > aSquared){
                    continue;
                }
                vec3 temp = accelerationBetweenObjects(inputMesh,i,k);
                accel = add_vec3(&accel,&temp);
                if (inputMesh->objects[k].next == -1){
                    break;
                }
            }

        }
        inputMesh->objects[i].acc = accel;
    }

 
}

// Calculate short range forces on each mesh using leapfrog method

void shortRangeForces(mesh* inputMesh){
    /* 
    1) Loop through every cell
    2) Check if each cell is on the edge or not and account for out of bound reads
    3) Otherwise Loop through every neighboring cell to take into account short range forces from them
    */

   const unsigned int maxCount = inputMesh->numMeshCellsPerSideLength;
   for (int i = 0; i < maxCount; i++){
        for (int j = 0; j < maxCount; j++){
            for (int k = 0; k < maxCount; k++){
                
                accelerationFromCell(inputMesh,i,j,k); // Does 2 and 3 lol

            }
        }
   }


}

void moveObject(const int iOld,const  int jOld,const int kOld,const int iNew,const int jNew,const int kNew,mesh* inputMesh, const int objectIndex){
    // printf("Attempting to move object %d from (%d,%d,%d) -> (%d,%d,%d)\n",objectIndex,iOld,jOld,kOld,iNew,jNew,kNew);
    static int numCalls = 0;
    numCalls++;
    object* objects = inputMesh->objects;
    objects[objectIndex].numTimesMoved++;
    meshCell* oldMeshCell = indexMesh(inputMesh,iOld,jOld,kOld);
    meshCell* newMeshCell = indexMesh(inputMesh,iNew,jNew,kNew);
    int prior = -1;

    // find object in old meshCell
    for (int i = oldMeshCell->head; i != -1; i = objects[i].next){
        if (i == objectIndex){
            
            // Remove from old list
            if (i == oldMeshCell->head){ // if I is the first entry in the list
                oldMeshCell->head = objects[objectIndex].next;
                // add to new list
                objects[objectIndex].next = newMeshCell->head;
                newMeshCell->head = objectIndex;
                
                return;
            }
            // Remove from old list
            objects[prior].next = objects[objectIndex].next;
            // Add to new list
            objects[objectIndex].next = newMeshCell->head;
            newMeshCell->head = objectIndex;
            return;

        }
        prior = i;
    }

    

 
    
    fprintf(stderr,"\nFailed to find object in meshcell\n Attempted to move object %d from (%d,%d,%d) -> (%d,%d,%d)\n Function has been called %d time(s).\n",objectIndex,iOld,jOld,kOld,iNew,jNew,kNew,numCalls);
    printf("Here are the objects contained in oldMeshCell:\n");
    printMeshCellObjects(oldMeshCell,inputMesh);
    printf("Here are the objects contained in newMeshCell:\n");
    printMeshCellObjects(newMeshCell,inputMesh);
    int actualI;
    int actualJ;
    int actualK;
    findObjectInMesh(inputMesh,objectIndex,&actualI,&actualJ,&actualK);
    printf("Object %d is located in mechCell (%d,%d,%d) and has been moved %d times\n",objectIndex,actualI,actualJ,actualK,objects[objectIndex].numTimesMoved);
    printf("There have been %d iterations\n",inputMesh->numIterationsPast);
    exit(1);


}

// Same a moveObject but for potential mesh
void moveObjectPotential(const int iOld,const  int jOld,const int kOld,const int iNew,const int jNew,const int kNew,mesh* inputMesh, const int objectIndex){
    // printf("Attempting to move object %d from (%d,%d,%d) PotentialMeshCell-> (%d,%d,%d)PotentialMeshCell\n",objectIndex,iOld,jOld,kOld,iNew,jNew,kNew);
    static int numCalls = 0;
    numCalls++;
    object* objects = inputMesh->objects;
    objects[objectIndex].numTimesMoved++;
    meshCell* oldMeshCell = indexPotentialMesh(inputMesh,iOld,jOld,kOld);
    meshCell* newMeshCell = indexPotentialMesh(inputMesh,iNew,jNew,kNew);
    int prior = -1;

    // find object in old meshCell
    for (int i = oldMeshCell->head; i != -1; i = objects[i].nextPotentialMesh){
        if (i == objectIndex){
            
            // Remove from old list
            if (i == oldMeshCell->head){ // if I is the first entry in the list
                oldMeshCell->head = objects[objectIndex].nextPotentialMesh;
                // add to new list
                objects[objectIndex].nextPotentialMesh = newMeshCell->head;
                newMeshCell->head = objectIndex;
                
                return;
            }
            // Remove from old list
            objects[prior].nextPotentialMesh = objects[objectIndex].nextPotentialMesh;
            // Add to new list
            objects[objectIndex].nextPotentialMesh = newMeshCell->head;
            newMeshCell->head = objectIndex;
            return;

        }
        prior = i;
    }

    fprintf(stderr,"\nFailed to find object in meshcell\n Attempted to move object %d from (%d,%d,%d) -> (%d,%d,%d)\n Function has been called %d time(s).\n",objectIndex,iOld,jOld,kOld,iNew,jNew,kNew,numCalls);
    printf("Here are the objects contained in oldMeshCell:\n");
    printMeshCellObjects(oldMeshCell,inputMesh);
    printf("Here are the objects contained in newMeshCell:\n");
    printMeshCellObjects(newMeshCell,inputMesh);
    int actualI;
    int actualJ;
    int actualK;
    findObjectInPotentialMesh(inputMesh,objectIndex,&actualI,&actualJ,&actualK); // Note: This will return the objects meshCell coordinates not nextMeshCell coordinates;
    printf("Object %d is located in mechCell (%d,%d,%d) and has been moved %d times\n",objectIndex,actualI,actualJ,actualK,objects[objectIndex].numTimesMoved);
    printf("There have been %d iterations\n",inputMesh->numIterationsPast);
    exit(1);


}

void periodicBoxBoundary(int obNum, mesh* inputMesh){
    object ob = inputMesh->objects[obNum];
    vec3 pos = ob.pos;
    vec3 vel = ob.vel;
    if (isnan(vel.x)){
        vel.x = 0.0;
    }
    if (isnan(vel.y)){
        vel.y = 0.0;
    }
    if (isnan(vel.z)){
        vel.z = 0.0;
    }
    double width = inputMesh->universeWidth;
    if (pos.x < 0.0 || (isnan(pos.x) && signbit(pos.x))){
        pos.x = width;
        vel.x = -1 * abs(vel.x);
        // vel.x = 0.0;
    }

   if (pos.x > width || isnan(pos.x)){
        pos.x = 0.0;
        vel.x =  abs(vel.x);
        // vel.x = 0.0;
        
    }

    if (pos.y < 0.0 || (isnan(pos.y) && signbit(pos.y))){
        pos.y = width;
        vel.y = -1 * abs(vel.y);
        // vel.y = 0.0;
    }

   if (pos.y > width || isnan(pos.y)){
        pos.y = 0;
        vel.y =  abs(vel.y);
        // vel.y = 0.0;
    }

    if (pos.z < 0.0 || (isnan(pos.x) && signbit(pos.x))){
        pos.z = width;
        vel.z = -1 * abs(vel.z);
        // vel.z = 0.0;
    }

   if (pos.z > width || isnan(pos.x)){
        pos.z = 0.0;
        vel.z = abs(vel.z);
        // vel.z = 0.0;
    }

    inputMesh->objects[obNum].pos = pos;
    inputMesh->objects[obNum].vel = vel;
    

    
}

void GreenKSpace(mesh* inputMesh){
    // Multiply every entry in rho_array by the greens function 1/k^2
    fftw_complex* rho_array = inputMesh->rho_array;
    int N = inputMesh->numPotentialMeshCellsPerSideLength;
    // Iterate over every rho_i
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                fftw_complex* rho_i = indexRhoArray(inputMesh,i,j,k);
                // Compute there x, y , and z positions


                vec3 kvec;
            
                kvec.x = 2.0*M_PI/inputMesh->numPotentialMeshCellsPerSideLength * i;
                kvec.y = 2.0*M_PI/inputMesh->numPotentialMeshCellsPerSideLength * j;
                kvec.z = 2.0*M_PI/inputMesh->numPotentialMeshCellsPerSideLength * k;
                double mag = vec3_mag_squared(kvec);
                if (mag != 0.0){
                *rho_i/=mag;
                }
                *rho_i /= (2*M_PI); // normalisation constant from fourier transform


            }
        }
    }
}

void longRangeForces(mesh* inputMesh){
    //TODO: This function
    // Calculate mass density in each potentialMeshCell and write it to the rho_array
    meshCell* potentialMeshCells = inputMesh->potentialMeshCells;
    int numPotentialMeshCells = inputMesh->numPotentialMeshCells;
    object* objects = inputMesh->objects;
    for (int i = 0; i < numPotentialMeshCells; i++){
        // Loop through every object in mesh Cell
        double totalMassInCell = 0.0;
        // printPotentialMeshCellObjects(&potentialMeshCells[i],inputMesh);
        for (int j = potentialMeshCells[i].head; j != -1; j = objects[j].nextPotentialMesh){
            totalMassInCell += objects[j].mass;

        }
        totalMassInCell /= potentialMeshCells[i].volume; // Calculate the mass density
        // Write to rho_array
        int iCoord = objects[potentialMeshCells[i].head].iPotential;
        int jCoord = objects[potentialMeshCells[i].head].jPotential;
        int kCoord = objects[potentialMeshCells[i].head].kPotential;
        // int iCoord;
        // int jCoord;
        // int kCoord;
        

        fftw_complex* ijkthDensity = indexRhoArray(inputMesh,iCoord,jCoord,kCoord); 
        *ijkthDensity = totalMassInCell + 0*I;
         
        
    }

    // Now to fourier transform
    fftw_execute(inputMesh->plan);
    // Multiply by Greens function in k space
    GreenKSpace(inputMesh);
    // multiply by ik
    int N = inputMesh->numPotentialMeshCellsPerSideLength;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                fftw_complex* rho_i = indexRhoArray(inputMesh,i,j,k);
                // Compute there x, y , and z positions


                vec3 ktemp;
            
                ktemp.x = 2.0*M_PI/inputMesh->numPotentialMeshCellsPerSideLength * i;
                ktemp.y = 2.0*M_PI/inputMesh->numPotentialMeshCellsPerSideLength * j;
                ktemp.z = 2.0*M_PI/inputMesh->numPotentialMeshCellsPerSideLength * k;
                const double magSquared = vec3_mag_squared(ktemp);
                const double mag = sqrt(magSquared); 
                fftw_complex ikx;
                fftw_complex iky;
                fftw_complex ikz;

                ikx = 0.0 + ktemp.x*I;
                iky = 0.0 + ktemp.y*I;
                ikz = 0.0 + ktemp.z*I;
                // Create vector of E field in k space at each point
                int index = indexArray(i,j,k,N); 
                // TODO: fill ikxArray, ikyArray, and ikzArray  with ik*G*rho
                const double g = 6.67e-11;
                const double constant = 4.0*M_PI*g; 
                inputMesh->kxArray[index] = ikx * (*rho_i) * constant; 
                inputMesh->kxArray[index] = iky * (*rho_i) * constant;
                inputMesh->kyArray[index] = ikz * (*rho_i) * constant;


            }
        }
    }

    fftw_execute(inputMesh->kxplan);
    fftw_execute(inputMesh->kyplan);
    fftw_execute(inputMesh->kzplan); // Now kx,ky,and kz contain the Ex, Ey, and Ez components of the field for each cell
    // Work out acceleration on each particle;
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                int index = indexArray(i,j,k,N);
                meshCell* potentialCell = &inputMesh->potentialMeshCells[index];
                // Calculate field strength at given point
                vec3 E;
                
                E.x = creal(inputMesh->kxArray[index]);
                E.y = creal(inputMesh->kyArray[index]);
                E.z = creal(inputMesh->kzArray[index]);
                //TODO: Interpolate Forces for better accuracy
                // update acceleration of each particle in cell
                for (int l = potentialCell->head; l != -1; l = inputMesh->objects[l].nextPotentialMesh){
                    object* ob = &inputMesh->objects[l];
                    // E is the acceleration, as F = mE = ma -> a = E
                    ob->acc = add_vec3(&ob->acc,&E);
              
                }
            }
        }
    }
}

void meshCellLeapFrogStep(mesh* inputMesh, double dt){
    shortRangeForces(inputMesh);
    static int numLongForcesComputes = 0;
    fprintf(stderr,"Comptuing long range forces for the %d'th time\n",numLongForcesComputes);
    numLongForcesComputes++;
    longRangeForces(inputMesh);
    for (int i = 0; i < inputMesh->numObjects; i++){
        vec3 temp = scalar_mul_vec3(dt,&inputMesh->objects[i].acc);
        inputMesh->objects[i].vel = add_vec3(&inputMesh->objects[i].vel,&temp);
        vec3 temp2 = scalar_mul_vec3(dt,&inputMesh->objects[i].vel);
        inputMesh->objects[i].pos = add_vec3(&inputMesh->objects[i].pos,&temp2);

        // Check if object has crossed into seperate cell
        object* objects = inputMesh->objects;

        int xIndex;
        int yIndex;
        int zIndex;

        int xIndexPotential;
        int yIndexPotential;
        int zIndexPotential;
        
        
        // Check if outside of universe box
        periodicBoxBoundary(i,inputMesh);
        assignObjectToMeshCell(i,inputMesh,&xIndex,&yIndex,&zIndex);
        assignObjectToPotentialMeshCell(i,inputMesh,&xIndexPotential,&yIndexPotential,&zIndexPotential);
        if (xIndex < 0){
            printf("xIndex = %d\n",xIndex);
        }
        if ((xIndex != objects[i].i) || (yIndex != objects[i].j) || (zIndex != objects[i].k)){
            
            moveObject(objects[i].i,objects[i].j,objects[i].k,xIndex,yIndex,zIndex,inputMesh,i);
            objects[i].i = xIndex;
            objects[i].j = yIndex;
            objects[i].k = zIndex;


        }
        if ((xIndexPotential != objects[i].iPotential) || (yIndexPotential != objects[i].jPotential) || (zIndexPotential != objects[i].kPotential)){

            moveObjectPotential(objects[i].iPotential,objects[i].jPotential,objects[i].kPotential,xIndexPotential,yIndexPotential,zIndexPotential,inputMesh,i);
            objects[i].iPotential = xIndexPotential;
            objects[i].jPotential = yIndexPotential; 
            objects[i].kPotential = zIndexPotential;
            

            
        }        


    }

    
}
















