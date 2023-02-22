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

typedef struct meshCell {
    int head; // head of chain
    vec3 pos; // position of center of cell;
    double rho;

} meshCell;

typedef struct mesh {
    unsigned int numMeshCells;
    unsigned int numMeshCellsPerSideLength;
    meshCell* meshCells;
    double cellWidth;
    double universeWidth;
    object* objects;
    unsigned int numObjects;
    double epsilon;
} mesh;

meshCell* indexMesh(const mesh* inputMesh, int i, int j, int k){
    const int N = inputMesh->numMeshCellsPerSideLength;
    return &(inputMesh->meshCells[((k * N*N) + (j*N) + i)]) ;
}




void assignObjectsToMesh(object* objects,const unsigned int numObjects, mesh* inputMesh){
    for (int i = 0; i < numObjects; i++){
        double fracX = objects[i].pos.x / inputMesh -> universeWidth;
        double fracY = objects[i].pos.y / inputMesh -> universeWidth;
        double fracZ = objects[i].pos.z / inputMesh -> universeWidth;

        int xIndex = fracX * inputMesh -> numMeshCellsPerSideLength;
        int yIndex = fracY * inputMesh -> numMeshCellsPerSideLength;
        int zIndex = fracZ * inputMesh -> numMeshCellsPerSideLength;

        meshCell* cell =  indexMesh(inputMesh,xIndex,yIndex,zIndex);
        if (cell -> head == -1){
            objects[i].next = -1; // appending to linked list
            cell -> head = i;
        } else {
            objects[i].next = cell->head;
            cell->head =  i;
        }
    }
}

mesh meshFrom(const double meshCellWidth, const double universeWidth, object* objects, int numObjects,double epsilon){
    mesh output;
    output.objects = objects;
    output.numObjects = numObjects;
    output.universeWidth = universeWidth;
    output.numMeshCellsPerSideLength = (unsigned int) (universeWidth / meshCellWidth);
    // printf("There are %d meshCells per side length\n",output.numMeshCellsPerSideLength);
    output.numMeshCells = output.numMeshCellsPerSideLength * output.numMeshCellsPerSideLength * output.numMeshCellsPerSideLength; // cube it 
    output.cellWidth = meshCellWidth;
    output.epsilon = epsilon;
    // printf("Num mesh cells per side length = %d\n",output.numMeshCellsPerSideLength);
    // Create array of meshcells
    output.meshCells = malloc(output.numMeshCells * sizeof(meshCell));
    if (output.meshCells == NULL){
        fprintf(stderr,"Failed to allocate meshCells\n");
        exit(1);
    }
    meshCell temp;
    for (int i = 0; i < output.numMeshCells; i++){ // initialise next pointers to each meshcell to be -1 (empty)
        
        output.meshCells[i] = temp; // Maybe undefined behaviour. Initialise heads to be empty (-1).
        output.meshCells[i].head = -1;
        
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
    printf("Got Past Initialisation of posititions of meshcells\n");

    assignObjectsToMesh(objects,numObjects,&output);

    return output;




}

vec3 accelerationBetweenObjects(mesh* inputMesh, int A, int B){
    const double g = 6.67e-11;
    const double epsilon = inputMesh->epsilon;
    vec3 outputAcc = inputMesh->objects[A].acc;
    const vec3 iPos = inputMesh->objects[A].pos;
    const double iMass = inputMesh->objects[A].mass;

    const vec3 jPos = inputMesh->objects[B].pos;
    const double jMass = inputMesh->objects[B].mass;

    const vec3 seperation = sub_vec3(&jPos,&iPos); // Should be the other way round but more computationally efficient to encode the minus sign here
    const double acc_mul = (g * jMass) / sqrt(pow(vec3_mag_squared(seperation) + epsilon*epsilon,3.0));
    const vec3 acc_mul_sep = scalar_mul_vec3(acc_mul,&seperation);
    outputAcc = add_vec3(&outputAcc,&acc_mul_sep);
    return outputAcc;
}

void getNeighbhors(int* inputArray, const int i, const int j, const int k, mesh* inputMesh){
    int counter = 0;
    const int numCells = inputMesh->numMeshCellsPerSideLength;
    // printf("Num mesh cells per side length = %d\n",numCells);
    for (int iLoop = -1; iLoop < 2; iLoop++){
         for (int jLoop = -1; jLoop < 2; jLoop++){
            for (int kLoop = -1; kLoop < 2; kLoop++){ // This gives the correct number of 
                
                // Logic goes here
                //TODO: Do this function
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
                // printf("Counter: %d | i,j,k = %d,%d,%d\n",counter,newI,newJ,newK);
                // printf(" | i , j , k = %d, %d, %d\n",newI,newJ,newK);
                
                meshCell* temp = indexMesh(inputMesh,newI,newJ,newK);
                const int head = temp->head;
                // printf("Read i,j,k = %d,%d,%d\n",newI,newJ,newK);
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
    // Loop through every object in cell
    // Calculate acceleration between object in the cell and every object in the cell plus it's neighbhors
    for (int iLoop = inputMesh->objects[indexMesh(inputMesh,i,j,k)->head]; iLoop != -1; iLoop = inputMesh->objects[iLoop].next){
        vec3 accel = vec3From(0.0,0.0,0.0);
        for (int jLoop = 0; jLoop < 27; jLoop++){
            int kLoop = neighbhors[jLoop];
            if (kLoop == -1){ //TODO: add cutof distance
                continue;
            }
            vec3 temp = accelerationBetweenObjects(inputMesh,iLoop,jLoop);
            accel = add_vec3(&accel,&temp);
        }
        inputMesh->objects[iLoop].acc = accel;

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














