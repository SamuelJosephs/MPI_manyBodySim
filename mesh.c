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
    meshCell*** meshCells;
    double cellWidth;
    double universeWidth;
    object* objects;
    unsigned int numObjects;
} mesh;

void assignObjectsToMesh(object* objects,const unsigned int numObjects, mesh* inputMesh){
    for (int i = 0; i < numObjects; i++){
        double fracX = objects[i].pos.x / inputMesh -> universeWidth;
        double fracY = objects[i].pos.y / inputMesh -> universeWidth;
        double fracZ = objects[i].pos.z / inputMesh -> universeWidth;

        int xIndex = fracX * inputMesh -> numMeshCellsPerSideLength;
        int yIndex = fracY * inputMesh -> numMeshCellsPerSideLength;
        int zIndex = fracZ * inputMesh -> numMeshCellsPerSideLength;

        meshCell* cell =  &(inputMesh -> meshCells[xIndex][yIndex][zIndex]);
        if (cell -> head == -1){
            objects[i].next = -1; // appending to linked list
            cell -> head = i;
        } else {
            objects[i].next = cell->head;
            cell->head =  i;
        }
    }
}

mesh meshFrom(const double meshCellWidth, const double universeWidth, object* objects, int numObjects){
    mesh output;
    output.objects = objects;
    output.numObjects = numObjects;
    output.universeWidth = universeWidth;
    output.numMeshCellsPerSideLength = (unsigned int) (universeWidth / meshCellWidth);
    output.numMeshCells = output.numMeshCellsPerSideLength * output.numMeshCellsPerSideLength * output.numMeshCellsPerSideLength; // cube it 
    output.cellWidth = meshCellWidth;
    // Create array of meshcells
    output.meshCells = malloc(output.numMeshCells * sizeof(meshCell));
    if (output.meshCells == NULL){
        fprintf(stderr,"Failed to allocate meshCells\n");
        exit(1);
    }
    meshCell temp;
    for (int i = 0; i < output.numMeshCells; i++){ // initialise next pointers to each meshcell to be -1 (empty)
        
        ((meshCell*)output.meshCells)[i] = temp; // Maybe undefined behaviour. Initialise heads to be empty (-1).
        ((meshCell*)output.meshCells)[i].head = -1;
        
    }

    for (int i = 0; i < output.numMeshCellsPerSideLength; i ++) {
        for (int j = 0; j < output.numMeshCellsPerSideLength; j++){
            for (int k = 0; k < output.numMeshCellsPerSideLength; k ++){
                vec3 pos;
                pos.x = (((double) i) / ((double) output.numMeshCellsPerSideLength)) * output.universeWidth - (0.5 * meshCellWidth); 
                pos.y = (((double) j) / ((double) output.numMeshCellsPerSideLength)) * output.universeWidth - (0.5 * meshCellWidth);
                pos.z = (((double) k) / ((double) output.numMeshCellsPerSideLength)) * output.universeWidth - (0.5 * meshCellWidth);
                output.meshCells[i][j][k].pos = pos;
            }
        }
    }

    assignObjectsToMesh(objects,numObjects,&output);




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
                accelerationFromCell(inputMesh,i,j,k);

            }
        }
   }
}

vec3 accelerationFromCell(mesh* inputMesh,const unsigned int i, const unsigned int j,const unsigned int k ){
    meshCell* cell = &inputMesh->meshCells[i][j][k];
    int neighbhors[27]; // max neighbhors = 27 = 3**3 
    neighbhors[26] = cell->head;
    getAdjacentCells(&neighbhors,inputMesh,i,j,k); // 0 - 25 are potentially adjacent cells, or 
    // Loop through every object in Cell
    for (int i = cell->head; i != -1; i = inputMesh->objects[i].next){
        // Loop through every object in this or neighbhoring cells
        vec3 accelForIthObject = vec3From(0.0,0.0,0.0);
        for (int j = 0; j < 27; j ++){
            int temp = neighbhors[j];
            if (temp == -1)
                continue;
            for (int k = temp, k != -1; k = inputMesh->objects[k].next){
                if (i == k)
                    continue;
                vec3 acc = accelerationBetweenObjects(inputMesh,i,k);
                accelForIthObject = add_vec3(&accelForIthObject,&acc); 
            }
        }
        
    }
}

void getAdjacentCells(int* neighbhorArray, mesh* inputMesh, int i, int j, int k){
    unsigned int counter = 0;
    for (int iLoop = i-1; iLoop <= i + 1; iLoop++){
        for (int jLoop = j-1; jLoop <= j + 1; jLoop++){
            for (int kLoop = k-1; iLoop <= k + 1; kLoop++){
                if ((iLoop == i) && (jLoop == j) && (kLoop == k))
                    continue;
                // Test to see if on the edges
                if (((iLoop >= inputMesh->numMeshCellsPerSideLength) || iLoop < 0) || ((jLoop >= inputMesh->numMeshCellsPerSideLength) || jLoop < 0) || ((kLoop >= inputMesh->numMeshCellsPerSideLength) || kLoop < 0)){
                    neighbhorArray[counter] = -1;
                    counter++;
                    continue;
                } else { // else write head to neighbhor array, this else isn't needed but I think it makes the code easier to follow.
                    neighbhorArray[counter] = inputMesh->meshCells[iLoop][jLoop][kLoop].head;
                    counter++;
                    continue;
                }

            }
        }   
    }
}

vec3 accelerationBetweenObjects(mesh* inputMesh, int A, int B){
    const double g = 6.67e-11;
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










