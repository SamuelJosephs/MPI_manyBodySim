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
// TODO: Check if object moves to a different cell and impose periodic boundary's 
// I.E.: IF an object leaves the domain of the "universe" it appears in the opposite side
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
    int numIterationsPast;
} mesh;



meshCell* indexMesh(const mesh* inputMesh, int i, int j, int k){
    const int N = inputMesh->numMeshCellsPerSideLength;
    return &(inputMesh->meshCells[((k * N*N) + (j*N) + i)]) ;
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
        // printf("(i,j,k) = (%d,%d,%d)\n",i,j,k);
        return;
    }

    meshCell* cell = indexMesh(inputMesh,i,j,k);
    // Prepend to linked list in Cell
    inputMesh->objects[obNum].next = cell->head;
    cell->head = obNum;
}


void printMeshCellObjects(meshCell* inputCell,mesh* inputMesh){
    
    for (int i = inputCell->head; i != -1; i = inputMesh->objects[i].next){
        printf("%i Index:(%d,%d,%d), ",i,inputMesh->objects[i].i,inputMesh->objects[i].j,inputMesh->objects[i].k);
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


    }
}

mesh meshFrom(const double meshCellWidth, const double universeWidth, object* objects, int numObjects,double epsilon){
    mesh output;
    output.numIterationsPast = 0;
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
    // Loop through every i'th object in cell
    for (int i = cell->head; i != -1; i = inputMesh->objects[i].next){
        vec3 accel = vec3From(0.0,0.0,0.0);
        for (int j = 0; j < 27; j++){
            if (neighbhors[j] == -1){
                continue;
            }
            // Calculate forces between i'th object and every neighbhoring cell
            for (int k = neighbhors[j]; k != -1; k = inputMesh->objects[k].next){
 
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
    static int numCalls = 0;
    numCalls++;
    object* objects = inputMesh->objects;
    objects[objectIndex].numTimesMoved++;
    meshCell* oldMeshCell = indexMesh(inputMesh,iOld,jOld,kOld);
    meshCell* newMeshCell = indexMesh(inputMesh,iNew,jNew,kNew);
    int prior = -1;

    objects[objectIndex].i = iNew;
    objects[objectIndex].j = jNew;
    objects[objectIndex].k = kNew;
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

void periodicBoxBoundary(int obNum, mesh* inputMesh){
    object ob = inputMesh->objects[obNum];
    vec3 pos = ob.pos;
    vec3 vel = ob.vel;
    double width = inputMesh->universeWidth;
    unsigned char xCheck = (pos.x < 0) || (pos.x > width);
    unsigned char yCheck = (pos.y < 0) || (pos.y > width);
    unsigned char zCheck = (pos.z < 0) || (pos.z > width);
    unsigned char changed = 0;
    if (pos.x < 0.0){
        pos.x = width;
        vel.x = -1 * vel.x;
    }

   if (pos.x > width){
        pos.x = 0.0;
        vel.x = -1 * vel.x;
    }

    if (pos.y < 0.0){
        pos.y = width;
        vel.y = -1 * vel.y;
    }

   if (pos.y > width){
        pos.y = 0;
        vel.y = -1 * vel.y;
    }

    if (pos.z < 0.0){
        pos.z = width;
        vel.z = -1 * vel.z;
    }

   if (pos.z > width){
        pos.z = 0.0;
        vel.z = -1 * vel.z;
    }

    inputMesh->objects[obNum].pos = pos;
    inputMesh->objects[obNum].vel = vel;
    

    
}

void meshCellLeapFrogStep(mesh* inputMesh, double dt){
    shortRangeForces(inputMesh);

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
        
        
        // Check if outside of universe box
        periodicBoxBoundary(i,inputMesh);
        assignObjectToMeshCell(i,inputMesh,&xIndex,&yIndex,&zIndex);
        if ((xIndex != objects[i].i) || (yIndex != objects[i].j) || (zIndex != objects[i].k)){
            
            moveObject(objects[i].i,objects[i].j,objects[i].k,xIndex,yIndex,zIndex,inputMesh,i);
            

        }
        


    }

    
}
















