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



typedef struct meshCell {
    int head; // head of chain
    vec3 pos; // position of center of cell;
    double rho;
    double volume;

} meshCell;



typedef struct mesh {
    unsigned int numChainMeshCells;
    unsigned int numChainMeshCellsPerSideLength;
    unsigned int numPotentialMeshCellsPerSideLength;
    unsigned int numPotentialMeshCells;
    meshCell* chainMeshCells;
    meshCell* potentialMeshCells;
    double chainCellWidth;
    double potentialCellWidth;
    double universeWidth;
    double G;
    object* objects;
    unsigned int numObjects;
    double epsilon;
    fftw_complex* rho_array;
    fftw_plan rhoPlanForward;
    fftw_plan rhoPlanBackward;
    fftw_complex* kxArray;
    fftw_plan kxplan;
    fftw_complex* kyArray;
    fftw_plan kyplan;
    fftw_complex* kzArray;
    fftw_plan kzplan;
} mesh;



int indexArray(int i, int j, int k, int N){
    return ((N*N*i) + N*j + k); // row major order for fftw
}

void assignObjectToMeshCell(int obNum,mesh* inputMesh,int* iOut, int* jOut, int* kOut){
    // Work out fraction across x,y, and z axis
    object ob = inputMesh->objects[obNum];
    double width = inputMesh->universeWidth;
    int numCells = inputMesh->numChainMeshCellsPerSideLength;
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
    int index = indexArray(i,j,k,inputMesh->numChainMeshCellsPerSideLength);
    meshCell* cell = &inputMesh->chainMeshCells[index];
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
    int index = indexArray(i,j,k,inputMesh->numPotentialMeshCellsPerSideLength);
    meshCell* cell = &inputMesh->potentialMeshCells[index];
    // Prepend to linked list in Cell
    fprintf(stderr, "Potential: Cell Head = %d, next: %d\n",cell->head,inputMesh->objects[cell->head].nextPotentialMesh);
    inputMesh->objects[obNum].nextPotentialMesh = cell->head;
    cell->head = obNum;
    if (cell->head == inputMesh->objects[obNum].nextPotentialMesh){
        fprintf(stderr,"Potential: Cell head == obnum.nextPotentialMesh, %d,%d This will cause an infinite loop\n",cell->head,inputMesh->objects[obNum].nextPotentialMesh);
        // exit(1);
    }
}

void calcIndex(double x, double y, double z, int numPerSide, int* i, int* j, int* k, mesh* inputMesh){
    double width = inputMesh->universeWidth;
    double fracX = x / width;
    double fracY = y / width;
    double fracZ = z / width;
    int iTemp = (int)(fracX * (double)numPerSide -0.5);
    int jTemp = (int)(fracY * (double)numPerSide - 0.5);
    int kTemp = (int)(fracZ * (double)numPerSide - 0.5);
    *i = iTemp;
    *j = jTemp;
    *k = kTemp;    
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

void assignChainMeshCellsPositions(mesh* inputMesh){
    int N = inputMesh->numChainMeshCellsPerSideLength;
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                double xPos = iLoop * (inputMesh->universeWidth / inputMesh->numChainMeshCellsPerSideLength) - 0.5 * inputMesh->chainCellWidth;
                double yPos = jLoop * (inputMesh->universeWidth / inputMesh->numChainMeshCellsPerSideLength) - 0.5 * inputMesh->chainCellWidth;
                double zPos = kLoop * (inputMesh->universeWidth / inputMesh->numChainMeshCellsPerSideLength) - 0.5 * inputMesh->chainCellWidth;
                int index = indexArray(iLoop,jLoop,kLoop,N);
                inputMesh->chainMeshCells[index].pos.x = xPos;
                inputMesh->chainMeshCells[index].pos.y = yPos;
                inputMesh->chainMeshCells[index].pos.z = zPos; 
            }
        }
    }
    return;
}
void assignPotentialMeshCellsPositions(mesh* inputMesh){
    int N = inputMesh->numPotentialMeshCellsPerSideLength;
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                double xPos = iLoop * (inputMesh->universeWidth / inputMesh->numPotentialMeshCellsPerSideLength) - 0.5 * inputMesh->potentialCellWidth;
                double yPos = jLoop * (inputMesh->universeWidth / inputMesh->numPotentialMeshCellsPerSideLength) - 0.5 * inputMesh->potentialCellWidth;
                double zPos = kLoop * (inputMesh->universeWidth / inputMesh->numPotentialMeshCellsPerSideLength) - 0.5 * inputMesh->potentialCellWidth;
                int index = indexArray(iLoop,jLoop,kLoop,N);
                inputMesh->potentialMeshCells[index].pos.x = xPos;
                inputMesh->potentialMeshCells[index].pos.y = yPos;
                inputMesh->potentialMeshCells[index].pos.z = zPos; 
            }
        }
    }
    return;
}



mesh meshFrom(object* objects, int numObjects, double universeWidth, double G, int numChainMeshCellsPerSideLength, int numPotentialMeshCellsPerChainCell){
    mesh output;
    output.objects = objects;
    output.numObjects = numObjects;
    output.numChainMeshCellsPerSideLength = numChainMeshCellsPerSideLength;
    output.numPotentialMeshCellsPerSideLength = output.numChainMeshCellsPerSideLength * numPotentialMeshCellsPerChainCell;
    output.universeWidth = universeWidth;
    output.chainCellWidth = universeWidth / (double) numChainMeshCellsPerSideLength;
    output.potentialCellWidth = output.chainCellWidth / numPotentialMeshCellsPerChainCell;
    output.G = G;

    output.numChainMeshCells = output.numChainMeshCellsPerSideLength *  output.numChainMeshCellsPerSideLength *  output.numChainMeshCellsPerSideLength;
    output.numPotentialMeshCells = output.numPotentialMeshCellsPerSideLength * output.numPotentialMeshCellsPerSideLength * output.numPotentialMeshCellsPerSideLength; 

    // allocate potentialMeshCells and ChainMeshCells
    output.chainMeshCells = malloc(output.numChainMeshCells * sizeof(meshCell));
    output.potentialMeshCells = malloc(output.numPotentialMeshCells * sizeof(meshCell));
    double a = output.chainCellWidth * output.chainCellWidth * output.chainCellWidth; 
    for (int i = 0; i < output.numChainMeshCells; i++){
        output.chainMeshCells[i].head = -1;
        output.chainMeshCells[i].volume = a; 
        
        
    }
    a = output.potentialCellWidth * output.potentialCellWidth * output.potentialCellWidth; 
    for (int i = 0; i < output.numPotentialMeshCells; i++){
        output.potentialMeshCells[i].head = -1;
        output.potentialMeshCells[i].volume = a;         
    }

    // Allocate fftw arrays
    int N = output.numPotentialMeshCellsPerSideLength;
    output.rho_array = fftw_malloc(output.numPotentialMeshCells * sizeof(complex));
    printf("Generating fftw3 plans\n");
    output.rhoPlanForward = fftw_plan_dft_3d(N,N,N,output.rho_array,output.rho_array,FFTW_FORWARD,FFTW_MEASURE);
    output.rhoPlanBackward = fftw_plan_dft_3d(N,N,N,output.rho_array,output.rho_array,FFTW_BACKWARD,FFTW_MEASURE);

    assignObjectsToMesh(output.objects,output.numObjects,&output); 
    assignChainMeshCellsPositions(&output);
    assignPotentialMeshCellsPositions(&output);
 
    return output; 

}


double CICw(double x,double x_p,double H){
    double difference = abs(x-x_p);
    if (difference > H){
        return 0.0;
    }
    return 1.0 - (difference / H);

}

double CICW(double x, double xp, double y, double yp, double z, double zp, double H){
    return CICw(x,xp,H) * CICw(y,yp,H) * CICw(z,zp,H);
}


