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



void prependObjectToChainMeshCell(meshCell* cell, int obnum, mesh* inputMesh){
    inputMesh->objects[obnum].next = cell->head;
    cell->head = obnum;
    return;
}
void prependObjectToPotentialMeshCell(meshCell* cell, int obnum, mesh* inputMesh){
    inputMesh->objects[obnum].nextPotentialMesh = cell->head;
    cell->head = obnum;
    return;
}
void assignObjectsToMesh(object* objects,const unsigned int numObjects, mesh* inputMesh){

    // First do a pass initialising all chainMeshCells and potentialMeshCells heads to be -1

    for (int i = 0; i < inputMesh->numChainMeshCells; i++){
        inputMesh->chainMeshCells[i].head = -1;
    }

    for (int i = 0; i < inputMesh->numPotentialMeshCells; i++){
        inputMesh->potentialMeshCells[i].head = -1;
    }

    for (int i = 0; i < numObjects; i++){
        int iChainPos = (int)(objects[i].pos.x / inputMesh->chainCellWidth);
        int iPotentialPos =(int)(objects[i].pos.x / inputMesh->potentialCellWidth);
        int jChainPos = (int)(objects[i].pos.y / inputMesh->chainCellWidth);
        int jPotentialPos = (int)(objects[i].pos.y / inputMesh->potentialCellWidth);
        int kChainPos = (int)(objects[i].pos.z / inputMesh->chainCellWidth);
        int kPotentialPos = (int)(objects[i].pos.z / inputMesh->potentialCellWidth);

        int chainIndex = indexArray(iChainPos,jChainPos,kChainPos,inputMesh->numChainMeshCellsPerSideLength);
        int potentialIndex = indexArray(iPotentialPos,jPotentialPos,kPotentialPos,inputMesh->numPotentialMeshCellsPerSideLength);

        meshCell* chainMeshCell = &inputMesh->chainMeshCells[chainIndex];
        meshCell* potentialMeshCell = &inputMesh->potentialMeshCells[potentialIndex];

        prependObjectToChainMeshCell(chainMeshCell,i,inputMesh);
        prependObjectToPotentialMeshCell(potentialMeshCell,i,inputMesh);


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

void getPotentialNeighbhors(int* inputArray, const int i, const int j, const int k, mesh* inputMesh){
    int counter = 0;
    const int numCells = inputMesh->numPotentialMeshCellsPerSideLength;
    
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

                int index = indexArray(newI,newJ,newK,numCells);
                meshCell* temp = &inputMesh->potentialMeshCells[index];
                const int head = temp->head;
 
                fflush(stdout);
                inputArray[counter] = head;
                
                // and ends before here
                counter++;

                
        
            }
        }   
    }
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

void assignCharge(mesh* inputMesh){
    int neighbhors[27];
    int N = inputMesh->numPotentialMeshCellsPerSideLength;
    int H = inputMesh->potentialCellWidth;
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                getPotentialNeighbhors(neighbhors,iLoop,jLoop,kLoop,inputMesh);
                // Loop through every object in neighbhors, then convolve with window function to assign charge
                int meshCellIndex = indexArray(iLoop,jLoop,kLoop,N);
                meshCell* cell = &inputMesh->potentialMeshCells[meshCellIndex];
                double xp = cell->pos.x;
                double yp = cell->pos.y;
                double zp = cell->pos.z;
                double rho = 0.0;
                for (int neighbhor = 0; neighbhor < 27; neighbhor++){
                    int head = neighbhors[neighbhor];
                    if (head == -1){
                        continue;
                    }

                    // Now loop through objects in chain
                    for (int c = head; c != -1; c = inputMesh->objects[c].nextPotentialMesh){
                        object* ob = &inputMesh->objects[c];
                        double x = ob->pos.x;
                        double y = ob->pos.y;
                        double z = ob->pos.z;
                        rho += CICW(x,xp,y,yp,z,zp,H);
                        
                    }
                    
                }

                rho /= cell->volume;
                // write rho to rho_array
                inputMesh->rho_array[meshCellIndex] = rho;
            }
        }
    }
}

double K(int i, int N){
    return 2.0*M_PI*((double) i / (double)N);
}

void D(double ki,double kj,double kk,fftw_complex* di, fftw_complex* dj, fftw_complex* dk, double H,double alpha){
    *di = I*alpha*sin(ki*H)/H + I*(1-alpha)*sin(2.0*ki*H)/(2.0*H);
    *dj = I*alpha*sin(kj*H)/H + I*(1-alpha)*sin(2.0*kj*H)/(2.0*H);
    *dk = I*alpha*sin(kk*H)/H + I*(1-alpha)*sin(2.0*kk*H)/(2.0*H);
}


// NGP: p = 0, CIC: p = 1, TSC: p = 2
fftw_complex U(fftw_complex kx, fftw_complex ky, fftw_complex kz, double H, double p){
    double Hcubed = H*H*H;
    fftw_complex a = csin(kx*H/2.0)/(kx*H/2.0);
    fftw_complex b = csin(ky*H/2.0)/(ky*H/2.0);
    fftw_complex c = csin(kz*H/2.0)/(kz*H/2.0);
    return cpow(a*b*c,p+1.0);
    
}

fftw_complex USquared(fftw_complex kx, fftw_complex ky, fftw_complex kz,double H){
    fftw_complex a = 1.0 - cpow(csin(kx*H/2.0),2.0) + (2.0/15.0)*cpow(csin(kx*H/2.0),4.0);
    fftw_complex b = 1.0 - cpow(csin(ky*H/2.0),2.0) + (2.0/15.0)*cpow(csin(ky*H/2.0),4.0);
    fftw_complex c = 1.0 - cpow(csin(kz*H/2.0),2.0) + (2.0/15.0)*cpow(csin(kz*H/2.0),4.0);
    return a*b*c;

}

double SkSpace(double k, double a){
   return (12.0/pow((k*a/2.0),4.0))*(2.0-2.0*cos(k*a/2.0)-(k*a/2.0)*sin(k*a/2.0))
}

void RkSpace(double kx, double ky, double kz, double a,fftw_complex* Rx, fftw_complex* Ry, fftw_complex* Rz){
    double kSquared = kx*kx + ky*ky + kz*kz;
    double k = sqrt(kSquared);
    double s = SkSpace(k,a);
    double sSquared = s*s;
    double product = -sSquared / kSquared;
    *Rx = -I*kx*product;
    *Ry = -I*ky*product;
    *Rz = -I*kz*product; 
}

void GkSpace(double p, double H, double alpha, fftw_complex* outputArray, int outputArrayLength ){
    for (int i = 0; i < outputArrayLength; i++){
        for (int j = 0; j < outputArrayLength; j++){
            for (int k = 0; k < outputArrayLength; k++){
                double kx = K(i,outputArrayLength);
                double ky = K(j,outputArrayLength);
                double kz = K(k,outputArrayLength);
                fftw_complex dx;
                fftw_complex dy;
                fftw_complex dz;
                D(kx,ky,kz,&dx,&dy,&dz,H,alpha);
                //TODO: compute proper greens function and write to outputArray, it will be used in all future calculations and only needs to be computed once.
            }
        }
    }
}



void longRangeForces(mesh* inputMesh){
   fftw_execute(inputMesh->rhoPlanForward);
    
}

