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
    double fx; // force components
    double fy;
    double fz;

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
    double* Garray;

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

void printAllPotentialMeshCellObjects(mesh* inputMesh){
    int N = inputMesh->numPotentialMeshCells;
    for (int i = 0; i < N; i++){
        int head = inputMesh->potentialMeshCells[i].head;
        printf("potentialMeshCell %d contains: ",i);
        for (int j = head; j != -1; j = inputMesh->objects[j].nextPotentialMesh){
            printf("%d, ",j);
        }
        printf("\n");
    }
}

int indexArray(int i, int j, int k, int N){
    return ((N*N*i) + N*j + k); // row major order for fftw
}



void prependObjectToChainMeshCell(meshCell* cell, int obnum, mesh* inputMesh){
    // printf("obnum = %d, cell->head = %d\n",obnum,cell->head);
    inputMesh->objects[obnum].next = cell->head;
    cell->head = obnum;
    return;
}
void prependObjectToPotentialMeshCell(meshCell* cell, int obnum, mesh* inputMesh){
    inputMesh->objects[obnum].nextPotentialMesh = cell->head;
    cell->head = obnum;
    return;
}
void assignObjectsToMesh(object* objects, mesh* inputMesh){
    int numObjects = inputMesh->numObjects;
    // First do a pass initialising all chainMeshCells and potentialMeshCells heads to be -1

    for (int i = 0; i < inputMesh->numChainMeshCells; i++){
        inputMesh->chainMeshCells[i].head = -1;
    }

    for (int i = 0; i < inputMesh->numPotentialMeshCells; i++){
        inputMesh->potentialMeshCells[i].head = -1;
    }

    for (int i = 0; i < numObjects; i++){
        int iChainPos = (int)((objects[i].pos.x / inputMesh->chainCellWidth));
        int iPotentialPos =(int)(objects[i].pos.x / inputMesh->potentialCellWidth);
        int jChainPos = (int)(objects[i].pos.y / inputMesh->chainCellWidth);
        int jPotentialPos = (int)(objects[i].pos.y / inputMesh->potentialCellWidth);
        int kChainPos = (int)(objects[i].pos.z / inputMesh->chainCellWidth);
        int kPotentialPos = (int)(objects[i].pos.z / inputMesh->potentialCellWidth);
        
        // printf("iChainPos, jChainPos, kChainPos = %d,%d,%d\n",iChainPos,jChainPos,kChainPos);
        int chainIndex = indexArray(iChainPos,jChainPos,kChainPos,inputMesh->numChainMeshCellsPerSideLength);
        int potentialIndex = indexArray(iPotentialPos,jPotentialPos,kPotentialPos,inputMesh->numPotentialMeshCellsPerSideLength);

        meshCell* chainMeshCell = &inputMesh->chainMeshCells[chainIndex];
        meshCell* potentialMeshCell = &inputMesh->potentialMeshCells[potentialIndex];
        // printf("i = %d, numObjects = %d\n",i,numObjects);
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



mesh meshFrom(object* objects, int numObjects, double universeWidth, double G, int numChainMeshCellsPerSideLength, int numPotentialMeshCellsPerChainCell,double epsilon
){
    mesh output;
    output.epsilon = epsilon;
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
    output.Garray = malloc(output.numPotentialMeshCells * sizeof(double));
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

    assignObjectsToMesh(output.objects,&output); 
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

//Populates inputArray with the potential indexes for neighbhoring potentialMeshCells
void getPotentialNeighbhorPotentialCells(int* inputArray, const int i, const int j, const int k, mesh* inputMesh){
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
                inputArray[counter] = index;
                
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
                        rho += CICW(x,xp,y,yp,z,zp,H) * ob->mass;
                        
                    }
                    
                }

                rho /= cell->volume;
                // write rho to rho_array
                inputMesh->rho_array[meshCellIndex] = rho;
            }
        }
    }
}
void checkIfRhoArrayIsZero(mesh* inputMesh,int count){
    const int N = inputMesh->numPotentialMeshCells;
    unsigned char all_zero = 1;
    for (int i = 0; i < N; i++){
        if (cabs(inputMesh->rho_array[i]) > 0.0){
            all_zero = 0;
        }
    }
    if (all_zero){
        printf("Rho Array is zero on count: %d\n",count);
    }
}
double K(int i, int N){
    if (i == 0){
        return 2*M_PI;
    }
    return 2.0*M_PI*((double) i / (double)N);
}

void D(double ki,double kj,double kk,fftw_complex* di, fftw_complex* dj, fftw_complex* dk, double H,double alpha){
    *di = I*alpha*sin(ki*H)/H + I*(1-alpha)*sin(2.0*ki*H)/(2.0*H);
    *dj = I*alpha*sin(kj*H)/H + I*(1-alpha)*sin(2.0*kj*H)/(2.0*H);
    *dk = I*alpha*sin(kk*H)/H + I*(1-alpha)*sin(2.0*kk*H)/(2.0*H);
    return;
}


// NGP: p = 0, CIC: p = 1, TSC: p = 2
fftw_complex U(fftw_complex kx, fftw_complex ky, fftw_complex kz, double H, double p){
    double Hcubed = H*H*H;
    fftw_complex a = csin(kx*H/2.0)/(kx*H/2.0);
    fftw_complex b = csin(ky*H/2.0)/(ky*H/2.0);
    fftw_complex c = csin(kz*H/2.0)/(kz*H/2.0);
    return cpow(a*b*c,p+1.0);
    
}

complex double USquared(fftw_complex kx, fftw_complex ky, fftw_complex kz,double H){
    complex double a = 1.0 - cpow(csin(kx*H/2.0),2.0) + (2.0/15.0)*cpow(csin(kx*H/2.0),4.0);
    complex double b = 1.0 - cpow(csin(ky*H/2.0),2.0) + (2.0/15.0)*cpow(csin(ky*H/2.0),4.0);
    complex double c = 1.0 - cpow(csin(kz*H/2.0),2.0) + (2.0/15.0)*cpow(csin(kz*H/2.0),4.0);
    return a*b*c;

}

double SkSpace(double k, double a){
   return (12.0/pow((k*a/2.0),4.0))*(2.0-2.0*cos(k*a/2.0)-(k*a/2.0)*sin(k*a/2.0));
}

void RkSpace(double kx, double ky, double kz, double a,fftw_complex* Rx, fftw_complex* Ry, fftw_complex* Rz){
    double kSquared = kx*kx + ky*ky + kz*kz;
    // printf("k^2 = %f\n",kSquared);
    double k = sqrt(kSquared);
    double s = SkSpace(k,a);
    // printf("s = %f\n",s);
    double sSquared = s*s;
    double product = -sSquared / kSquared;
    *Rx = -I*kx*product;
    *Ry = -I*ky*product;
    *Rz = -I*kz*product; 
    return;
}

void GkSpace(double p, double H, double alpha,double a, double* outputArray, int outputArrayLength,mesh* inputMesh ){
    printf("Computing G(k)\n");
    // Compute D(k)
    fftw_complex dx;
    fftw_complex dy;
    fftw_complex dz;
    complex double sumProductX = 0.0 + 0.0*I;
    complex double sumProductY = 0.0 + 0.0*I;
    complex double sumProductZ = 0.0 + 0.0*I;
    complex double sumUSquared = 0.0 + 0.0*I;
    const int N = inputMesh->numPotentialMeshCellsPerSideLength;
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                const double kx = K(iLoop,N);
                const double ky = K(jLoop,N);
                const double kz = K(jLoop,N);
                if (isnan(kx) || isnan(ky) || isnan(kz)){
                    exit(2);
                }




                // Compute U^2(k_n)R(K)
                complex double rx, ry, rz;
                RkSpace(kx,ky,kz,a,&rx,&ry,&rz);
                // printf("rx,ry,rz = %f + %f I, %f + %f I, %f + %f I\n",creal(rx),cimag(rx),creal(ry),cimag(ry),creal(rz),cimag(rz));
                if (isnan(creal(rx)) || isnan(cimag(rx)) || isnan(creal(ry)) || isnan(cimag(ry))|| isnan(creal(rz)) || isnan(cimag(rz))){
                    printf("rx real: %f,rx imag: %f, ry real: %f, ry imag: %f, rz real: %f, rz imgag:%f\n",creal(rx),cimag(rx),creal(ry),cimag(ry),creal(rz),cimag(rz) );
                    printf("(kx,ky,kz) = (%f,%f,%f)\n",kx,ky,kz);
                    exit(3);
                }
                const complex double usquared = USquared(kx,ky,kz,H);
                 if (isnan(creal(usquared)) || isnan(cimag(usquared)) || isnan(creal(usquared)) || isnan(cimag(usquared))|| isnan(creal(usquared)) || isnan(cimag(usquared))){
                    exit(4);
                }              
                const complex double sumProdTempx = usquared * rx;
                const complex double sumProdTempy = usquared * ry; 
                const complex double sumProdTempz = usquared * rz;
                sumProductX += sumProdTempx;
                sumProductY += sumProdTempy;
                sumProductZ += sumProdTempz;
                sumUSquared += usquared;
                // printf("sumproduct x,y,z, uSquared = %f + %f I,%f + %f I,%f + %f I,%f + %f I\n",creal(sumProductX),cimag(sumProductX),creal(sumProductY),cimag(sumProductY),creal(sumProductZ),cimag(sumProductZ),creal(sumUSquared),cimag(sumUSquared));



            }
        }
    }
    // printf("Got past first loop\n");
    fflush(stdout);
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                const double kx = K(iLoop,N);
                const double ky = K(jLoop,N);
                const double kz = K(kLoop,N);
                complex double dx;
                complex double dy;
                complex double dz;
                D(kx,ky,kz,&dx,&dy,&dz,H,alpha);
                // printf("dx, dy, dz = %f + %f I, %f + %f I, %f + %f I\n",creal(dx),cimag(dx),creal(dy),cimag(dy),creal(dz),cimag(dz));

                if (isnan(creal(dx)) || isnan(cimag(dx)) || isnan(creal(dy)) || isnan(cimag(dy))|| isnan(creal(dz)) || isnan(cimag(dz))){
                    exit(5);
                }

                // Work out D(k) squared
                complex double dkSquared = dx*conj(dx) + dy*conj(dy) + dz*conj(dz);
                if (isnan(creal(dkSquared)) || isnan(cimag(dkSquared))){
                    exit(6);
                }
                const int index = indexArray(iLoop,jLoop,kLoop,N);
                dx *= sumProductX;
                dy *= sumProductY;
                dz *= sumProductZ;
                if (isnan(creal(dx)) || isnan(cimag(dx)) || isnan(creal(dy)) || isnan(cimag(dy))|| isnan(creal(dz)) || isnan(cimag(dz))){
                    exit(7);
                }

                complex double numerator = dx + dy + dz;
                complex double denominator = dkSquared * sumUSquared;
                
                complex double result = numerator / denominator;
                outputArray[index] = result;
                // printf("Result: %f + %f I\n",creal(result),cimag(result));
                // printf("numerator  = %f + %f I, demoninator = %f + %f I\n",creal(numerator),cimag(numerator),creal(denominator),cimag(denominator));
                if (isnan(creal(result)) || isnan(cimag(result))){
                    exit(8);
                }

            }
        }
    }




    return;
}

void scanRhoArrayForNan(mesh* inputMesh){
    int count = 0;
    for (int i = 0; i < inputMesh->numPotentialMeshCells; i++){
        if (isnan(creal(inputMesh->rho_array[i])) || isnan(cimag(inputMesh->rho_array[i]))){
            count++;
        }
    }
    // printf("Encountered %d nan values in rho_array\n",count);
}

void scanGarrayForNan(mesh* inputMesh){
    int count = 0;
    for (int i = 0; i < inputMesh->numPotentialMeshCells; i++){
        if (isnan(inputMesh->Garray[i])){
            count++;
        }
    }
    printf("Encountered %d nan values in Garray\n",count);
}

void longRangeForces(mesh* inputMesh,double dt){
    static unsigned int numTimesCalled = 0;
    numTimesCalled++;
    checkIfRhoArrayIsZero(inputMesh,numTimesCalled);
    fftw_execute(inputMesh->rhoPlanForward);
    scanRhoArrayForNan(inputMesh);
   // multiply every element in rho_array by it's greens function
   const int N = inputMesh->numPotentialMeshCellsPerSideLength; 
   for (int iLoop = 0; iLoop < N; iLoop++){
    for (int jLoop = 0; jLoop < N; jLoop++){
        for (int kLoop = 0; kLoop < N; kLoop++){
            const int index = indexArray(iLoop,jLoop,kLoop,N);
            if (inputMesh->Garray[index] == 0.0){
                continue;
            }
            // if (inputMesh->Garray[index] != 0.0){
            //     printf("G(%d,%d,%d) = %f + %f I")
            // }
            inputMesh->rho_array[index] *= inputMesh->Garray[index];
            }
        }
    }
    // printf("After multiplying by greens function:\n");
    // scanRhoArrayForNan(inputMesh);
    fftw_execute(inputMesh->rhoPlanBackward);
    // printf("Backwards pass:\n");
    // scanRhoArrayForNan(inputMesh);
    fftw_complex* const rho_array = inputMesh->rho_array;
    double H = inputMesh->potentialCellWidth;
   // compute force at each cell
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                const int index = indexArray(iLoop,jLoop,kLoop,N);
                const double fx = (rho_array[iLoop] - rho_array[iLoop - 1]) /H; // - grad phi 
                const double fy = (rho_array[jLoop] - rho_array[jLoop - 1]) /H; 
                const double fz = (rho_array[kLoop] - rho_array[kLoop - 1]) /H; 
                inputMesh->potentialMeshCells[index].fx = fx;
                inputMesh->potentialMeshCells[index].fy = fy;
                inputMesh->potentialMeshCells[index].fz = fz;
            }
        }
   }

    int neighbhors[27];
    for (int iLoop = 0; iLoop < N; iLoop++){
        for (int jLoop = 0; jLoop < N; jLoop++){
            for (int kLoop = 0; kLoop < N; kLoop++){
                getPotentialNeighbhorPotentialCells(neighbhors,iLoop,jLoop,kLoop,inputMesh);
                int index = indexArray(iLoop,jLoop,kLoop,N);
                for (int obnum = inputMesh->potentialMeshCells[index].head; obnum != -1; obnum = inputMesh->objects[obnum].nextPotentialMesh ){
                    vec3 obForce = vec3From(0.0,0.0,0.0);
                    const vec3 objectPos = inputMesh->objects[obnum].pos;
                    for (int n = 0; n < 27; n++){
                        if (neighbhors[n] == -1){
                            continue;
                        }
                        int head = inputMesh->potentialMeshCells[neighbhors[n]].head;
                        double meshcellfx = inputMesh->potentialMeshCells[neighbhors[n]].fx; 
                        double meshcellfy = inputMesh->potentialMeshCells[neighbhors[n]].fy; 
                        double meshcellfz = inputMesh->potentialMeshCells[neighbhors[n]].fz;

                        vec3 meshCellPos =  inputMesh->potentialMeshCells[neighbhors[n]].pos;
                        meshcellfx *= CICW(meshCellPos.x,objectPos.x,meshCellPos.y,objectPos.y,meshCellPos.z,objectPos.z,H);
                        meshcellfy *= CICW(meshCellPos.x,objectPos.x,meshCellPos.y,objectPos.y,meshCellPos.z,objectPos.z,H);
                        meshcellfz *= CICW(meshCellPos.x,objectPos.x,meshCellPos.y,objectPos.y,meshCellPos.z,objectPos.z,H);
 
                        obForce.x += meshcellfx;
                        obForce.y += meshcellfy;
                        obForce.z += meshcellfz;
                           
                        
                    }
                    // obForce = scalar_mul_vec3(1.0/inputMesh->objects[obnum].mass,&obForce);
                    inputMesh->objects[obnum].acc = add_vec3(&obForce,&inputMesh->objects[obnum].acc);
                    
                    printf("(fx,fy,fz) = %f,%f,%f | num times called = %d\n",obForce.x,obForce.y,obForce.z,numTimesCalled); 
    
                                                // Now do Euler step on object
                    // vec3 newVel = scalar_mul_vec3(dt,&obForce);
                    vec3 newVel = scalar_mul_vec3(dt,&inputMesh->objects[obnum].acc);
                    newVel = add_vec3(&newVel,&inputMesh->objects[obnum].vel);
                    vec3 newPos = scalar_mul_vec3(dt,&newVel);
                    newPos = add_vec3(&newPos,&inputMesh->objects[obnum].pos); 
                    if (isnan(newVel.x) || isnan(newVel.y) || isnan(newVel.z) || isnan(newVel.z) || isnan(newPos.x) || isnan(newPos.y) || isnan(newPos.z)){
                        printf("Encountered nan value at (%d,%d,%d)\n",iLoop,jLoop,kLoop);
                        continue;
                    }
                    // printf("Force %d = (%f,%f,%f)\n",index,obForce.x,obForce.y,obForce.z);
                    inputMesh->objects[obnum].vel = newVel;
                    
                    inputMesh->objects[obnum].pos = newPos;
                }

            }
        }
   }
  


}

void periodicBox(mesh* inputMesh,double offset){
    double universeWidth = inputMesh->universeWidth;
    
    for (int i = 0; i < inputMesh->numObjects; i++){
        vec3 pos = inputMesh->objects[i].pos;

        if (pos.x >= universeWidth){
            pos.x = 0.0 + offset;
        }
        if (pos.x < 0.0){
            pos.x = universeWidth - offset;
        }
        if (pos.y >= universeWidth){
            pos.y = 0.0 + offset;
        }
        if (pos.y < 0.0){
            pos.y = universeWidth - offset;
        }
        if (pos.z >= universeWidth){
            pos.z = 0.0 + offset;
        }
        if (pos.z < 0.0){
            pos.z = universeWidth - offset;
        }
        inputMesh->objects[i].pos = pos;
    }


}

vec3 accelerationBetweenObjects(mesh* inputMesh, int A, int B){
    // const double g = 6.67e-11;
    const double g = 1.0 / (4*M_PI);
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

void accelerationFromCell(mesh* inputMesh,const unsigned int i, const unsigned int j,const unsigned int k ){
    int index = indexArray(i,j,k,inputMesh->numChainMeshCellsPerSideLength);
    // meshCell* cell = indexArray(i,j,k,inputMesh->numPotentialMeshCellsPerSideLength);
    meshCell* cell = &inputMesh->chainMeshCells[index];
    int neighbhors[27]; // max neighbhors = 27 = 3**3 
    getPotentialNeighbhors(neighbhors,i,j,k,inputMesh);  
    const double a = 0.7*inputMesh->chainCellWidth;
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

void shortRangeForces(mesh* inputMesh){
    /* 
    1) Loop through every cell
    2) Check if each cell is on the edge or not and account for out of bound reads
    3) Otherwise Loop through every neighboring cell to take into account short range forces from them
    */

   const unsigned int maxCount = inputMesh->numChainMeshCellsPerSideLength;
   for (int i = 0; i < maxCount; i++){
        for (int j = 0; j < maxCount; j++){
            for (int k = 0; k < maxCount; k++){
                
                accelerationFromCell(inputMesh,i,j,k); // Does 2 and 3 lol

            }
        }
   }


}