#include <stdio.h>
#include "mesh2.c"
#include "systems.c"
const static double epsilon = 1e1; 
const static double dt = 1.0e-2;
const double universeWidth = 2000.0;
const int numObjectsPerSideLength = 5;
const int numObjects = numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength; 
double g = 6.67e-11;
int main(int argc, char** argv){
    int nsteps;
    // getnsteps
    char* ptr;
    if (argc == 2){
        nsteps = strtol(argv[1],&ptr,10);
    } else {
        nsteps = 150;
    }
    object* objects = uniform(numObjectsPerSideLength,epsilon,dt,universeWidth);
    mesh test = meshFrom(objects,numObjects,universeWidth,g,5,15,epsilon);
    double H = test.potentialCellWidth;
    // assignCharge(&test);
    GkSpace(1.0,H,0.5,0.7*test.chainCellWidth,test.Garray,test.numPotentialMeshCellsPerSideLength,&test);
    
    scanGarrayForNan(&test);
    FILE* output = fopen("Test2_data.csv","w+");
    int counter = 0;
    for (int i = 0; i < nsteps; i++){
        periodicBox(&test,0.5*test.potentialCellWidth);
        assignObjectsToMesh(objects,&test);
        // printAllPotentialMeshCellObjects(&test);
        assignCharge(&test);
        shortRangeForces(&test);
        longRangeForces(&test,dt);
        if (counter == 10){
            writeObjectsToFile(output,&test);
            counter =  0;
            continue;
        }
        counter++;
    }
    printf("Succsess!");
}