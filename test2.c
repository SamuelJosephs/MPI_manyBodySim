#include <stdio.h>
#include "mesh2.c"
#include "systems.c"
const static double epsilon = 1e1; 
const static double dt = 1.0e0;
const double universeWidth = 2000.0;
const int numObjectsPerSideLength = 10;
const int numObjects = numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength; 
double g = 6.67e-11;
int main(int argc, char** argv){
    FILE* output = fopen("Test2_data.csv","w+");

    int nsteps;
    // getnsteps
    char* ptr;
    if (argc == 2){
        nsteps = strtol(argv[1],&ptr,10);
    } else {
        nsteps = 150;
    }
    object* objects = uniform(numObjectsPerSideLength,epsilon,dt,universeWidth);

    mesh test = meshFrom(objects,numObjects,universeWidth,g,5,12,epsilon);
    
    double H = test.potentialCellWidth;
    // assignCharge(&test);
    GkSpace(2.0,H,4.0/3.0,0.7*test.chainCellWidth,test.Garray,test.numPotentialMeshCellsPerSideLength,&test);
    
    scanGarrayForNan(&test);
    int counter = 0;
    // writeObjectsToFile(output,&test);
 // set all x components of velocity to 0
    for (int i = 0; i < test.numObjects; i++){
        test.objects[i].vel.x = 0.0;
    }
    writeObjectsToFile(output,&test);
    for (int i = 0; i < nsteps; i++){
        periodicBox(&test,0.001*test.potentialCellWidth);
        assignObjectsToMesh(objects,&test);
        // printAllPotentialMeshCellObjects(&test);
        assignCharge(&test);
        shortRangeForces(&test);
        longRangeForces(&test,dt);
        if (counter == 2){
            writeObjectsToFile(output,&test);
            counter =  0;
            continue;
        }
        counter++;
        periodicBox(&test,0.001*test.potentialCellWidth);
    }
    printf("Succsess!");
}