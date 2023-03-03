#include <stdio.h>
#include "mesh2.c"
#include "systems.c"
const static double epsilon = 1e-3; 
const static double dt = 1.0e-5;
const double universeWidth = 200.0;
const int numObjectsPerSideLength = 5;
const int numObjects = numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength; 
double g = 6.67e-11;
int main(int argc, char** argv){
    object* objects = uniform(numObjectsPerSideLength,epsilon,dt,universeWidth);
    mesh test = meshFrom(objects,numObjects,universeWidth,g,10,15);
    double H = test.potentialCellWidth;
    // assignCharge(&test);
    GkSpace(1.0,H,0.5,0.7*test.chainCellWidth,test.Garray,test.numPotentialMeshCellsPerSideLength,&test);
    
    scanGarrayForNan(&test);
    FILE* output = fopen("Test2_data.csv","w+");
    for (int i = 0; i < 200; i++){
        periodicBox(&test,0.5*test.potentialCellWidth);
        assignObjectsToMesh(objects,&test);
        // printAllPotentialMeshCellObjects(&test);
        assignCharge(&test);
        shortRangeForces(&test);
        longRangeForces(&test,dt);
        
        writeObjectsToFile(output,&test);
    }
    printf("Succsess!");
}