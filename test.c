#include <stdio.h>
#include "mesh.c"
#include "systems.c"
const static double epsilon = 1e3; 
const static double dt = 1e-2;
const double universeWidth = 20.0 * 778e9;


int main(){
    int numObjects = 100;
    object* distribution = uniformRandom(numObjects,epsilon,dt,universeWidth);
    leapFrogSetup(distribution,numObjects,epsilon,dt);
    mesh domainMesh = meshFrom(universeWidth/20.0,universeWidth,distribution,numObjects,epsilon);
    FILE* outputFile = fopen("Test_Data.csv","w+");
    printf("Num mesh cells per side length from main = %d\n",domainMesh.numMeshCellsPerSideLength);
    int counter = 0;
    for (int i = 0; i < domainMesh.numMeshCells; i++){
        printMeshCellObjects(&domainMesh.meshCells[i],&domainMesh);
        
    }
    
    printMeshCellObjects(indexMesh(&domainMesh,0,4,0),&domainMesh);
    printMeshCellObjects(indexMesh(&domainMesh,2,4,1),&domainMesh);
    for (int i = 0; i < 100; i++){
        meshCellLeapFrogStep(&domainMesh,dt);


        if (counter == 10){
            printf("numobjects: %d\n",numObjects);
            writeOutputToFile(outputFile,distribution,numObjects);
            counter = 0;
            domainMesh.numIterationsPast++;
        }
        counter++;
        
    }
    // writeOutputToFile(outputFile,distribution,numObjects);
    printf("Succsess!\n");
    free(domainMesh.meshCells);
    free(distribution);
    

}