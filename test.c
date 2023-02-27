#include <stdio.h>
#include "mesh.c"
#include "systems.c"
const static double epsilon = 1e3; 
const static double dt = 1e-2;
const double universeWidth = 20.0 * 778e9;


int main(){
    int numObjects = 200;
    object* distribution = uniformRandom(numObjects,epsilon,dt,universeWidth);
    leapFrogSetup(distribution,numObjects,epsilon,dt);
    mesh domainMesh = meshFrom(universeWidth/2.0,universeWidth,2,distribution,numObjects,epsilon);
    FILE* outputFile = fopen("Test_Data.csv","w+");
    printf("Num mesh cells per side length from main = %d\n",domainMesh.numMeshCellsPerSideLength);
    int counter = 0;
    printf("Num Mesh Cells Per side length: %d\n",domainMesh.numMeshCellsPerSideLength);
    for (int i = 0; i < domainMesh.numMeshCells; i++){
        printf("Accsessing the %d'th meshCell out of %d\n",i,domainMesh.numMeshCells);
        printMeshCellObjects(&domainMesh.meshCells[i],&domainMesh);
        
    }
    

    for (int i = 0; i < 300; i++){
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
    fflush(stdout);
    

}