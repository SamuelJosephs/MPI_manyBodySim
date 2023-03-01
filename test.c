#include <stdio.h>
#include "mesh.c"
#include "systems.c"
const static double epsilon = 1e3; 
const static double dt = 1e-2;
const double universeWidth = 20.0 * 778e9;


int main(){
    int numObjectsPerSideLength = 3;
    int numObjects = numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength;
    object* distribution = uniform(numObjectsPerSideLength,epsilon,dt,universeWidth);
    leapFrogSetup(distribution,numObjects,epsilon,dt);
    mesh domainMesh = meshFrom(universeWidth/5,universeWidth,10,distribution,numObjects,epsilon);
    FILE* outputFile = fopen("Test_Data.csv","w+");
    printf("Num mesh cells per side length from main = %d\n",domainMesh.numMeshCellsPerSideLength);
    int counter = 0;
    printf("Num Mesh Cells Per side length: %d\n",domainMesh.numMeshCellsPerSideLength);
    for (int i = 0; i < domainMesh.numPotentialMeshCells; i++){
        // printf("Accsessing the %d'th meshCell out of %d\n",i,domainMesh.numPotentialMeshCells);
        // printPotentialMeshCellObjects(&domainMesh.potentialMeshCells[i],&domainMesh);
        if (domainMesh.potentialMeshCells[i].head != -1){
            printf("Accsessing the %d'th meshCell out of %d\n",i,domainMesh.numPotentialMeshCells);
            printPotentialMeshCellObjects(&domainMesh.potentialMeshCells[i],&domainMesh);            
        }       
    }
    

    for (int i = 0; i < 100; i++){
        meshCellLeapFrogStep(&domainMesh,dt);


        if (counter == 10){
            printf("numobjects: %d\n",numObjects);
            writeOutpuDensityProfileToFile(&domainMesh,outputFile);
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

    // for (int i = 0; i < domainMesh.numPotentialMeshCells; i++){
    //     // printf("Accsessing the %d'th meshCell out of %d\n",i,domainMesh.numPotentialMeshCells);
    //     // printPotentialMeshCellObjects(&domainMesh.potentialMeshCells[i],&domainMesh);
    //     if (domainMesh.potentialMeshCells[i].head != -1){
    //         printf("Accsessing the %d'th meshCell out of %d\n",i,domainMesh.numPotentialMeshCells);
    //         printPotentialMeshCellObjects(&domainMesh.potentialMeshCells[i],&domainMesh);            
    //     }       
    // }
        
    

}