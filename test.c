#include <stdio.h>
#include "mesh.c"
#include "systems.c"
const static double epsilon = 1e3; 
const static double dt = 1e-9;
const double universeWidth = 200.0 * 778e9;


int main(int argc, char** argv){
    int nsteps;
    // getnsteps
    char* ptr;
    if (argc == 2){
        nsteps = strtol(argv[1],&ptr,10);
    } else {
        nsteps = 150;
    }
    int numObjectsPerSideLength = 5;
    int numObjects = numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength;
    object* distribution = uniform(numObjectsPerSideLength,epsilon,dt,universeWidth);
    leapFrogSetup(distribution,numObjects,epsilon,dt);
    mesh domainMesh = meshFrom(universeWidth/5,universeWidth,50,distribution,numObjects,epsilon);
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
    

    for (int i = 0; i < nsteps; i++){
        meshCellLeapFrogStep(&domainMesh,dt);


        // if (counter == 10){
            printf("numobjects: %d\n",numObjects);
            // writeOutpuDensityProfileToFile(&domainMesh,outputFile);
            writeObjectsToFile(outputFile,&domainMesh);
            counter = 0;
            domainMesh.numIterationsPast++;
        // }
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