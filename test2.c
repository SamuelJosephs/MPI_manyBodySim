#include <stdio.h>
#include "mesh2.c"
#include "systems.c"
const static double epsilon = 1e3; 
const static double dt = 1e-6;
const double universeWidth = 200.0 * 778e9;
const int numObjectsPerSideLength = 10;
const int numObjects = numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength; 
double g = 6.67e-11;
int main(int argc, char** argv){
    object* objects = uniform(10,epsilon,dt,universeWidth);
    mesh test = meshFrom(objects,numObjects,universeWidth,g,10,5);
    printf("Succsess!");
}