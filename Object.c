#ifndef SIMPLE_VEC
#define SIMPLE_VEC
#include <simpleVec.c>
#endif


typedef struct {
    vec3 acc;
    vec3 vel;
    vec3 pos;
    double mass;
    double KE;
    double GPE;
    int next; // next object in head, -1 means you are at the ned of the list
    int nextPotentialMesh;
    int i;
    int j;
    int k;
    int iPotential;
    int jPotential;
    int kPotential;
    int numTimesMoved;
} object;

object obFrom(vec3 pos, vec3 vel, vec3 acc, double mass, double KE, double GPE){
    object output;
    output.pos = pos;
    output.vel = vel;
    output.acc = acc;
    output.mass = mass;
    output.KE = KE;
    output.GPE = GPE;
    output.next = 0;
    return output;
}