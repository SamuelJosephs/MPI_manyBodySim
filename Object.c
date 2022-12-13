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
} object;

object obFrom(vec3 pos, vec3 vel, vec3 acc, double mass, double KE, double GPE){
    object output;
    output.pos = pos;
    output.vel = vel;
    output.acc = acc;
    output.mass = mass;
    output.KE = KE;
    output.GPE = GPE;
    return output;
}