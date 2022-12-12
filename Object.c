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