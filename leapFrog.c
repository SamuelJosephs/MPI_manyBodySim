#ifndef SIMPLE_VEC
#define SIMPLE_VEC
#include "simpleVec.c"
#endif

#ifndef OBJECT
#define OBJECT
#include "Object.c"
#endif

#include <math.h>

vec3 acceleration(const object* restrict const inputArray,const unsigned int inputArraySize, const unsigned int ithObject, const double epsilon){
    vec3 outputAcc = vec3From(0.0,0.0,0.0);
    const vec3 iPos = inputArray[ithObject].pos;
    
    for (unsigned int j = 0; j < inputArraySize; j++){
        if (j == ithObject){
            continue;
        }
        const vec3 jPos = inputArray[j].pos;
        const double jMass = inputArray[j].mass;

        const vec3 seperation = sub_vec3(&jPos,&iPos); // Should be the other way round but more computationally efficient to encode the minus sign here
        const double acc_mul = jMass * sqrt(pow(vec3_mag_squared(seperation) + epsilon*epsilon,3.0));
        const vec3 acc_mul_sep = scalar_mul_vec3(acc_mul,&seperation);
        outputAcc = add_vec3(&outputAcc,&acc_mul_sep);
    }

    return outputAcc;
}

// updates posArray
void leapFrogStep( object* restrict const inputArray, const unsigned int inputArraySize, const unsigned int procRank, const unsigned int worldSize, const double epsilon, const double dt){
    // First Calculate which elements of acceleration we are going to work on
    const unsigned int batchSize = inputArraySize / worldSize;
    const unsigned int startIndex = batchSize * procRank;
    
    for (int i = startIndex; i < (startIndex + batchSize); i++){ // possible segfault if startIndex + batchSize > inputArraySize
        vec3 acc = acceleration(inputArray,inputArraySize,i,epsilon);
        const vec3 acc_dt = scalar_mul_vec3(dt,&acc);
        inputArray[i].vel = add_vec3(&inputArray[i].vel,&acc_dt);
    }
    // Now update positions

    for (int i = startIndex; i < (startIndex + batchSize); i++){
        const vec3 v_dt = scalar_mul_vec3(dt,&inputArray[i].vel);
        inputArray[i].pos = add_vec3(&inputArray[i].pos,&v_dt);
    }

}