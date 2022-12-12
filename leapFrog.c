#ifndef SIMPLE_VEC
#define SIMPLE_VEC
#include "simpleVec.c"
#endif

#ifndef OBJECT
#define OBJECT
#include "Object.c"
#endif

#include <math.h>

vec3 acceleration_and_GPE(object* restrict const inputArray,const unsigned int inputArraySize, const unsigned int ithObject, const double epsilon){
    vec3 outputAcc = vec3From(0.0,0.0,0.0);
    const vec3 iPos = inputArray[ithObject].pos;
    const double iMass = inputArray[ithObject].mass;
    double GPE = 0;
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

        const double ith_GPE = (iMass * jMass) / sqrt(vec3_mag_squared(seperation));
        inputArray[ithObject].GPE += ith_GPE;

    }

    return outputAcc;
}

// updates pos, vel, KE, and GPE of every object in inputArray
void leapFrogStep( object* restrict const inputArray, const unsigned int inputArraySize, const unsigned int procRank, const unsigned int worldSize, const double epsilon, const double dt){
    unsigned int ws;
    unsigned int batchSize;
    unsigned int startIndex;
    if (worldSize > inputArraySize){ // edge case
        if (procRank >= inputArraySize){
            return;
        }
        ws = inputArraySize;
    } else {
        ws = worldSize;
    }

    if ((inputArraySize % worldSize != 0 )) {
        if (procRank == worldSize - 1){
            batchSize = worldSize % inputArraySize;
            startIndex = inputArraySize - batchSize;
        } else {
            batchSize = inputArraySize / ws;
            startIndex = batchSize * procRank;
        }
    } else {
            batchSize = inputArraySize / ws;
            startIndex = batchSize * procRank;
    }
    // First Calculate which elements of acceleration we are going to work on
    
   

    
    for (int i = startIndex; i < (startIndex + batchSize); i++){ // possible segfault if startIndex + batchSize > inputArraySize
        vec3 acc = acceleration_and_GPE(inputArray,inputArraySize,i,epsilon);
        const vec3 acc_dt = scalar_mul_vec3(dt,&acc);
        const vec3 acc_half_dt = scalar_mul_vec3(0.5*dt,&acc);
        const vec3 velSyncedToPos = add_vec3(&inputArray[i].vel,&acc_half_dt); // to calculate KE
        inputArray[i].vel = add_vec3(&inputArray[i].vel,&acc_dt);

        inputArray[i].KE = 0.5 * inputArray[i].mass * vec3_mag_squared(velSyncedToPos);

        // Now compute potential energy

        
    }
    // Now update positions

    for (int i = startIndex; i < (startIndex + batchSize); i++){
        const vec3 v_dt = scalar_mul_vec3(dt,&inputArray[i].vel);
        inputArray[i].pos = add_vec3(&inputArray[i].pos,&v_dt);
    }

}