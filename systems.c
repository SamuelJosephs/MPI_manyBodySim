#ifndef LEAPFROG
#define LEAPFROG
#include "leapFrog.c"
#endif

#include <string.h>


#define MASS_Earth 5.9722e24
#define MASS_SUN 1.989e30
#define MASS_JUPITER 1.898e27

#define RAD_EARTH 1.4959787e11
#define RAD_JUPITER 778e9
#define G 6.67e-11
#define V_RAD(R) sqrt((G*MASS_SUN)/R)

object* EarthSunJupiter(int* NUM_OBJECTS,double epsilon, double dt){
    *NUM_OBJECTS = 3;
    object* ARRAY = malloc(3*sizeof(object));
    if (ARRAY == NULL){
        fprintf(stderr,"Failed to allocate memory\n");
        exit(1);
    }

	// Set up system
	vec3 earthPos = vec3From(RAD_EARTH,0.0,0.0);
	vec3 earthVel = vec3From(0.0,V_RAD(RAD_EARTH),0.0);
	vec3 earthAcc = vec3From(0.0,0.0,0.0);
	object earth = obFrom(earthPos,earthVel,earthAcc,MASS_Earth,0.0,0.0);

	vec3 jupiterPos = vec3From(-RAD_JUPITER,0.0,0.0);
	vec3 jupiterVel = vec3From(0.0,-V_RAD(RAD_JUPITER),0.0);
	vec3 jupiterAcc = vec3From(0.0,0.0,0.0);
	object jupiter = obFrom(jupiterPos,jupiterVel,jupiterAcc,MASS_JUPITER,0.0,0.0);

	vec3 sunPos = vec3From(0.0,0.0,0.0);
	vec3 sunVel = vec3From(0.0,0.0,0.0);
	vec3 sunAcc = vec3From(0.0,0.0,0.0);
	object sun = obFrom(sunPos,sunVel,sunAcc,MASS_SUN,0.0,0.0);

	memcpy(&ARRAY[0],&earth,sizeof(object));
    memcpy(&ARRAY[1],&sun,sizeof(object));
    memcpy(&ARRAY[2],&jupiter,sizeof(object));

    return ARRAY;
	
}