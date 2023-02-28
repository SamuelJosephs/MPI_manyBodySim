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

object* sphericalDistribution(int* NUM_OBJECTS,double epsilon, double dt){
	*NUM_OBJECTS = 100;
	srand(1);
	int count = 0;
	object* outArray = malloc(*NUM_OBJECTS * sizeof(object));
	const double a = 10 * RAD_JUPITER;
	while (count < *NUM_OBJECTS) {
		double xrand = (((float) rand())/(float) RAND_MAX) * a;
		double yrand = (((float) rand())/(float) RAND_MAX) * a;
		double zrand = (((float) rand())/(float) RAND_MAX) * a;

		if (xrand*xrand + yrand*yrand + zrand*zrand > a*a)
			continue;
		
		vec3 position = vec3From(xrand,yrand,zrand);
		vec3 vel = vec3From(0.,0.,0.);
		vec3 acc = vec3From(0.,0.,0.);
		double mass = MASS_Earth;
		outArray[count] = obFrom(position,vel,acc,mass,0.,0.);
		count++;
	}
	printf("Count = %d\n",count);
	return outArray;
}

object* uniformRandom(const int NUM_OBJECTS,double epsilon, double dt,const double universeWidth){
	object* objects = malloc(NUM_OBJECTS * sizeof(object));
	srand(1);
	const vec3 zeroVec = vec3From(0.0,0.0,0.0);
	for (int i = 0; i < NUM_OBJECTS; i++){
		object temp;
		vec3 pos;
		temp.acc = zeroVec;
		temp.GPE = 0.0;
		temp.KE = 0.0;
		temp.mass = 10*MASS_Earth;
		temp.next = -1;
		temp.nextPotentialMesh = -1;

		pos.x = (((float) rand())/(float) RAND_MAX) * universeWidth;
		pos.y = (((float) rand())/(float) RAND_MAX) * universeWidth;
		pos.z = (((float) rand())/(float) RAND_MAX) * universeWidth;

		temp.pos = pos;
		temp.vel = zeroVec;
		temp.numTimesMoved = 0;
		objects[i] = temp;	
	}
	return objects;
}


object* uniform(const int numObjectsPerSideLength,double epsilon, double dt,const double universeWidth){
	object* output = malloc(numObjectsPerSideLength * numObjectsPerSideLength * numObjectsPerSideLength * sizeof(object));
	if (output == NULL){
		fprintf(stderr, "Failed to allocate objects\n");
		exit(1);
	}
	vec3 zeroVec = vec3From(0.0,0.0,0.0);
	for (int i = 0; i < numObjectsPerSideLength; i++){
		for (int j = 0; j < numObjectsPerSideLength; j++){
			for (int k = 0; k < numObjectsPerSideLength; k++){
				// Work out position of object
				vec3 pos;
				pos.x = (double) i * (universeWidth / (double) numObjectsPerSideLength);
				pos.y = (double) j * (universeWidth / (double) numObjectsPerSideLength);
				pos.z = (double) k * (universeWidth / (double) numObjectsPerSideLength);
				object ob;
				ob.pos = pos;
				ob.vel = zeroVec;
				ob.acc = zeroVec;
				ob.next = -1;
				ob.nextPotentialMesh = -1;
				ob.mass = MASS_Earth;
				ob.KE = 0.0;
				ob.GPE = 0.0;
				int index = (numObjectsPerSideLength*numObjectsPerSideLength*i + numObjectsPerSideLength*j + k );
				output[index] = ob;

			}
		}
	}

	return output;
}