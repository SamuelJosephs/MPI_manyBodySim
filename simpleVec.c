typedef struct {
	double x;
	double y;
	double z;
} vec3;

inline vec3 vec3From(const double x, const double y,const double z){
    vec3 output;
    output.x = x;
    output.y = y;
    output.z = z;
    return output;
}

vec3 add_vec3(const vec3* restrict const left,const vec3* restrict const right){
	vec3 output;
	output.x = left->x + right->x;
	output.y = left->y + right->y;
	output.z = left->z + right->z;
    return output;
}

vec3 sub_vec3(const vec3* restrict const left,const vec3* restrict const right){
	vec3 output;
	output.x = left->x - right->x;
	output.y = left->y - right->y;
	output.z = left->z - right->z;
    return output;
}

vec3 scalar_mul_vec3(const double mul,const vec3* restrict const right){
	vec3 output;
	output.x = mul * right->x;
	output.y = mul * right->y;
	output.z = mul * right->z; 
    return output;
}

double vec3_mag_squared(const vec3 input){
    return (input.x*input.x + input.y*input.y + input.z*input.z);
}


