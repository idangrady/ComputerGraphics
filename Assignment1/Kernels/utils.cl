#include "template/common.h"
#define SCRWIDTH	1280
#define SCRHEIGHT	720

typedef struct __attribute__((aligned(64))) 
{
	// Access all float4 as float3 by using var_name.xyz
	float4 O;	// Ray origin
	float4 D;	// Ray direction (normalized?)
	float4 rD_t;	// Reciprocal ray direction, and time distance of hit
	uint pixel;	// Pixel corresponding to the ray
	uint primIdx; // Type indedx (2 bit), instance index (10 bit) and primitive index (20 bit), if needed
	float2 uv;	// barycentric uv coordinates of hit
} Ray;

typedef struct __attribute__((aligned(64))) 
{
	float4 camPos;
	float4 topLeft;
	float4 topRight;
	float4 bottomLeft;
} Camera;

typedef struct __attribute__((aligned(64)))
{
	float4 vertex0; //Stores C0 as 4th float
	float4 vertex1; //Stores C1 as 4th float
	float4 vertex2; //Stores C2 as 4th float
	float N0, N1, N2;
	uint id;
} Triangle;

typedef struct __attribute__((aligned(8))){
	int matId;
	int textureID; // -1 if no texture
} TriEx;

typedef struct __attribute__((aligned(64)))
{
	float4 albedoSpecularity; // Albedo in first 3, specularity last float
	float4 absorption; // If higher than zero, then the material is glass
	bool isEmissive; // If this is a light
	uint medium;
} Material;

inline uint GetObjectType(uint idx){
    const uint mask = ~0 << 30;
    return (idx & mask) >> 30;
}

inline uint GetObjectIndex(uint idx){
    const uint mask = ~((~0 << 30) | 0xFFFFF);
    return (idx & mask) >> 20;
}

inline uint GetTriangleIndex(uint idx){
    const uint mask = ~(~0 << 20);
    return idx & mask;
}

inline uint makeId(uint type, uint id, uint tri){
    return (type << 30) + (id << 20) + (tri);
}

// Marsaglia xor32, hope this is fine
inline uint randomUint(uint* seed)
{
	*seed ^= *seed << 13;
	*seed ^= *seed >> 17;
	*seed ^= *seed << 5;
	return *seed;
}

inline ulong randomULong(ulong* seed){
	*seed ^= *seed << 13;
	*seed ^= *seed >> 17;
	*seed ^= *seed << 5;
	return *seed;
}

inline float randomFloat(uint* seed)
{
	uint r = randomUint(seed);
	return r * 2.3283064365387e-10f;
	//ret = 2.f * ret - 1.f;
	//printf("random float: %f\n", ret)
}

inline float randomUFloat(ulong* seed){
	ulong r = randomULong(seed);
	return r * 5.42101086242752217e-20f;
}

float3 GetDiffuseReflection(float3 N, ulong* seed){
	float3 random_vec;
	//printf("Entering while loop\n");
	bool gen = false;
	int i = 0;
	while(!gen){
		i++;
		float x = 2.f * randomUFloat(seed) - 1.f;
		float y = 2.f * randomUFloat(seed) - 1.f;
		float z = 2.f * randomUFloat(seed) - 1.f;
		//printf("Seed 3 = %u\n",seed);
		//printf("x = %f, y = %f, z = %f\n", x, y, z);

		random_vec = (float3)(x, y, z);
		if(length(random_vec) < 1.f){
			random_vec = normalize(random_vec);
			if(dot(N, random_vec) < 0) random_vec *= -1.f;
			gen = true;;
		}
		if(i > 100){
			printf("100 iterations? Vec = %f, %f, %f\n Seed = %u\n", x, y, z, *seed);
		}
	}
	//printf("Exited diffuse reflection while loop\n");
	return random_vec;
}

float3 CosineWeightedDiffuseReflection(ulong* seed){
	float r0 = randomFloat(seed);
	float r1 = randomFloat(seed);
	float r = sqrt(r0);
	float theta = 2 * PI * r1;
	float x = r * cos(theta);
	float y = r * sin(theta);
	return (float3)(x, y, sqrt(1 - r0));
}

Ray reflectRay(Ray ray, float3 I, float3 N){
	float3 dir = ray.D.xyz;
	float3 reflected = normalize(dir - 2.0f * (dot(dir, N) * N));
	Ray newRay = {
		(float4)(I + (0.0002f * reflected), 1),
		(float4)(reflected, 1),
		(float4)(1.0f / reflected, 1e34f),
		ray.pixel,
		0,
		(float2)(0,0)
	};
	return newRay;
}