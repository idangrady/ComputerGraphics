#ifndef UTILS_IMPORT
#define UTILS_IMPORT
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
	uint primType; // 0 = sky, 1 = triangle, 2 = sphere
} Triangle;

typedef struct __attribute__((aligned(64))){
	float4 N;
	float2 uv0, uv1, uv2;
	int matId;
	int textureID; // -1 if no texture
	float A;
} TriEx;

typedef struct __attribute__((aligned(64)))
{
	float4 albedoSpecularity; // Albedo in first 3, specularity last float
	float4 absorption; // If higher than zero, then the material is glass
	bool isEmissive; // If this is a light
	uint medium;
} Material;

typedef struct __attribute__((aligned(8)))
{
	int width;
	int height;
} TextureData;

typedef struct __attribute__((aligned(32)))
{
	union{
		float4 aabbMin;
		struct{float minx, miny, minz; uint leftFirst;};
	};
	union{
		float4 aabbMax;
		struct{float maxx, maxy, maxz; uint triCount;};
	};
	//float4 aabbMin, aabbMax;			// boundary
	//uint leftFirst, triCount;			// count and start
} BVHNode; 

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
	while(!gen){
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
	}
	//printf("Exited diffuse reflection while loop\n");
	return random_vec;
}

float3 CosineWeightedDiffuseReflection(ulong* seed){
	float r0 = randomUFloat(seed);
	float r1 = randomUFloat(seed);
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

float4 getSkyBox(float3 dir, int width, int height, __constant float4* skybox) {
	//return (float4)(0.3, 0.3, 0.3, 1);
	float phi = atan2(dir.x, dir.z);
	float theta = atan2(hypot(dir.x, dir.z), dir.y);
	if (phi < 0) phi += 2.0f * PI;
	if (theta < 0) printf("Theta < 0?\n");
	float x = phi / (2.0f * PI);
	float y = theta / PI;
	int x_ = (int)(x * width);
	int y_ = (int)(y * height);
	// float4 pixel = skybox[y_ * width + x_];
	// printf("R: %f\n", pixel.x);
	return skybox[y_ * width + x_];
}

float IntersectAABB( Ray* ray, BVHNode* node )
{
    float tx1 = (node->minx - ray->O.x) * ray->rD_t.x, tx2 = (node->maxx - ray->O.x) * ray->rD_t.x;
    float tmin = min( tx1, tx2 ), tmax = max( tx1, tx2 );
    float ty1 = (node->miny - ray->O.y) * ray->rD_t.y, ty2 = (node->maxy - ray->O.y) * ray->rD_t.y;
    tmin = max( tmin, min( ty1, ty2 ) ), tmax = min( tmax, max( ty1, ty2 ) );
    float tz1 = (node->minz - ray->O.z) * ray->rD_t.z, tz2 = (node->maxz - ray->O.z) * ray->rD_t.z;
    tmin = max( tmin, min( tz1, tz2 ) ), tmax = min( tmax, max( tz1, tz2 ) );
    if (tmax >= tmin && tmin < ray->rD_t.w && tmax > 0) return tmin; else return 1e30f;
}

void IntersectTri(Ray* ray, __global Triangle* tri, int id) 
{
    const float3 edge1 = tri->vertex1.xyz - tri->vertex0.xyz;
    const float3 edge2 = tri->vertex2.xyz - tri->vertex0.xyz;
    const float3 h = cross(ray->D.xyz, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O.xyz - tri->vertex0.xyz;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D.xyz, q);
    if (v < 0 || u + v > 1) return;
    const float t = f * dot(edge2, q);
    if (t > 0.0001f && t < ray->rD_t.w){
        ray->rD_t.w = t;
        ray->uv = (float2)(u, v);
        ray->primIdx = makeId(1, 0, id);
    }
}

int DoesIntersectTri(Ray* ray, __global Triangle* tri){
    const float3 edge1 = tri->vertex1.xyz - tri->vertex0.xyz;
    const float3 edge2 = tri->vertex2.xyz - tri->vertex0.xyz;
    const float3 h = cross(ray->D.xyz, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return 0; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O.xyz - tri->vertex0.xyz;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return 0;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D.xyz, q);
    if (v < 0 || u + v > 1) return 0;
    const float t = f * dot(edge2, q);
    if(t > 0.0001f) return 1;
    else return 0;
}

int ShadowRayIntersect(Ray* ray, __global Triangle* tri){
    const float3 edge1 = tri->vertex1.xyz - tri->vertex0.xyz;
    const float3 edge2 = tri->vertex2.xyz - tri->vertex0.xyz;
    const float3 h = cross(ray->D.xyz, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return 0; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O.xyz - tri->vertex0.xyz;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return 0;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D.xyz, q);
    if (v < 0 || u + v > 1) return 0;
    const float t = f * dot(edge2, q);
    if(t > 0.0001f && t < ray->rD_t.w) return 1;
    else return 0;
}

uint RandomLight(uint lightCount, ulong* seed){
	return randomULong(seed) % lightCount;
}

float3 RandomPointOnLight(Triangle tri, ulong* seed){
	if(tri.primType == 1){
		float a = randomUFloat(seed);
		float b = randomUFloat(seed);
		float sqrt_a = sqrt(a);
		// Taken from Shape Distributions, Osada et al.
		float3 p = (1.f - sqrt_a) * tri.vertex0.xyz + sqrt_a * (1.f - b) * tri.vertex1.xyz + sqrt_a * b * tri.vertex2.xyz;
		return p;
	}
	printf("Warning: Unidentified Primitive Type in func: RandomPointOnLight()\n");
	return (float3)(0, 0, 0);
}

#endif