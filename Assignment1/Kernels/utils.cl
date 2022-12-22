#define SCRWIDTH	1280
#define SCRHEIGHT	720



typedef struct __attribute__((aligned(64))) 
{
	// Access all float4 as float3 by using var_name.xyz
	float4 O;	// Ray origin
	float4 D;	// Ray direction (normalized?)
	float4 rD;	// Reciprocal ray direction
	float t;	// Time/distance of hit
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
	uint objId;
} Triangle;

typedef struct __attribute__((aligned(64)))
{
	float3 aabbMin, aabbMax;			// boundary
	uint leftFirst, triCount;			// count and start
} BVHNode;

