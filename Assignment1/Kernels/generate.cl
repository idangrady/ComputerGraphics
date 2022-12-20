#include "template/common.h"
#include "Kernels/utils.cl"

Ray createPrimaryRay(__constant Camera* c, int i) {
	int x = i % SCRWIDTH;
	int y = i / SCRWIDTH;
	// calculate pixel position on virtual screen plane
	float u = (float)x * (1.0f / SCRWIDTH);
	float v = (float)y * (1.0f / SCRHEIGHT);
	float3 P = c->topLeft.xyz + u * (c->topRight.xyz - c->topLeft.xyz) + v * (c->bottomLeft.xyz - c->topLeft.xyz);
	float3 D = normalize(P - c->camPos.xyz);
	Ray ray = { 
		(float4)(c->camPos.xyz, 1), 
		(float4)(D, 1),
		(float4)(1.0f / D, 1),
		1e34f,
		0,
		(float2)(0, 0)
	};
	return ray;
}

//__attribute__ ((reqd_work_group_size(32, 1, 1)))
__kernel void generate(__global Ray* rays, __constant Camera* c) {
	int threadIdx = get_global_id(0);
	//int workerIdx = get_global_id(1);
	Ray ray = createPrimaryRay(c, threadIdx);
	rays[threadIdx] = ray;
}