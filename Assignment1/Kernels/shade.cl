#include "template/common.h"
#include "Kernels/utils.cl"

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

__kernel void shade(__global Ray* rays, __constant Triangle* triangles, __constant float4* triangleColors, __global float4* accumulator){
	int threadIdx = get_global_id(0);
    uint I = rays[threadIdx].primIdx;
    if(GetObjectType(I) > 0)
    {
        uint mI = GetTriangleIndex(I);
        accumulator[threadIdx] = triangleColors[mI];
    }
    else{
        accumulator[threadIdx] = (float4)(0, 0, 0, 0);
    }
}