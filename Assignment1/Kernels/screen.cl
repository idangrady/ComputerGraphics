#include "Kernels/utils.cl"

__kernel void renderToScreen(__global float4* intermediate, __global float4* accumulator, write_only image2d_t target, __constant int* framesSinceMoved){
    int threadIdx = get_global_id(0);
    int x = threadIdx % SCRWIDTH;
	int y = threadIdx / SCRWIDTH;
    accumulator[threadIdx].xyz += intermediate[threadIdx].xyz;
    write_imagef( target, (int2)(x, y), (float4)(accumulator[threadIdx].xyz / framesSinceMoved[0], 1));
}