#include "Kernels/utils.cl"

__kernel void renderToScreen(__global float4* accumulator, write_only image2d_t target){
    int threadIdx = get_global_id(0);
    int x = threadIdx % SCRWIDTH;
	int y = threadIdx / SCRWIDTH;
    write_imagef( target, (int2)(x, y), accumulator[threadIdx]);
}