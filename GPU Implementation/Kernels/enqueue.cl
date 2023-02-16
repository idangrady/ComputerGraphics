#include "Kernels/utils.cl"

__kernel void copy(__global Ray* rays, __global Ray *newRays)
{
    int threadIdx = get_global_id(0);
    rays[threadIdx] = newRays[threadIdx];
}