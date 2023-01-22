#include "Kernels/utils.cl"
#pragma OPENCL EXTENSION cl_intel_printf : enable



__kernel void constructBVH(__global float3* vertices, __global uint3* indices, __global BVH_GPU* bvh)
    {
        uint i = get_global_id(0);
        float3 v1 = vertices[indices[i].x];
        float3 v2 = vertices[indices[i].y];
        float3 v3 = vertices[indices[i].z];

        // Compute the bounding box of the triangle
        float3 min = v1;
        float3 max = v1;
        min.x = fmin(min.x, v2.x);
        min.y = fmin(min.y, v2.y);
        min.z = fmin(min.z, v2.z);
        max.x = fmax(max.x, v2.x);
        max.y = fmax(max.y, v2.y);
        max.z = fmax(max.z, v2.z);
        min.x = fmin(min.x, v3.x);
        min.y = fmin(min.y, v3.y);
        min.z = fmin(min.z, v3.z);
        max.x = fmax(max.x, v3.x);
        max.y = fmax(max.y, v3.y);
        max.z = fmax(max.z, v3.z);


        // Store the bounding box in the BVH node
        bvh[i].aabbMin = (float4)(min.x, min.y, min.z, 0);
        bvh[i].aabbMax = (float4)(max.x, max.y, max.z, 0);
        bvh[i].leftFirst = i;
        bvh[i].triCount = 1;

    }
