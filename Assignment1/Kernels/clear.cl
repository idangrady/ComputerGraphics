__kernel void clear(__global float4* intermediate){
    int threadIdx = get_global_id(0);
    intermediate[threadIdx] = (float4)(1, 1, 1, 1);
}

__kernel void resetAccumulator(__global float4* accumulator){
    int threadIdx = get_global_id(0);
    accumulator[threadIdx] = (float4)(0, 0, 0, 0);
}