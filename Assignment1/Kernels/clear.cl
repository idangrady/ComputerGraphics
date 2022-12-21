__kernel void clear(__global float4* accumulator){
    int threadIdx = get_global_id(0);
    accumulator[threadIdx] = (float4)(1, 1, 1, 1);
}