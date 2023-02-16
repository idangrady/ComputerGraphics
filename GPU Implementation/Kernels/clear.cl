__kernel void clear(__global float4* intermediate, __global float* oldPDFs){
//__kernel void clear(__global float4* intermediate, __global uint* crossedBuffer, __global uint* intersectedTriBuffer){
    int threadIdx = get_global_id(0);
    intermediate[threadIdx] = (float4)(1, 1, 1, 1);
    oldPDFs[threadIdx] = 1.0f;
    //crossedBuffer[threadIdx] = 0;
    //intersectedTriBuffer[threadIdx] = 0;
}

__kernel void resetAccumulator(__global float4* accumulator){
    int threadIdx = get_global_id(0);
    accumulator[threadIdx] = (float4)(0, 0, 0, 0);
}