#include "Kernels/utils.cl"

__kernel void shade(__global Ray *rays, __constant Triangle *triangles,
                    __constant Material *materials,
                    __global float4 *accumulator,
                    volatile __global int *counter, __global Ray *newRays,
                    __constant uint *seedDepth) {
  int threadIdx = get_global_id(0);
  //   if(seedDepth[1] < 1 ) accumulator[rays[threadIdx].pixel] *= (float4)(0.8,
  //   0.8, 0.8, 1); else accumulator[rays[threadIdx].pixel] *= (float4)(1, 0.8,
  //   0, 1); if(threadIdx < 20000 && seedDepth[1] < 10){
  //     Ray newray = rays[threadIdx];
  //     uint old = atomic_add(counter, 1);
  //     newRays[old] = newray;
  //   }
  //   return;
  if (seedDepth[1] > 20) {
    // Max depth
    accumulator[rays[threadIdx].pixel] = (float4)(0, 0, 0, 0);
    //printf("Max depth.\n");
    return;
  }
  uint I = rays[threadIdx].primIdx;
  uint unique_seed = seedDepth[0] + threadIdx;
  //printf("Unique seed: %u\n", unique_seed);
  if (I > 0) {
    uint mI = GetTriangleIndex(I);
    Material m = materials[mI];
    if (m.isEmissive) {
        //printf("Hit light.\n");
      accumulator[rays[threadIdx].pixel] *=
          (float4)(m.albedoSpecularity.xyz, 1);
      return;
    } else {
      Triangle tri = triangles[mI];
      float3 I_loc =
          rays[threadIdx].O.xyz + rays[threadIdx].rD_t.w * rays[threadIdx].D.xyz;
      float3 N = (float3)(tri.N0, tri.N1, tri.N2);
      float3 BRDF = m.albedoSpecularity.xyz * PI;
      //printf("Hit triangle. Getting diffuse reflection.\n");
      float3 random_dir = GetDiffuseReflection(N, &unique_seed);
      Ray newray = {
          (float4)(I_loc + 0.0002f * random_dir, 1),
          (float4)(random_dir, 1),
          (float4)(1.0f / random_dir, 1e34f),
          rays[threadIdx].pixel,
          0,
          (float2)(0, 0)};
      newray.D = (float4)(random_dir, 1);
      newray.pixel = rays[threadIdx].pixel;
      uint old = atomic_add(counter, 1);
      newRays[old] = newray;
      float3 val = dot(N, random_dir) * BRDF * INVPI;
      accumulator[rays[threadIdx].pixel] *= (float4)(val, 1);
      return;
    }
  } else {
    //printf("Hit Nothing.\n");
    accumulator[rays[threadIdx].pixel] = (float4)(0, 0, 0, 0);
    return;
  }
}