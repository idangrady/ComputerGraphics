#include "Kernels/utils.cl"

float FresnelReflection(float cos_theta_i, float n1, float n2, float refr,
                        float3 N, float3 D) {
  float sin_theta_i = length(cross(N, -D));
  float refr_sin = (refr * sin_theta_i);
  float cos_theta_t = sqrt(1.0f - (refr_sin * refr_sin));
  float first_term = ((n1 * cos_theta_i) - (n2 * cos_theta_t)) /
                     ((n1 * cos_theta_i) + (n2 * cos_theta_t));
  float second_term = ((n1 * cos_theta_t) - (n2 * cos_theta_i)) /
                      ((n1 * cos_theta_t) + (n2 * cos_theta_i));
  return 0.5f * ((first_term * first_term) + (second_term * second_term));
}

float GetSnell(float refr, float cos_theta_i) {
  float k = 1.0f - (refr * refr) * (1.0f - (cos_theta_i * cos_theta_i));
  return k;
}

__kernel void shade(__global Ray *rays, __constant Triangle *triangles,
                    __constant TriEx *triExes, __constant Material *materials,
                    __global float4 *intermediate,
                    volatile __global int *counter, __global Ray *newRays,
                    __global ulong *seeds, __constant uint *depth,
                    __constant float4* skybox, int width, int height) {
  int threadIdx = get_global_id(0);
  if (depth[0] > 20) {
    // Max depth
    intermediate[rays[threadIdx].pixel] = (float4)(0, 0, 0, 0);
    return;
  }
  ulong unique_seed = seeds[threadIdx];
  // printf("Unique seed: %u\n", unique_seed);
  uint I = rays[threadIdx].primIdx;
  if (I == 0) {
    // printf("Hit Nothing.\n");
    //intermediate[rays[threadIdx].pixel] = (float4)(0, 0, 0, 0);
    intermediate[rays[threadIdx].pixel] *= getSkyBox(rays[threadIdx].D.xyz, width, height, skybox);
    return;
  }
  uint mI = GetTriangleIndex(I);
  Triangle tri = triangles[mI];
  float3 I_loc =
      rays[threadIdx].O.xyz + rays[threadIdx].rD_t.w * rays[threadIdx].D.xyz;
  Material m = materials[triExes[mI].matId];
  float d = 1.0f - m.albedoSpecularity.w;
  float3 N = (float3)(tri.N0, tri.N1, tri.N2);
  float r = randomFloat(&unique_seed);
  bool hit_back = false;
  if (dot(N, rays[threadIdx].D.xyz) > 0) {
    hit_back = true;
    N = -N;
  }
  if (m.isEmissive) {
    // printf("Hit light.\n");
    intermediate[rays[threadIdx].pixel] *= (float4)(m.albedoSpecularity.xyz, 1);
  } else if (m.medium > 0) {
    float n1, n2;
    if (hit_back) {
      n1 = 1.52f;
      n2 = 1.0f;
    } else {
      n1 = 1.0f;
      n2 = 1.52f;
    }
    float refr = n1 / n2;
    float traveled = rays[threadIdx].rD_t.w;
    float cos_theta_i = dot(N, -rays[threadIdx].D.xyz);
    float k = GetSnell(refr, cos_theta_i);
    float R = 1.0f;
    float T = 0.f;
    if (k > 0.00001f) {
      R = FresnelReflection(cos_theta_i, n1, n2, refr, N,
                            rays[threadIdx].D.xyz);
      if (r > R) { // Randomly refracted
        float3 t_dir = (refr * rays[threadIdx].D.xyz) +
                       (N * (refr * cos_theta_i - sqrt(k)));
        t_dir = normalize(t_dir);
        Ray newray = {(float4)(I_loc + (0.0002f * t_dir), 1),
                      (float4)(t_dir, 1),
                      (float4)(1.0f / t_dir, 1e34f),
                      rays[threadIdx].pixel,
                      0,
                      (float2)(0, 0)};
        uint old = atomic_add(counter, 1);
        if(old > SCRWIDTH * SCRHEIGHT) printf("Huge counter!: %f", old);
        newRays[old] = newray;
      } else { // Randomly reflected
        Ray newray = reflectRay(rays[threadIdx], I_loc, N);
        int old = atomic_add(counter, 1);
        newRays[old] = newray;
      }
    } else { // just reflect
      Ray newray = reflectRay(rays[threadIdx], I_loc, N);
      int old = atomic_add(counter, 1);
      if(old > SCRWIDTH * SCRHEIGHT) printf("Huge counter!: %f", old);
      newRays[old] = newray;
    }
    if (hit_back) {
      intermediate[rays[threadIdx].pixel].x *= exp(-m.absorption.x * traveled);
      intermediate[rays[threadIdx].pixel].y *= exp(-m.absorption.y * traveled);
      intermediate[rays[threadIdx].pixel].z *= exp(-m.absorption.z * traveled);
    }
  } else {
    if (r > d) {
      Ray newray = reflectRay(rays[threadIdx], I_loc, N);
      uint old = atomic_add(counter, 1);
      if(old > SCRWIDTH * SCRHEIGHT) printf("Huge counter!: %f", old);
      newRays[old] = newray;
    } else {
      float3 BRDF = m.albedoSpecularity.xyz * PI;
      // printf("Hit triangle. Getting diffuse reflection.\n");
      float3 random_dir = GetDiffuseReflection(N, &unique_seed);
      Ray newray = {(float4)(I_loc + 0.0002f * random_dir, 1),
                    (float4)(random_dir, 1),
                    (float4)(1.0f / random_dir, 1e34f),
                    rays[threadIdx].pixel,
                    0,
                    (float2)(0, 0)};
      int old = atomic_add(counter, 1);
      if(old > SCRWIDTH * SCRHEIGHT) printf("Huge counter!: %f", old);
      newRays[old] = newray;
      float3 val = dot(N, random_dir) * BRDF * INVPI;
      intermediate[rays[threadIdx].pixel] *= (float4)(val, 1);
    }
  }
  seeds[threadIdx] = unique_seed;
}