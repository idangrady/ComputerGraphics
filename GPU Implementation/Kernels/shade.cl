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

float3 GetTexel(uint* texture, TextureData textureData, float2 uv_I, TriEx tri){
  int width = textureData.width;
  int height = textureData.height;
  float2 uv = uv_I.x * tri.uv1 + uv_I.y * tri.uv2 + (1.0f - (uv_I.x + uv_I.y)) * tri.uv0;
  int iu = (int)(uv.x * width) % width;
  int iv = (int)(uv.y * height) % height;
  uint texel = texture[iu + iv * width];
  uchar x = (texel >> 16);
  uchar y = (texel >> 8);
  uchar z = texel;
  return (float3)(x / 255.f, y / 255.f, z / 255.f);
}

__kernel void shade(__global Ray *rays, __global Triangle* tris,
                    __global TriEx *triExes, __constant Material *materials,
                    __global float4 *intermediate, __global float* oldPDFs,
                    volatile __global int *counter, __global Ray *newRays, __global Ray *shadowRays,
                    __global ulong *seeds, __global uint *depth,
                    __constant float4* skybox, __private int skyBoxWidth, __private int skyBoxHeight,
                    //__global uint* textures, __constant TextureData* textureData, __constant int* textureIndices, __global uint* crossedBuffer,  __global uint* intersectedTriBuffer) {
                    __global uint* textures, __constant TextureData* textureData, __constant int* textureIndices,
                    __constant uint* lightIndices, uint lightCount) {
  int threadIdx = get_global_id(0);
  Ray* ray = &rays[threadIdx];
  //intermediate[ray->pixel] = (float4)(2* intersectedTriBuffer[threadIdx], 2 * crossedBuffer[threadIdx] / 255.f, 0, 0);
  //return;
  if (depth[0] > 10) {
    // Max depth
    intermediate[ray->pixel] = (float4)(0, 0, 0, 0);
    return;
  }
  float oldPDF = oldPDFs[ray->pixel];
  float INV_OldPDF = 1.0f / oldPDFs[ray->pixel];
  ulong unique_seed = seeds[threadIdx];
  // printf("Unique seed: %u\n", unique_seed);
  uint I = ray->primIdx;
  if (I == 0) {
    // printf("Hit Nothing.\n");
    //intermediate[ray->pixel] = (float4)(0, 0, 0, 0);
    // Postponed PDF division done here
    intermediate[ray->pixel] *= INV_OldPDF * getSkyBox(ray->D.xyz, skyBoxWidth, skyBoxHeight, skybox);
    return;
  }
  uint mI = GetTriangleIndex(I);
  //printf("%i\r\n", mI);
  TriEx tri = triExes[mI];
  float3 I_loc =
      ray->O.xyz + ray->rD_t.w * ray->D.xyz;
  Material m = materials[tri.matId];
  float d = 1.0f - m.albedoSpecularity.w;
  float3 N = tri.N.xyz;
  //if(mI == 0) printf("Left wall normal: {%f, %f, %f}\n", tri.N.x, tri.N.y, tri.N.z);
  float r = randomFloat(&unique_seed);
  bool hit_back = false;
  if (dot(N, ray->D.xyz) > 0) {
    hit_back = true;
    N = -N;
  }
  if (m.isEmissive) {
    if(hit_back) intermediate[ray->pixel] = (float4)(0, 0, 0, 0);
    // printf("Hit light.\n");
    else{
      float solidAngle = (dot(N, -ray->D.xyz) * tri.A) / (ray->rD_t.w * ray->rD_t.w);
      float lightPDF = 1.f / solidAngle;
      float INV_combinedPDF = 1.f / (oldPDF + lightPDF);
      //if(oldPDF == 1.0f) printf("lightPDF:\t%f\r\n", lightPDF);
      // Postponed PDF division done here
      intermediate[ray->pixel] *= INV_combinedPDF * (float4)(m.albedoSpecularity.xyz, 1);
    } 
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
    float traveled = ray->rD_t.w;
    float cos_theta_i = dot(N, -ray->D.xyz);
    float k = GetSnell(refr, cos_theta_i);
    float R = 1.0f;
    float T = 0.f;
    if (k > 0.00001f) {
      R = FresnelReflection(cos_theta_i, n1, n2, refr, N,
                            ray->D.xyz);
      if (r > R) { // Randomly refracted
        float3 t_dir = (refr * ray->D.xyz) +
                       (N * (refr * cos_theta_i - sqrt(k)));
        t_dir = normalize(t_dir);
        Ray newray = {(float4)(I_loc + (0.0002f * t_dir), 1),
                      (float4)(t_dir, 1),
                      (float4)(1.0f / t_dir, 1e34f),
                      ray->pixel,
                      0,
                      (float2)(0, 0)};
        uint old = atomic_add(counter, 1);
        //if(old > skyBoxWIDTH * SCRHEIGHT) printf("Error: Huge counter!: %i\n", old);
        newRays[old] = newray;
      } else { // Randomly reflected
        Ray newray = reflectRay(*ray, I_loc, N);
        int old = atomic_add(counter, 1);
        newRays[old] = newray;
      }
    } else { // just reflect
      Ray newray = reflectRay(*ray, I_loc, N);
      int old = atomic_add(counter, 1);
      //if(old > skyBoxWIDTH * SCRHEIGHT) printf("Error: Huge counter!: %i\n", old);
      newRays[old] = newray;
    }
    if (hit_back) {
      intermediate[ray->pixel].x *= exp(-m.absorption.x * traveled);
      intermediate[ray->pixel].y *= exp(-m.absorption.y * traveled);
      intermediate[ray->pixel].z *= exp(-m.absorption.z * traveled);
    }
  } else {
    if (r > d) {
      Ray newray = reflectRay(rays[threadIdx], I_loc, N);
      uint old = atomic_add(counter, 1);
      //if(old > skyBoxWIDTH * SCRHEIGHT) printf("Error: Huge counter!: %i\n", old);
      newRays[old] = newray;
    } else {
      int old_shadow = atomic_add(counter + 1, 1);
      uint randomLightID = lightIndices[RandomLight(lightCount, &unique_seed)];
      Triangle randomLight = tris[randomLightID];
      float3 random_point = RandomPointOnLight(randomLight, &unique_seed);
      //printf("Random Lightpoint:\t%2.2v3hlf\n", random_point);
      float3 l_dir = normalize(random_point - I_loc);
      Ray shadowray = {
        (float4)(I_loc + 0.0002f * l_dir, 1),
        (float4)(l_dir, 1),
        (float4)(1.0f / l_dir, length(random_point - I_loc + 0.0002f * l_dir)),
        ray->pixel,
        makeId(1, 0, randomLightID),
        (float2)(dot(N, l_dir), 0)
      };
      shadowRays[old_shadow] = shadowray;
      // We can finally do the postponed PDF division from previous bounce
      intermediate[ray->pixel] *= INV_OldPDF;
      float3 BRDF = INVPI;
      if(tri.textureID < 0) BRDF *= m.albedoSpecularity.xyz;
      else {
        BRDF *= GetTexel(&(textures[textureIndices[tri.textureID]]), textureData[tri.textureID], ray->uv, tri);
      }
      // printf("Hit triangle. Getting diffuse reflection.\n");
      float3 random_dir = GetDiffuseReflection(N, &unique_seed);
      Ray newray = {(float4)(I_loc + 0.0002f * random_dir, 1),
                    (float4)(random_dir, 1),
                    (float4)(1.0f / random_dir, 1e34f),
                    ray->pixel,
                    0,
                    (float2)(0, 0)};
      int old = atomic_add(counter, 1);
      //if(old > skyBoxWIDTH * SCRHEIGHT) printf("Error: Huge counter!: %i\n", old);
      newRays[old] = newray;
      // Insert the hemiPDF in the old PDF array for postponing
      oldPDFs[ray->pixel] = 1.f / (PI * 2.f);
      float3 val = dot(N, random_dir) * BRDF;// * INVPDF; PDF postponed!
      intermediate[ray->pixel] *= (float4)(val, 1);
    }
  }
  seeds[threadIdx] = unique_seed;
}