#include "Kernels/utils.cl"

void IntersectTri(Ray* ray, __constant Triangle* tri) 
{
    const float3 edge1 = tri->vertex1.xyz - tri->vertex0.xyz;
    const float3 edge2 = tri->vertex2.xyz - tri->vertex0.xyz;
    const float3 h = cross(ray->D.xyz, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O.xyz - tri->vertex0.xyz;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D.xyz, q);
    if (v < 0 || u + v > 1) return;
    const float t = f * dot(edge2, q);
    if (t > 0.0001f && t < ray->rD_t.w){
        ray->rD_t.w = t;
        ray->uv = (float2)(u, v);
        ray->primIdx = makeId(1, 0, tri->id);
    }
}

__kernel void extend(__global Ray* rays, __constant Triangle* triangles,__constant int*arrPrimitivesIdx , __constant BVHNode* bvhnodes, int triangleCount) 
 {
	int threadIdx = get_global_id(0);
    Ray* ray = &rays[threadIdx];
    for(int i = 0; i < triangleCount; i++){
        IntersectTri(ray, &triangles[i]);
    }
}
//__kernel void extend(__global Ray* rays, __constant Triangle* triangles, int triangleCount)