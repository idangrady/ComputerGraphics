#include "template/common.h"
#include "Kernels/utils.cl"

inline uint makeId(uint type, uint id, uint tri){
    return (type << 30) + (id << 20) + (tri);
}

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
    if (t > 0.0001f && t < ray->t){
        ray->t = t;
        ray->uv = (float2)(u, v);
        ray->primIdx = makeId(1, 0, tri->id);
    }
}

__kernel void extend(__global Ray* rays, __constant Triangle* triangles, int triangleCount) {
	int threadIdx = get_global_id(0);
    Ray* ray = &rays[threadIdx];
    for(int i = 0; i < 2; i++){
        IntersectTri(ray, &triangles[i]);
    }
}
