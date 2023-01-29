#include "Kernels/utils.cl"

void IntersectTri(Ray* ray, __global Triangle* tri, int id) 
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
        ray->primIdx = makeId(1, 0, id);
    }
}

int DoesIntersectTri(Ray* ray, __global Triangle* tri){
    const float3 edge1 = tri->vertex1.xyz - tri->vertex0.xyz;
    const float3 edge2 = tri->vertex2.xyz - tri->vertex0.xyz;
    const float3 h = cross(ray->D.xyz, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return 0; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O.xyz - tri->vertex0.xyz;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return 0;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D.xyz, q);
    if (v < 0 || u + v > 1) return 0;
    const float t = f * dot(edge2, q);
    if(t > 0.0001f) return 1;
    else return 0;
}

// __kernel void extend(__global Ray* rays, __constant Triangle* triangles, int triangleCount) {
// 	int threadIdx = get_global_id(0);
//     Ray* ray = &rays[threadIdx];
//     for(int i = 0; i < triangleCount; i++){
//         IntersectTri(ray, &triangles[i], i);
//     }
// }

float IntersectAABB( Ray* ray, BVHNode* node )
{
    float tx1 = (node->minx - ray->O.x) * ray->rD_t.x, tx2 = (node->maxx - ray->O.x) * ray->rD_t.x;
    float tmin = min( tx1, tx2 ), tmax = max( tx1, tx2 );
    float ty1 = (node->miny - ray->O.y) * ray->rD_t.y, ty2 = (node->maxy - ray->O.y) * ray->rD_t.y;
    tmin = max( tmin, min( ty1, ty2 ) ), tmax = min( tmax, max( ty1, ty2 ) );
    float tz1 = (node->minz - ray->O.z) * ray->rD_t.z, tz2 = (node->maxz - ray->O.z) * ray->rD_t.z;
    tmin = max( tmin, min( tz1, tz2 ) ), tmax = min( tmax, max( tz1, tz2 ) );
    if (tmax >= tmin && tmin < ray->rD_t.w && tmax > 0) return tmin; else return 1e30f;
}


//__kernel void extend(__global Ray* rays, __global Triangle* triangles, __global BVHNode* bvhNodes, __global uint* indicesIntoTriangles){
__kernel void extend(__global Ray* rays, __global Triangle* triangles, __global BVHNode* bvhNodes, __global uint* indicesIntoTriangles, __global uint* crossedBuffer, __global uint* intersectedTriBuffer){
    int threadIdx = get_global_id(0);
    Ray* ray = &rays[threadIdx];
    BVHNode* node = &bvhNodes[0], *stack[32];
    uint stackPtr = 0;
    while(1)
    {
        if(node->triCount > 0) //leaf
        {
            for(uint i = 0; i < node->triCount; i++)
            {
                //IntersectTri(ray, &triangles[indicesIntoTriangles[node->leftFirst + i]], indicesIntoTriangles[node->leftFirst + i]);
                int does = DoesIntersectTri(ray, &triangles[indicesIntoTriangles[node->leftFirst + i]]);
                intersectedTriBuffer[threadIdx] += does;
            }
            if(stackPtr == 0) break; else node = stack[--stackPtr];
            continue;
        }
        BVHNode* child1 = &bvhNodes[node->leftFirst];
        BVHNode* child2 = &bvhNodes[node->leftFirst + 1];
        float dist1 = IntersectAABB(ray, child1);
        float dist2 = IntersectAABB(ray, child2);
        if(dist1 > dist2)
        {
            float d = dist1; dist1 = dist2; dist2 = d;
            BVHNode* c = child1; child1 = child2; child2 = c;
        }
        if(dist1 == 1e30f)
        {
            if(stackPtr == 0) break; else node = stack[--stackPtr];   
        }
        else
        {
            crossedBuffer[threadIdx] += 1;
            node = child1;
            if(dist2 != 1e30f){
                crossedBuffer[threadIdx] += 1;
                stack[stackPtr++] = child2;
            } 
        }
    }
}