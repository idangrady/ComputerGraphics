#include "Kernels/utils.cl"

__kernel void connect(__global Ray* shadowRays, __global Triangle* triangles, __global TriEx *triExes,
                      __constant Material *materials, __global BVHNode* bvhNodes, __global uint* indicesIntoTriangles,
                      __global float4* intermediate, __global float4* accumulator, __global float* oldPDFs){
//__kernel void extend(__global Ray* rays, __global Triangle* triangles, __global BVHNode* bvhNodes, __global uint* indicesIntoTriangles, __global uint* crossedBuffer, __global uint* intersectedTriBuffer){
    int threadIdx = get_global_id(0);
    Ray* ray = &shadowRays[threadIdx];
    float tLight = ray->rD_t.w;
    BVHNode* node = &bvhNodes[0], *stack[32];
    Triangle* lightPrim = &triangles[GetTriangleIndex(ray->primIdx)];
    TriEx* lightEx = &triExes[GetTriangleIndex(ray->primIdx)];
    if(!(ray->uv.x > 0 && dot(lightEx->N.xyz, -ray->D.xyz))) return;
    __constant Material* lightMat = &materials[lightEx->matId];
    uint stackPtr = 0;
    while(1)
    {
        if(node->triCount > 0) //leaf
        {
            for(uint i = 0; i < node->triCount; i++)
            {
                //IntersectTri(ray, &triangles[indicesIntoTriangles[node->leftFirst + i]], indicesIntoTriangles[node->leftFirst + i]);
                if(ShadowRayIntersect(ray, &triangles[indicesIntoTriangles[node->leftFirst + i]])){
                    return;
                };
                //intersectedTriBuffer[threadIdx] += does;
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
           // crossedBuffer[threadIdx] += 1;
            node = child1;
            if(dist2 != 1e30f){
                //crossedBuffer[threadIdx] += 1;
                stack[stackPtr++] = child2;
            } 
        }
    }
    float summedPDF = oldPDFs[ray->pixel];
    float solidAngle = (dot(lightEx->N.xyz, -ray->D.xyz) * lightEx->A) / (ray->rD_t.w * ray->rD_t.w);
    float lightPDF = 1.f / solidAngle;
    summedPDF += lightPDF;
    // We hide the N dot L from the intersected triangle in the uv coordinates of the ray
    accumulator[ray->pixel] += intermediate[ray->pixel] * (float4)((ray->uv.x * (1.f / summedPDF) * lightMat->albedoSpecularity.xyz), 1);
}