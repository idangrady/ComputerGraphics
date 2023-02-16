#include "Kernels/utils.cl"

// __kernel void extend(__global Ray* rays, __constant Triangle* triangles, int triangleCount) {
// 	int threadIdx = get_global_id(0);
//     Ray* ray = &rays[threadIdx];
//     for(int i = 0; i < triangleCount; i++){
//         IntersectTri(ray, &triangles[i], i);
//     }
// }

__kernel void extend(__global Ray* rays, __global Triangle* triangles, __global BVHNode* bvhNodes, __global uint* indicesIntoTriangles){
//__kernel void extend(__global Ray* rays, __global Triangle* triangles, __global BVHNode* bvhNodes, __global uint* indicesIntoTriangles, __global uint* crossedBuffer, __global uint* intersectedTriBuffer){
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
                IntersectTri(ray, &triangles[indicesIntoTriangles[node->leftFirst + i]], indicesIntoTriangles[node->leftFirst + i]);
                //int does = DoesIntersectTri(ray, &triangles[indicesIntoTriangles[node->leftFirst + i]]);
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
}