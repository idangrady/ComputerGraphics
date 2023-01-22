#include "Kernels/utils.cl"

typedef struct {
    float4 min;
    float4 max;
} BoundingBox;


inline void UpdateNodeBounds(int nodeIdx,  uint* arrPrimitiveIdx,  BVH_GPU* bvh,  Triangle* arrPrimitive) {
    BVH_GPU node = bvh[nodeIdx];
    node.aabbMin = (float4)(1e30f);
    node.aabbMax = (float4)(-1e30f);
    for (uint first = node.leftFirst, i = 0; i < node.triCount; i++) {
        uint curr_val = arrPrimitiveIdx[first + i];
        Triangle* ssass = arrPrimitive[curr_val];
        pair<float4, float4> aabb_minMax = createAABB(ssass);
        for(int j = 0; j < 4; j++) {
            node.aabbMin.s[j] = fmin(node.aabbMin.s[j], aabb_minMax.first.s[j]);
            node.aabbMax.s[j] = fmax(node.aabbMax.s[j], aabb_minMax.second.s[j]);
        }
    }
    bvh[nodeIdx] = node;
}

__kernel void BuildBVH(
    __global BVH_GPU* bvh,
    __global uint* arrPrimitiveIdx,
    __global Triangle* arrPrimitive,
    int count
    )
{
    int rootNodeIdx = 0;
    int nodesUsed = 1;
    // assign all triangles to root node
    BVH_GPU root = bvh[rootNodeIdx];
    root.leftFirst = 0;
    root.triCount = count; // amount of primitives

    UpdateNodeBounds(rootNodeIdx);
    // create a stack to store the indices of the nodes
    __global int* stack = arrPrimitiveIdx;
    int top = -1;
    // push the root node index to the stack
    stack[++top] = rootNodeIdx;
    // while the stack is not empty
    while (top >= 0) {
        
        // pop the top node index from the stack
        int nodeIdx = stack[top--];
        BVH_GPU node = bvh[nodeIdx];
        // if the node contains only one triangle
        if (node.triCount <= 1) {
            continue;
        }
        // determine split axis and position
        float4 extent = node.aabbMax - node.aabbMin;
        int axis = 0;
        if (extent.y > extent.x) axis = 1;
        if (extent.z > extent.w) axis = 2;
        float splitPos = node.aabbMin.w + extent.w * 0.5f;
        // in-place partition
        int i = node.leftFirst;
        int j = i + node.triCount - 1;

          //      while (i <= j)
        //{
        float4 center = (float4)(arrPrimitive[arrPrimitiveIdx[i].vertex0.w,arrPrimitive[arrPrimitiveIdx[i].vertex1.w, arrPrimitive[arrPrimitiveIdx[i].vertex2.w,1)
            //if (center[axis] < splitPos)
             //   i++;
           // else
          //      swap(arrPrimitiveIdx[i], arrPrimitiveIdx[j--]);
        //}
        // create left and right child nodes
        int leftChildIdx = nodesUsed++;
        int rightChildIdx = nodesUsed++;
        bvh[leftChildIdx].leftFirst = node.leftFirst;
        bvh[leftChildIdx].triCount = i - node.leftFirst;
        bvh[rightChildIdx].leftFirst = i;
        bvh[rightChildIdx].triCount = node.triCount - (i - node.leftFirst);
        UpdateNodeBounds(leftChildIdx);
        UpdateNodeBounds(rightChildIdx);

        // push the left and right child nodes to the stack
        stack[++top] = leftChildIdx;
        stack[++top] = rightChildIdx;
    }
 }