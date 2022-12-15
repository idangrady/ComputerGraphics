//#pragma once
//#include "precomp.h"
//
//
//#define N 30
//
//
//static inline BVHNode* pool[N];
//static inline primitives* arrPrimitive[N];
//int arrPrimitiveIdx[N];
//
//
//
//struct BVHNode
//{
//    float3 aabbMin, aabbMax;
//    uint leftFirst, triCount;
//}; 
//
//
//class BVH {
//    BVHNode *bvhNode[N * 2 - 1];
//    uint rootNodeIdx = 0, nodesUsed = 1;
//
//    BVH() = default;
//    ~BVH() {};
//
//    BVH(BVHNode* arr[]) 
//    {
//        for (int i = N; i > 0; i--) 
//        {
//            // assigning the poiters for the 
//            bvhNode[i] = arr[i];
//        }
//    }
//    void BuildBVH()
//    {
//        //for (int i = 0; i < N; i++) arrPrimitive[i]->centroid =
//        //    (arrPrimitive[i][i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.3333f;
//
//        // assign all triangles to root node
//        BVHNode& root = *bvhNode[rootNodeIdx];
//        root.leftFirst  = 0;
//        root.triCount = N;
//        UpdateNodeBounds(rootNodeIdx);
//        // subdivide recursively
//        Subdivide(rootNodeIdx);
//    }
//
//    void UpdateNodeBounds(uint nodeIdx)
//    {
//        BVHNode& node = *bvhNode[nodeIdx];
//        node.aabbMin = float3(1e30f);
//        node.aabbMax = float3(-1e30f);
//        for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
//        {
//            primitives& leafTri = *arrPrimitive[first + i];
//            node.aabbMin = fminf(node.aabbMin, leafTri.trig->vertex0);
//            node.aabbMin = fminf(node.aabbMin, leafTri.trig->vertex1);
//            node.aabbMin = fminf(node.aabbMin, leafTri.trig->vertex2);
//            node.aabbMax = fmaxf(node.aabbMax, leafTri.trig->vertex0);
//            node.aabbMax = fmaxf(node.aabbMax, leafTri.trig->vertex1);
//            node.aabbMax = fmaxf(node.aabbMax, leafTri.trig->vertex2);
//        }
//    }
//
//    void Subdivide(uint nodeIdx)
//    {
//        // terminate recursion
//        BVHNode& node = *bvhNode[nodeIdx];
//        if (node.triCount <= 2) return;
//        // determine split axis and position
//        float3 extent = node.aabbMax - node.aabbMin;
//        int axis = 0;
//        if (extent.y > extent.x) axis = 1;
//        if (extent.z > extent[axis]) axis = 2;
//        float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
//        // in-place partition
//        int i = node.leftFirst;
//        int j = i + node.triCount - 1;
//        while (i <= j)
//        {
//            if (arrPrimitive[arrPrimitiveIdx[i]]->trig->centroid[axis] < splitPos)
//                i++;
//            else
//                swap(arrPrimitiveIdx[i], arrPrimitiveIdx[j--]);
//        }
//        // abort split if one of the sides is empty
//        int leftCount = i - node.leftFirst;
//        if (leftCount == 0 || leftCount == node.triCount) return;
//        // create child nodes
//        int leftChildIdx = nodesUsed++;
//        int rightChildIdx = nodesUsed++;
//        node.leftFirst = leftChildIdx;
//        bvhNode[leftChildIdx]->leftFirst = node.leftFirst;
//        bvhNode[leftChildIdx]->triCount = leftCount;
//        bvhNode[rightChildIdx]->leftFirst = i;
//        bvhNode[rightChildIdx]->triCount = node.triCount - leftCount;
//        node.triCount = 0;
//        UpdateNodeBounds(leftChildIdx);
//        UpdateNodeBounds(rightChildIdx);
//        // recurse
//        Subdivide(leftChildIdx);
//        Subdivide(rightChildIdx);
//    }
//
//};
//
