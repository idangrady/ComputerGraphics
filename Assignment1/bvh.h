//#pragma once
//
//
//#include<precomp.h>
//
//const int BINS = 10;
//
//static inline BVHNode* pool[20];
//static inline primitives* arrPrimitive[20];
//static inline int triIdx[20];
//
//
//inline float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
//{
//
//	float tx1 = (bmin.x - ray.O.x) / ray.D.x, tx2 = (bmax.x - ray.O.x) / ray.D.x;
//	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
//	float ty1 = (bmin.y - ray.O.y) / ray.D.y, ty2 = (bmax.y - ray.O.y) / ray.D.y;
//	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
//	float tz1 = (bmin.z - ray.O.z) / ray.D.z, tz2 = (bmax.z - ray.O.z) / ray.D.z;
//	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
//	if (tmax >= tmin && tmin < ray.I.t && tmax > 0) return tmin; else return 1e30f;
//}
//
//struct BVHNode
//	// total of 32 bytes. two nodes fit in cashe
//{
//	aabb AABB;
//	uint leftFirst;												// if count >0 its Leafm else-> leftFirst is the left node, right is leftFirst++
//	uint triCount;												// if count >0 the primitives are from idx leftFirst to leftFirst+triCount
//	bool isLeaf() { return triCount > 0; }
//
//	void Traverse(Ray& r)
//	{
//		//if (IntersectAABB(r, AABB.bmin3, AABB.bmax3) == 1e30f)return;
//		if (isLeaf) intersectPrimitives(r);
//		else {
//			pool[leftFirst]->Traverse(r);							// traverse left node
//			pool[leftFirst + 1]->Traverse(r);						// traverse right node
//		}
//	}
//	float3 intersectPrimitives(Ray& r) {
//
//		for (int i = triCount; i >= 0; i--)
//		{
//			arrPrimitive[leftFirst + i]->Intersect(r);
//		}
//	}
//};
//
//struct Bin { aabb bounds; int triCount = 0; };
//
//
//class BVH 
//{
//public:
//
//	
//	static inline primitives* arrPrimitive[20];
//
//	static uint N;
//	BVHNode* bvhNode[20];
//	uint rootNodeIdx = 0, nodesUsed = 1;
//
//
//
//	BVH() = default;
//	~BVH() {};																		// later
//
//	void RefitBVH()
//	{
//		for (int i = nodesUsed - 1; i >= 0; i--) if (i != 1)
//		{
//			BVHNode& node = *bvhNode[i];
//			if (node.isLeaf())
//			{
//				// leaf node: adjust bounds to contained triangles
//				UpdateNodeBounds(i);
//				continue;
//			}
//			// interior node: adjust bounds to child node bounds
//			BVHNode& leftChild = *bvhNode[node.leftFirst];
//			BVHNode& rightChild = *bvhNode[node.leftFirst + 1];
//			node.AABB.bmin3 = fminf(leftChild.AABB.bmin3, rightChild.AABB.bmin3);
//			node.AABB.bmax3 = fmaxf(leftChild.AABB.bmax3, rightChild.AABB.bmax3);
//		}
//	}
//
//	void BuildBVH()
//	{
//		//for (int i = 0; i < N; i++) tri[i].centroid =								Perhas we could 
//		//	(tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.3333f;
//		// assign all triangles to root node
//		BVHNode& root = *bvhNode[rootNodeIdx];
//		root.leftFirst = 0;
//		root.triCount = 0;
//		// subdivide recursively
//		Subdivide(rootNodeIdx);
//	}
//
//	void Subdivide(int nodeIdx);
//	void UpdateNodeBounds(uint nodeIdx)
//	{
//		BVHNode& node = *bvhNode[nodeIdx];
//		node.AABB.bmin3 = float3(1e30f);
//		node.AABB.bmax3 = float3(-1e30f);
//		for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
//		{
//			primitives& prim = *arrPrimitive[first + i];
//			node.AABB.bmin3 = fminf(node.AABB.bmin3, prim.trig->vertex0);				// This needs to change to be adapted to more peimitives. maybe center and vertical and horizontal
//			node.AABB.bmin3 = fminf(node.AABB.bmin3, prim.trig->vertex1);				// This needs to change to be adapted to more peimitives. maybe center and vertical and horizontal
//			node.AABB.bmin3 = fminf(node.AABB.bmin3, prim.trig->vertex2);				// This needs to change to be adapted to more peimitives. maybe center and vertical and horizontal
//			node.AABB.bmax3 = fmaxf(node.AABB.bmax3, prim.trig->vertex0);				// This needs to change to be adapted to more peimitives. maybe center and vertical and horizontal
//			node.AABB.bmax3 = fmaxf(node.AABB.bmax3, prim.trig->vertex1);				// This needs to change to be adapted to more peimitives. maybe center and vertical and horizontal
//			node.AABB.bmax3 = fmaxf(node.AABB.bmax3, prim.trig->vertex2);				// This needs to change to be adapted to more peimitives. maybe center and vertical and horizontal
//		}
//	}	
//	
//	float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos)
//	{
//		float bestCost = 1e30f;
//		for (int a = 0; a < 3; a++)
//		{
//			float boundsMin = 1e30f, boundsMax = -1e30f;
//			for (int i = 0; i < node.triCount; i++)
//			{
//				primitives& triangle = *arrPrimitive[triIdx[node.leftFirst + i]];
//				boundsMin = min(boundsMin, triangle.centroid[a]);
//				boundsMax = max(boundsMax, triangle.centroid[a]);
//			}
//			if (boundsMin == boundsMax) continue;
//			// populate the bins
//			Bin bin[BINS];
//			float scale = BINS / (boundsMax - boundsMin);
//			for (uint i = 0; i < node.triCount; i++)
//			{
//				primitives& primitive = *arrPrimitive[triIdx[node.leftFirst + i]];
//				int binIdx = min(BINS - 1,
//					(int)((primitive.centroid[a] - boundsMin) * scale));
//				bin[binIdx].triCount++;
//				bin[binIdx].bounds.Grow(primitive.trig->vertex0);				// I need to find a better way to do it for all primitives
//				bin[binIdx].bounds.Grow(primitive.trig->vertex1);				// I need to find a better way to do it for all primitives
//				bin[binIdx].bounds.Grow(primitive.trig->vertex2);				// I need to find a better way to do it for all primitives
//			}
//			// gather data for the 7 planes between the 8 bins
//			float leftArea[BINS - 1], rightArea[BINS - 1];
//			int leftCount[BINS - 1], rightCount[BINS - 1];
//			aabb leftBox, rightBox;
//			int leftSum = 0, rightSum = 0;
//			for (int i = 0; i < BINS - 1; i++)
//			{
//				leftSum += bin[i].triCount;
//				leftCount[i] = leftSum;
//				leftBox.Grow(bin[i].bounds);
//				leftArea[i] = leftBox.Area();
//				rightSum += bin[BINS - 1 - i].triCount;
//				rightCount[BINS - 2 - i] = rightSum;
//				rightBox.Grow(bin[BINS - 1 - i].bounds);
//				rightArea[BINS - 2 - i] = rightBox.Area();
//			}
//			// calculate SAH cost for the 7 planes
//			scale = (boundsMax - boundsMin) / BINS;
//			for (int i = 0; i < BINS - 1; i++)
//			{
//				float planeCost =
//					leftCount[i] * leftArea[i] + rightCount[i] * rightArea[i];
//				if (planeCost < bestCost)
//					axis = a, splitPos = boundsMin + scale * (i + 1),
//					bestCost = planeCost;
//			}
//		}
//		return bestCost;
//	}
//
//}
//
//;
//
