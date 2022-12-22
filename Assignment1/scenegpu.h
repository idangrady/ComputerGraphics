#pragma 
#include <../lib/stb_image.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <config.h>

namespace Tmpl8 {
	__declspec(align(64)) struct TriGPU {
		cl_float4 vertex0, vertex1, vertex2;
		cl_float N0, N1, N2;
		cl_uint id;
	};


	__declspec(align(64)) struct BVH_GPU
	{
		float4 aabbMin, aabbMax;			// boundary
		uint leftFirst, triCount;			// count and start
	};

	__declspec(align(64)) struct Primitive_GPU {
		float4 vertex0, vertex1, vertex2;
		cl_float N0, N1, N2;
		cl_uint id;
		cl_uint objId;
	};

	float4 getCentroid(Primitive_GPU primitive)
	{
		if (primitive.objId == 1)
		{
			auto x = (primitive.vertex0.x + primitive.vertex1.x + primitive.vertex2.x) / 3;
			auto y = (primitive.vertex0.y + primitive.vertex1.y + primitive.vertex2.y) / 3;
			auto z = (primitive.vertex0.z + primitive.vertex1.z + primitive.vertex2.z) / 3;
			return float4(x, y, z, 1);
		}

	}

	pair<float4, float4>  createAABB(Primitive_GPU* primitive_GPU)
	{
		float4 aabbMin = float4(1e30f);
		float4 aabbMax = float4(-1e30f);
		Primitive_GPU* curr_primitive = primitive_GPU;
		if (primitive_GPU->objId == 0) {
			// Triangle
			aabbMin = fminf(aabbMin, primitive_GPU->vertex0);
			aabbMin = fminf(aabbMin, primitive_GPU->vertex1);
			aabbMin = fminf(aabbMin, primitive_GPU->vertex2);
			aabbMax = fmaxf(aabbMax, primitive_GPU->vertex0);
			aabbMax = fmaxf(aabbMax, primitive_GPU->vertex1);
			aabbMax = fmaxf(aabbMax, primitive_GPU->vertex2);
		}
		else if (primitive_GPU->objId == 1) {
			// sphere
			auto radius = primitive_GPU->vertex1.x;
			aabbMin = getCentroid(*primitive_GPU) - float4(radius, radius, radius, radius);
			aabbMax = getCentroid(*primitive_GPU) + float4(radius, radius, radius, radius);
		}
		else if (primitive_GPU->objId == 2) {
			// Cube
		}
		return std::make_pair(aabbMin, aabbMax);
	}


	class SceneGPU {
	public:


		//TriGPU* tris;
		cl_float4* triColors;
		Primitive_GPU* arrPrimitive;
		cl_uint arrPrimitiveIdx[NumberOfElements];
		BVH_GPU* bvhNode;
		cl_uint rootNodeIdx;
		cl_uint nodesUsed;
		SceneGPU() {
			// init
			bvhNode = new BVH_GPU[2 * NumberOfElements - 1];
			arrPrimitive = new Primitive_GPU[NumberOfElements];
			triColors = new cl_float4[NumberOfElements];
			rootNodeIdx = 0; nodesUsed = 0;

			// First test triangle
			float3 v0_i = float3(3.f, -3.f, 3.f);
			float3 v1_i = float3(3.f, 3.f, 3.f);
			float3 v2_i = float3(3.f, -3.f, -3.f);
			float3 c_i = (v0_i + v1_i + v2_i) / 3.f;
			float3 N_i = normalize(cross((v1_i - v0_i), (v2_i - v0_i)));
			arrPrimitive[0] = {
				{ v0_i.x, v0_i.y, v0_i.z, c_i.x },
				{ v1_i.x, v1_i.y, v1_i.z, c_i.y },
				{ v2_i.x, v2_i.y, v2_i.z, c_i.z },
				N_i.x,
				N_i.y,
				N_i.z,
				0,
				0
			};
			// Second test triangle
			float3 v0_j = float3(3.f, 3.f, -3.f);
			float3 v1_j = float3(3.f, -3.f, -3.f);
			float3 v2_j = float3(3.f, 3.f, 3.f);
			float3 c_j = (v0_i + v1_i + v2_i) / 3.f;
			float3 N_j = normalize(cross((v1_i - v0_i), (v2_i - v0_i)));
			arrPrimitive[1] = {
				{v0_j.x, v0_j.y, v0_j.z, c_j.x},
				{v1_j.x, v1_j.y, v1_j.z, c_j.y},
				{v2_j.x, v2_j.y, v2_j.z, c_j.z},
				N_j.x,
				N_j.y,
				N_j.z,
				1,
				0
			};

			float3 v0_d = float3(6.f, 3.f, -3.f);
			float3 v1_d = float3(6.f, -3.f, -3.f);
			float3 v2_d = float3(6.f, 3.f, 3.f);
			float3 c_d = (v0_d + v1_d + v2_d) / 3.f;
			float3 N_d = normalize(cross((v1_d - v0_d), (v2_d - v0_d)));

			arrPrimitive[2] = {
					{v0_d.x, v0_d.y, v0_d.z, c_d.x},
					{v1_d.x, v1_d.y, v1_d.z, c_d.y},
					{v2_d.x, v2_d.y, v2_d.z, c_d.z},
					N_d.x,
					N_d.y,
					N_d.z,
					2,
					1
			};

			float3 v0_3 = float3(3.f, 3.f, -3.f);
			float3 v1_3 = float3(3.f, -3.f, -3.f);
			float3 v2_3 = float3(3.f, 3.f, 3.f);
			float3 c_3 = (v0_3 + v1_3 + v2_3) / 3.f;
			float3 N_3 = normalize(cross((v1_3 - v0_3), (v2_3 - v0_3)));
			arrPrimitive[3] = {
				{v0_3.x, v0_3.y, v0_3.z, c_3.x},
				{v1_3.x, v1_3.y, v1_3.z, c_3.y},
				{v2_3.x, v2_3.y, v2_3.z, c_3.z},
				N_3.x,
				N_3.y,
				N_3.z,
				1,
				2
			};

			for (int i = NumberOfElements - 1; i >= 0; i--) { ; triColors[i] = { 1.f, 1.f, 0.f, 1.f }; arrPrimitiveIdx[i] = i; }
			//triColors[0] = { 1.f, 0.f, 0.f, 1.f };
			//triColors[1] = { 0.f, 1.f, 0.f, 1.f };
			//triColors[2] = { 0.5f, 1.f, 0.5f, 1.f };

			BuildBVH();
		}


		void BuildBVH()
		{
			// assign all triangles to root node
			BVH_GPU* root = &bvhNode[rootNodeIdx];
			root->leftFirst = 0;
			root->triCount = NumberOfElements;
			UpdateNodeBounds(rootNodeIdx);
			// subdivide recursively
			Subdivide(rootNodeIdx);
		}

		void Subdivide(cl_uint nodeIdx)
		{
			// terminate recursion
			BVH_GPU* node = &bvhNode[nodeIdx];
			if (node->triCount <= 1) return;
			// determine split axis and position
			float4 extent = node->aabbMax - node->aabbMin;
			int axis = 0;
			if (extent.y > extent.x) axis = 1;
			if (extent.z > extent[axis]) axis = 2;
			float splitPos = node->aabbMin[axis] + extent[axis] * 0.5f;

			int i = node->leftFirst;
			int j = i + node->triCount - 1;
			while (i <= j)
			{
				if (getCentroid(arrPrimitive[arrPrimitiveIdx[i]])[axis] < splitPos)
					i++;
				else
					swap(arrPrimitiveIdx[i], arrPrimitiveIdx[j--]);
			}

			int leftCount = i - node->leftFirst;
			if (leftCount == 0 || leftCount == node->triCount) return;

			int leftChildIdx = nodesUsed++;
			int rightChildIdx = nodesUsed++;
			bvhNode[leftChildIdx].leftFirst = node->leftFirst;
			bvhNode[leftChildIdx].triCount = leftCount;
			bvhNode[rightChildIdx].leftFirst = i;
			bvhNode[rightChildIdx].triCount = node->triCount - leftCount;
			node->leftFirst = leftChildIdx;
			node->triCount = 0;
			UpdateNodeBounds(leftChildIdx);
			UpdateNodeBounds(rightChildIdx);
			// recurse
			Subdivide(leftChildIdx);
			Subdivide(rightChildIdx);
		}
		void UpdateNodeBounds(cl_uint nodeIdx)
		{
			BVH_GPU* node = &bvhNode[nodeIdx];

			node->aabbMin = float4(1e30f);
			node->aabbMax = float4(-1e30f);
			for (uint first = node->leftFirst, i = 0; i < node->triCount; i++)
			{
				//createAABB(&arrPrimitive[*arrPrimitiveIdx[first + i]]);
				pair<float4, float4> aabb_minMax = createAABB(&arrPrimitive[arrPrimitiveIdx[first + i]]);
				node->aabbMin = fminf(node->aabbMin, aabb_minMax.first);
				node->aabbMax = fmaxf(node->aabbMax, aabb_minMax.second);
			}
		}

		void traverse(Ray& ray)
		{
			int stack_size = 0;
			BVH_GPU* stack[10];									// stack for traversing the tree
			stack[stack_size++] = &bvhNode[0];							// init with the first bvhnode

			while (stack_size > 0 )											// iterate over all values -> would run at most SCRWIDTH2*SCRHEIGHT2
			{
				BVH_GPU node =*stack[--stack_size];						// check if init is currect

				float tx1 = (node.aabbMin.x - ray.O.x) / ray.D.x, tx2 = (node.aabbMax.x - ray.O.x) / ray.D.x;
				float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
				float ty1 = (node.aabbMin.y - ray.O.y) / ray.D.y, ty2 = (node.aabbMax.y - ray.O.y) / ray.D.y;
				tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
				float tz1 = (node.aabbMin.z - ray.O.z) / ray.D.z, tz2 = (node.aabbMax.z - ray.O.z) / ray.D.z;
				tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));


				if (tmax >= tmin && tmin < ray.I.t && tmax > 0) {}
				else {
					int i = 0;

					if (node.triCount > 0)										// if root
					{

						for (uint i = 0; i < node.triCount; i++)
						{
							cout << "In" << endl;
							//if (triangles[primitives_idx[node.leftFirst + i]].objId == 0) IntersectTri(ray, &triangles[primitives_idx[node.leftFirst + i]]);
							//else { IntersectSphere(ray, &triangles[primitives_idx[node.leftFirst + i]], 2); }
						}
					}
					else												// add to the stack 
					{
						stack[stack_size++] = &bvhNode[arrPrimitiveIdx[node.leftFirst]];
						stack[stack_size++] = &bvhNode[arrPrimitiveIdx[node.leftFirst + 1]];
					}
				}
			}
		}
	};


};
