#pragma 
#include <../lib/stb_image.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <config.h>

namespace Tmpl8 {
	__declspec(align(8)) struct TriExGPU {
		cl_int matId;
		cl_int textureId;
	};
	__declspec(align(64)) struct Primitive_GPU {
		cl_float4 vertex0, vertex1, vertex2;
		cl_float N0, N1, N2;
		cl_uint id;
		cl_uint objId;
	};

	__declspec(align(64)) struct MaterialGPU {
		cl_float4 albedoSpecularity;
		cl_float4 absorption;
		cl_bool isEmissive;
		cl_uint medium;
	};

	__declspec(align(64)) struct BVH_GPU
	{
		float4 aabbMin, aabbMax;			// boundary
		uint leftFirst, triCount;			// count and start
	};

	float4 getCentroid(Primitive_GPU *primitive)
	{
		if (primitive->objId = 0)
		{
			auto x = (primitive->vertex0.x + primitive->vertex1.x + primitive->vertex2.x) / 3;
			auto y = (primitive->vertex0.y + primitive->vertex1.y + primitive->vertex2.y) / 3;
			auto z = (primitive->vertex0.z + primitive->vertex1.z + primitive->vertex2.z) / 3;
			return float4(x, y, z,1);
		}
		else { return float4(primitive->vertex0.x, primitive->vertex0.y, primitive->vertex0.z, 1); }
	}

	pair<float4, float4>  createAABB(Primitive_GPU* primitive)
	{
		float4 aabbMin = float4(1e30f);
		float4 aabbMax = float4(-1e30f);
		if (primitive->objId == 0) {
			// Triangle
			float4 ver_1 = float4(primitive->vertex0.x, primitive->vertex0.y, primitive->vertex0.z,1);
			float4 ver_2 = float4(primitive->vertex1.x, primitive->vertex1.y, primitive->vertex1.z, 1);
			float4 ver_3 = float4(primitive->vertex2.x, primitive->vertex2.y, primitive->vertex2.z,1);

			aabbMin = fminf(aabbMin, ver_1);
			aabbMin = fminf(aabbMin, ver_2);
			aabbMin = fminf(aabbMin, ver_3);
			aabbMax = fmaxf(aabbMax, ver_1);
			aabbMax = fmaxf(aabbMax, ver_2);
			aabbMax = fmaxf(aabbMax, ver_3);
		}
		else if (primitive->objId == 1) {
			// sphere
			auto radius = primitive->vertex1.x;
			aabbMin = getCentroid(primitive) - float4(radius, radius, radius, radius);
			aabbMax = getCentroid(primitive) + float4(radius, radius, radius, radius);
		}
		else if (primitive->objId == 2) {
			// Cube
		}
		return std::make_pair(aabbMin, aabbMax);
	}


	class simpleBVHGPU 
	{
	public:
		simpleBVHGPU() = default;

		void add_primitive(Primitive_GPU* arrPrimitive, int count)
		{
			this->count = count;
			for (int i = count - 1; i >= 0; i--)
			{
				// assigning the poiters for the 
				this->arrPrimitive[i] = &arrPrimitive[i];
				arrPrimitiveIdx[i] = i;
			}
			for (int j = 0; j < 2 * count - 1; j++) {
				bvhNode[j] = new BVH_GPU();

			}
		}
		void BuildBVH()
		{
			// assign all triangles to root node
			BVH_GPU& root = *bvhNode[rootNodeIdx];
			root.leftFirst = 0;
			root.triCount = count;
			UpdateNodeBounds(rootNodeIdx);
			// subdivide recursively
			Subdivide(rootNodeIdx);
		}
		void UpdateNodeBounds(uint nodeIdx)
		{
			BVH_GPU& node = *bvhNode[nodeIdx];

			node.aabbMin = float4(1e30f);
			node.aabbMax = float4(-1e30f);
			for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
			{
				cout << "Current: " << arrPrimitiveIdx[first + i] << endl;
				auto curr_val = arrPrimitiveIdx[first + i];
				Primitive_GPU *ssass = arrPrimitive[curr_val];
				cout << ssass << endl;
				pair<float4, float4> aabb_minMax = createAABB(ssass);
				node.aabbMin = fminf(node.aabbMin, aabb_minMax.first);
				node.aabbMax = fmaxf(node.aabbMax, aabb_minMax.second);
			}
		}
		void Subdivide(uint nodeIdx)
		{
			// terminate recursion
			BVH_GPU& node = *bvhNode[nodeIdx];
			if (node.triCount <= 1) return;
			// determine split axis and position
			float3 extent = node.aabbMax - node.aabbMin;
			int axis = 0;
			if (extent.y > extent.x) axis = 1;
			if (extent.z > extent[axis]) axis = 2;
			float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
			// in-place partition
			int i = node.leftFirst;
			int j = i + node.triCount - 1;
			while (i <= j)
			{
				if (getCentroid(arrPrimitive[arrPrimitiveIdx[i]])[axis] < splitPos)
					i++;
				else
					swap(arrPrimitiveIdx[i], arrPrimitiveIdx[j--]);
			}
			// abort split if one of the sides is empty
			int leftCount = i - node.leftFirst;
			if (leftCount == 0 || leftCount == node.triCount) return;

			int leftChildIdx = nodesUsed++;
			int rightChildIdx = nodesUsed++;
			bvhNode[leftChildIdx]->leftFirst = node.leftFirst;
			bvhNode[leftChildIdx]->triCount = leftCount;
			bvhNode[rightChildIdx]->leftFirst = i;
			bvhNode[rightChildIdx]->triCount = node.triCount - leftCount;
			node.leftFirst = leftChildIdx;
			node.triCount = 0;
			UpdateNodeBounds(leftChildIdx);
			UpdateNodeBounds(rightChildIdx);
			// recurse
			Subdivide(leftChildIdx);
			Subdivide(rightChildIdx);
		}
		int count = 18;

		BVH_GPU* bvhNode[2*18-1];
		uint rootNodeIdx = 0, nodesUsed = 1;
		Primitive_GPU* arrPrimitive[18];
		int arrPrimitiveIdx[18];
		int ss_ = 0;
	};



	class SceneGPU {
	public:
		SceneGPU() {
			arrPrimitive = new Primitive_GPU[tri_count];
			triExs = new TriExGPU[tri_count];
			mats = new MaterialGPU[mat_count];

			bvhNode = new BVH_GPU[2 * tri_count - 1];
			arrPrimitiveIdx = new int[tri_count];

			MaterialGPU red = {
				{1, 0, 0, 0},
				{0, 0, 0, 0},
				false,
				0,
			};
			MaterialGPU green = {
				{0, 1, 0, 0},
				{0, 0, 0, 0},
				false,
				0,
			};
			MaterialGPU white = {
				{0.8f, 0.8f, 0.8f, 0.f},
				{0, 0, 0, 0},
				false,
				0,
			};
			MaterialGPU lamp = {
				{24, 24, 24 , 0},
				{0, 0, 0, 0},
				true,
				0,
			};
			MaterialGPU redglass = {
				{0, 0, 1, 0},
				{0, 5, 5, 0},
				false,
				1,
			};
			mats[0] = red;
			mats[1] = green;
			mats[2] = white;
			mats[3] = lamp;
			mats[4] = redglass;
			// Left wall 1
			MakeTriangle(float3(-3.f, -3.f, 3.f), float3(-3.f, -3.f, -3.f), float3(-3.f, 3.f, 3.f), 0);
			// Left wall 2
			MakeTriangle(float3(-3.f, 3.f, -3.f), float3(-3.f, 3.f, 3.f), float3(-3.f, -3.f, -3.f), 0);
			// Right wall 1
			MakeTriangle(float3(3.f, -3.f, 3.f), float3(3.f, 3.f, 3.f), float3(3.f, -3.f, -3.f), 1);
			// Right wall 2
			MakeTriangle(float3(3.f, 3.f, -3.f), float3(3.f, -3.f, -3.f), float3(3.f, 3.f, 3.f), 1);
			// Ceiling 1
			MakeTriangle(float3(3.f, 3.f, 3.f), float3(-3.f, 3.f, 3.f), float3(3.f, 3.f,- 3.f), 2);
			// Ceiling 2
			MakeTriangle(float3(-3.f, 3.f, -3.f), float3(3.f, 3.f, -3.f), float3(-3.f, 3.f, 3.f), 2);
			// Back wall 1
			MakeTriangle(float3(-3.f, -3.f, 3.f), float3(-3.f, 3.f, 3.f), float3(3.f, -3.f, 3.f), 2);
			// Back wall 2
			MakeTriangle(float3(3.f, 3.f, 3.f), float3(3.f, -3.f, 3.f), float3(-3.f, 3.f, 3.f), 2);
			// Floor 1
			MakeTriangle(float3(3.f, -3.f, 3.f), float3(3.f, -3.f, -3.f), float3(-3.f, -3.f, 3.f), 2);
			// Floor 2
			MakeTriangle(float3(-3.f, -3.f, -3.f), float3(-3.f, -3.f, 3.f), float3(3.f, -3.f, -3.f), 2);
			// Lamp in the air
			MakeTriangle(float3(1.5f, 2.95f, 1.5f), float3(-1.5f, 2.95f, 1.5f), float3(1.5f, 2.95f, -1.5f), 3);
			MakeTriangle(float3(-1.5f, 2.95f, -1.5f), float3(1.5f, 2.95f, -1.5f), float3(-1.5f, 2.95f, 1.5f), 3);
			// Pyramid floor
			MakeTriangle(float3(0.f, -2.5f, -1.f), float3(1.f, -2.5f, 0.f), float3(-1.f, -2.5f, 0.f), 4);
			MakeTriangle(float3(0.f, -2.5f, 1.f), float3(-1.f, -2.5f, 0.f), float3(1.f, -2.5f, 0.f), 4);
			// Pyramid walls
			MakeTriangle(float3(0.f, -1.f, 0.f), float3(0.f, -2.5f, -1.f), float3(1.f, -2.5f, 0.f), 4);
			MakeTriangle(float3(0.f, -1.f, 0.f), float3(1.f, -2.5f, 0.f), float3(0.f, -2.5f, 1.f), 4);
			MakeTriangle(float3(0.f, -1.f, 0.f), float3(0.f, -2.5f, 1.f), float3(-1.f, -2.5f, 0.f), 4);
			MakeTriangle(float3(0.f, -1.f, 0.f), float3(-1.f, -2.5f, 0.f), float3(0.f, -2.5f, -1.f), 4);

			//SimpleBVHGPUInstance = new simpleBVHGPU();
			//SimpleBVHGPUInstance->add_primitive(arrPrimitive, tri_count);		// initlize the BVH tree
			//SimpleBVHGPUInstance->BuildBVH();

			bvhNode = new BVH_GPU[18 * 2 - 1];
			arrPrimitiveIdx = new int[18];
			for (int i = 0; i < 2*tri_count-1; i++) 
			{
				bvhNode[i] = BVH_GPU();
				if (i < 18) arrPrimitiveIdx[i] = i;
			}

			BuildBVH();
		}
 
		void MakeTriangle(float3 v0, float3 v1, float3 v2, int mat, int obj = 0) {
			static int id = 0;
			float3 N = normalize(cross((v1 - v0), (v2 - v0)));
			float3 C = (v0 + v1 + v2) / 3.f;
			TriExGPU triEx = {
				mat,
				-1,
			};
			Primitive_GPU tri = {
				{v0.x, v0.y, v0.z, C.x},
				{v1.x, v1.y, v1.z, C.y},
				{v2.x, v2.y, v2.z, C.z},
				N.x,
				N.y,
				N.z,
				id,
				0
			};
			triExs[id] = triEx;
			arrPrimitive[id++] = tri;
		}

		void BuildBVH()
		{
			// assign all triangles to root node
			BVH_GPU& root = bvhNode[rootNodeIdx];
			root.leftFirst = 0;
			root.triCount = N_bvh;
			UpdateNodeBounds(rootNodeIdx);
			// subdivide recursively
			Subdivide(rootNodeIdx);
		}

		void UpdateNodeBounds(uint nodeIdx)
		{
			BVH_GPU& node = bvhNode[nodeIdx];

			node.aabbMin = float4(1e30f);
			node.aabbMax = float4(-1e30f);
			for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
			{
				cout << "Current: " << arrPrimitiveIdx[first + i] << endl;
				auto curr_val = arrPrimitiveIdx[first + i];
				Primitive_GPU& ssass = arrPrimitive[curr_val];
				cout << &ssass << endl;
				pair<float4, float4> aabb_minMax = createAABB(&ssass);
				node.aabbMin = fminf(node.aabbMin, aabb_minMax.first);
				node.aabbMax = fmaxf(node.aabbMax, aabb_minMax.second);
			}
		}


		void Subdivide(uint nodeIdx)
		{
			// terminate recursion
			BVH_GPU& node = bvhNode[nodeIdx];
			if (node.triCount <= 1) return;
			// determine split axis and position
			float3 extent = node.aabbMax - node.aabbMin;
			int axis = 0;
			if (extent.y > extent.x) axis = 1;
			if (extent.z > extent[axis]) axis = 2;
			float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
			// in-place partition
			int i = node.leftFirst;
			int j = i + node.triCount - 1;
			while (i <= j)
			{
				auto centroid = getCentroid(&arrPrimitive[arrPrimitiveIdx[i]]);
				if (getCentroid(&arrPrimitive[arrPrimitiveIdx[i]])[axis] < splitPos)
					i++;
				else
					swap(arrPrimitiveIdx[i], arrPrimitiveIdx[j--]);
			}
			// abort split if one of the sides is empty
			int leftCount = i - node.leftFirst;
			if (leftCount == 0 || leftCount == node.triCount) return;

			//------------------------------------------------------------------Start here ---Check this part please
			// create child nodes

			// create child nodes	
			int leftChildIdx = nodesUsed++;
			int rightChildIdx = nodesUsed++;
			bvhNode[leftChildIdx].leftFirst = node.leftFirst;
			bvhNode[leftChildIdx].triCount = leftCount;
			bvhNode[rightChildIdx].leftFirst = i;
			bvhNode[rightChildIdx].triCount = node.triCount - leftCount;
			node.leftFirst = leftChildIdx;
			node.triCount = 0;
			UpdateNodeBounds(leftChildIdx);
			UpdateNodeBounds(rightChildIdx);
			// recurse
			Subdivide(leftChildIdx);
			Subdivide(rightChildIdx);
		}




		TriExGPU* triExs;
		MaterialGPU* mats;

		//simpleBVHGPU* SimpleBVHGPUInstance;
		BVH_GPU* bvhNode;
		Primitive_GPU *arrPrimitive;
		int* arrPrimitiveIdx;

		int tri_count = 18;
		int mat_count = 5;
		int nodesUsed = 0;
		int rootNodeIdx = 0;
};
}