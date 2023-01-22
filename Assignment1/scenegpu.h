#pragma 
#include <../lib/stb_image.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <config.h>


static const int numTrigs = 18;

namespace Tmpl8 {
	__declspec(align(8)) struct TriExGPU {
		cl_int matId;
		cl_int textureId;
	};
	__declspec(align(64)) struct TriGPU {
		cl_float4 vertex0, vertex1, vertex2;
		cl_float N0, N1, N2;
		cl_uint id;
		cl_uint ObjId;

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

	pair<float4, float4>  createAABB(TriGPU &primitive)
	{
		float4 aabbMin = float4(1e30f);
		float4 aabbMax = float4(-1e30f);
		TriGPU& tri = primitive;
		if (tri.ObjId == 0) {
			// Triangle
			float4 ver_1 = float4(primitive.vertex0.x, primitive.vertex0.y, primitive.vertex0.z, 1);
			float4 ver_2 = float4(primitive.vertex1.x, primitive.vertex1.y, primitive.vertex1.z, 1);
			float4 ver_3 = float4(primitive.vertex2.x, primitive.vertex2.y, primitive.vertex2.z, 1);

			aabbMin = fminf(aabbMin, ver_1);
			aabbMin = fminf(aabbMin, ver_2);
			aabbMin = fminf(aabbMin, ver_3);
			aabbMax = fmaxf(aabbMax, ver_1);
			aabbMax = fmaxf(aabbMax, ver_2);
			aabbMax = fmaxf(aabbMax, ver_3);
		}
		//else if (primitive->ObjId == 1) {
		//	// sphere
		//	auto radius = primitive->vertex1.x;
		//	aabbMin = getCentroid(primitive) - float4(radius, radius, radius, radius);
		//	aabbMax = getCentroid(primitive) + float4(radius, radius, radius, radius);
		//}
		//else if (primitive->objId == 2) {
		//	// Cube
		//}
		return std::make_pair(aabbMin, aabbMax);
	}


	class simpleBVHGPU
	{
	public:
		simpleBVHGPU() {
			this->arrPrimitive = new TriGPU[numTrigs];
		};

		void add_primitive(vector<TriGPU>& tris) {
			memcpy(this->arrPrimitive, tris.data(), tris.size() * sizeof(TriGPU));
			this->count = tris.size();

			for (int j = 0; j < 2 * count - 1; j++) {
				bvhNode[j] = new BVH_GPU();
				if (j < count) { arrPrimitiveIdx[j] = j; }
			}


			//for (int i = 0; i < 18; i++)
			//{
			//	TriGPU current = arrPrimitive[i];
			//	cout << current.ObjId << endl;
			//	// or
			//	cout << (arrPrimitive + i)->ObjId << endl;
			//	//or
			//	cout << (*(arrPrimitive + i)).ObjId << endl;
			//}
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
				TriGPU &ssass = arrPrimitive[curr_val];

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
				auto curNode = arrPrimitive[arrPrimitiveIdx[i]];
				float3 centroid = float3(curNode.vertex0.w, curNode.vertex1.w, curNode.vertex2.w);
				if (centroid[axis] < splitPos)
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
		int count;
		BVH_GPU* bvhNode[2 * numTrigs - 1];
		uint rootNodeIdx = 0, nodesUsed = 1;
		TriGPU* arrPrimitive;
		int arrPrimitiveIdx[numTrigs];
		int ss_ = 0;
	};



	class SceneGPU {
	public:
		SceneGPU() {
			//tris = new TriGPU[tri_count];
			//triExs = new TriExGPU[tri_count];
			//mats = new MaterialGPU[mat_count];
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
				{24, 24, 24, 0},
				{0, 0, 0, 0},
				true,
				0,
			};
			MaterialGPU redglass = {
				{0, 0, 1, 0},
				{0, 0, 0, 0},
				false,
				1,
			};
			mats.push_back(red);
			mats.push_back(green);
			mats.push_back(white);
			mats.push_back(lamp);
			mats.push_back(redglass);
			/*mats[0] = red;
			mats[1] = green;
			mats[2] = white;
			mats[3] = lamp;
			mats[4] = redglass;*/
			// Left wall 1
			MakeTriangle(float3(-300001.f, -3.f, 3.f), float3(-30001.f, -30001.f, -3.f), float3(-300001.f, 3.f, 3.f), 0);
			// Left wall 2
			MakeTriangle(float3(-3.f, 3.f, -3.f), float3(-3.f, 3006.f, 30008.f), float3(-3008.f, -3.f, -29999.f), 0);
			// Right wall 1
			MakeTriangle(float3(3.0005f, -3.f, 30002.f), float3(3.f, 3.f, 3.f), float3(3.f, -3.f, -3.f), 1);
			// Right wall 2
			MakeTriangle(float3(3.f, 30007.f, -30007.f), float3(3.f, -3.f, -3.f), float3(3.f, 3.f, 3.f), 1);
			// Ceiling 1
			MakeTriangle(float3(3.f, 3.f, 3.f), float3(-3.f, 3.f, 3.f), float3(3.f, 3.f, -3.f), 2);
			// Ceiling 2
			MakeTriangle(float3(-3.f, 30007.f, -3.f), float3(3.f, 3.f, -3.f), float3(-3.f, 3.f, 3.f), 2);
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

			// Load skybox
			LoadSkyBox("assets/clarens_midday_4k.hdr");

			BVHclass = new simpleBVHGPU();

			BVHclass->add_primitive(tris);
			BVHclass->BuildBVH();
		}
		void MakeTriangle(float3 v0, float3 v1, float3 v2, int mat) {
			static int id = 0;
			float3 N = normalize(cross((v1 - v0), (v2 - v0)));
			float3 C = (v0 + v1 + v2) / 3.f;
			TriExGPU triEx = {
				mat,
				-1,
			};
			TriGPU tri = {
				{v0.x, v0.y, v0.z, C.x},
				{v1.x, v1.y, v1.z, C.y},
				{v2.x, v2.y, v2.z, C.z},
				N.x,
				N.y,
				N.z,
				id++,
				0
			};
			trigCount++;
			triExs.push_back(triEx);
			tris.push_back(tri);
		}
		void LoadSkyBox(const char* filename) {
			float* data = stbi_loadf(filename, &width, &height, &nrChannels, 0);
			if (data) {
				skybox = new float4[width * height];
				if (nrChannels == 3) //rgb
				{
					for (int i = 0; i < width * height; i++) {
						skybox[i] = float4(data[3 * i], data[3 * i + 1], data[3 * i + 2], 1);
					}
				}
				else if (nrChannels == 4) //rgba
				{
					for (int i = 0; i < width * height; i++) {
						skybox[i] = float4(data[4 * i], data[4 * i + 1], data[4 * i + 2], data[4 * i + 3]);
					}
				}
				else {
					throw exception("Skybox image needs either 3 or 4 channels.");
				}
				cout << "Skybox loaded. W: " << width << " H: " << height << "\n";
				loadedSkybox = true;
				stbi_image_free(data);
			}
			else {
				cout << stbi_failure_reason() << endl;
				throw exception("Failed to load Skybox.");
			}
		}

		vector<TriExGPU> triExs;
		vector<TriGPU> tris;
		vector<MaterialGPU> mats;
		
		// BVH
		simpleBVHGPU* BVHclass;
		int trigCount =0;

		//Skybox
		float4* skybox;
		int width, height, nrChannels;
		bool loadedSkybox;
	};
}