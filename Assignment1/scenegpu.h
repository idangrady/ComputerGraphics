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
	__declspec(align(64)) struct TriGPU {
		cl_float4 vertex0, vertex1, vertex2;
		cl_float N0, N1, N2;
		cl_uint id;
	};

	__declspec(align(64)) struct MaterialGPU {
		cl_float4 albedoSpecularity;
		cl_float4 absorption;
		cl_bool isEmissive;
		cl_uint medium;
	};

	class SceneGPU {
	public:
		SceneGPU() {
			tris = new TriGPU[tri_count];
			triExs = new TriExGPU[tri_count];
			mats = new MaterialGPU[mat_count];
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
				{0.2, 0.8, 0.8, 0},
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
				id,
			};
			triExs[id] = triEx;
			tris[id++] = tri;
		}
		TriExGPU* triExs;
		TriGPU* tris;
		MaterialGPU* mats;
		int tri_count = 18;
		int mat_count = 5;
	};
}