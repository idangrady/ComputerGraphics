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

	__declspec(align(64)) struct MaterialGPU {
		cl_float4 albedoSpecularity;
		cl_float4 absorption;
		cl_bool isEmissive;
		cl_int textureId;
	};

	class SceneGPU {
	public:
		SceneGPU() {
			tris = new TriGPU[4];
			mats = new MaterialGPU[4];
			// First test triangle
			float3 v0_i = float3(3.f, -3.f, 3.f);
			float3 v1_i = float3(3.f, 3.f, 3.f);
			float3 v2_i = float3(3.f, -3.f, -3.f);
			float3 c_i = (v0_i + v1_i + v2_i) / 3.f;
			float3 N_i = normalize(cross((v1_i - v0_i), (v2_i - v0_i)));
			tris[0] = {
				{ v0_i.x, v0_i.y, v0_i.z, c_i.x },
				{ v1_i.x, v1_i.y, v1_i.z, c_i.y },
				{ v2_i.x, v2_i.y, v2_i.z, c_i.z },
				N_i.x,
				N_i.y,
				N_i.z,
				0,
			};
			// Second test triangle
			float3 v0_j = float3(3.f, 3.f, -3.f);
			float3 v1_j = float3(3.f, -3.f, -3.f);
			float3 v2_j = float3(3.f, 3.f, 3.f);
			float3 c_j = (v0_j + v1_j + v2_j) / 3.f;
			float3 N_j = normalize(cross((v1_j - v0_j), (v2_j - v0_j)));
			tris[1] = {
				{v0_j.x, v0_j.y, v0_j.z, c_j.x},
				{v1_j.x, v1_j.y, v1_j.z, c_j.y},
				{v2_j.x, v2_j.y, v2_j.z, c_j.z},
				N_j.x,
				N_j.y,
				N_j.z,
				1,
			};
			// Lamp in the air
			float3 v0_k = float3(1.5f, 2.95f, 1.5f);
			float3 v1_k = float3(-1.5f, 2.95f, 1.5f);
			float3 v2_k = float3(1.5f, 2.95f, -1.5f);
			float3 c_k = (v0_k + v1_k + v2_k) / 3.f;
			float3 N_k = normalize(cross((v1_k - v0_k), (v2_k - v0_k)));
			tris[2] = {
				{v0_k.x, v0_k.y, v0_k.z, c_k.x},
				{v1_k.x, v1_k.y, v1_k.z, c_k.y},
				{v2_k.x, v2_k.y, v2_k.z, c_k.z},
				N_k.x,
				N_k.y,
				N_k.z,
				2,
			};
			float3 v0_l = float3(-1.5f, 2.95f, -1.5f);
			float3 v1_l = float3(1.5f, 2.95f, -1.5f);
			float3 v2_l = float3(-1.5f, 2.95f, 1.5f);
			float3 c_l = (v0_l + v1_l + v2_l) / 3.f;
			float3 N_l = normalize(cross((v1_l - v0_l), (v2_l - v0_l)));
			tris[3] = {
				{v0_l.x, v0_l.y, v0_l.z, c_l.x},
				{v1_l.x, v1_l.y, v1_l.z, c_l.y},
				{v2_l.x, v2_l.y, v2_l.z, c_l.z},
				N_l.x,
				N_l.y,
				N_l.z,
				3,
			};
			MaterialGPU mat0 = {
				{1, 0, 0, 0},
				{0, 0, 0, 0},
				false,
				-1,
			};
			MaterialGPU mat1 = {
				{0, 1, 0, 0},
				{0, 0, 0, 0},
				false,
				-1,
			};
			MaterialGPU mat2 = {
				{24, 24, 24, 0},
				{0, 0, 0, 0},
				true,
				-1,
			};
			mats[0] = mat0;
			mats[1] = mat1;
			mats[2] = mat2;
			mats[3] = mat2;
		}
		TriGPU* tris;
		MaterialGPU* mats;
		int tri_count = 4;
	};
}