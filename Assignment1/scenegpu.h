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

	class SceneGPU {
	public:
		SceneGPU() {
			tris = new TriGPU[2];
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
			float3 c_j = (v0_i + v1_i + v2_i) / 3.f;
			float3 N_j = normalize(cross((v1_i - v0_i), (v2_i - v0_i)));
			tris[1] = {
				{v0_j.x, v0_j.y, v0_j.z, c_j.x},
				{v1_j.x, v1_j.y, v1_j.z, c_j.y},
				{v2_j.x, v2_j.y, v2_j.z, c_j.z},
				N_j.x,
				N_j.y,
				N_j.z,
				1,
			};
			triColors = new cl_float4[2];
			triColors[0] = { 1.f, 0.f, 0.f, 1.f };
			triColors[1] = { 0.f, 1.f, 0.f, 1.f };
		}
		TriGPU* tris;
		cl_float4* triColors;
	};
}