#pragma 
#include <../lib/stb_image.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <config.h>

namespace Tmpl8 {
	__declspec(align(64)) struct TriExGPU {
		cl_float4 N;
		cl_float2 uv0, uv1, uv2;
		cl_int matId;
		cl_int textureId;
	};
	__declspec(align(64)) struct TriGPU {
		cl_float4 vertex0, vertex1, vertex2; // + centroid hidden in the last floats of the float4
	};

	__declspec(align(64)) struct MaterialGPU {
		cl_float4 albedoSpecularity;
		cl_float4 absorption;
		cl_bool isEmissive;
		cl_uint medium;
	};

	__declspec(align(8)) struct TextureData {
		cl_int width;
		cl_int height;
	};

	class MeshGPU {
	public:
		MeshGPU(uint objId, mat4 transform = mat4::Identity()) {
			objIdx = objId;
			M = transform;
			invM = transform.FastInvertedTransformNoScale();
		};

		void Translate(float3 d, vector<TriGPU>* tris) {
			for (int i = index; i < index + triangleCount; i++) {
				(*tris)[i].vertex0.y += d.y;
				(*tris)[i].vertex0.z += d.z;
				(*tris)[i].vertex0.w += d.x;

				(*tris)[i].vertex1.x += d.x;
				(*tris)[i].vertex1.y += d.y;
				(*tris)[i].vertex1.z += d.z;
				(*tris)[i].vertex1.w += d.y;

				(*tris)[i].vertex2.x += d.x;
				(*tris)[i].vertex2.y += d.y;
				(*tris)[i].vertex2.z += d.z;
				(*tris)[i].vertex2.w += d.z;
			}
		}

		void Scale(float3 d, vector<TriGPU>* tris) {
			for (int i = index; i < index + triangleCount; i++) {
				(*tris)[i].vertex0.x *= d.x;
				(*tris)[i].vertex0.y *= d.y;
				(*tris)[i].vertex0.z *= d.z;
				(*tris)[i].vertex0.w *= d.x;

				(*tris)[i].vertex1.x *= d.x;
				(*tris)[i].vertex1.y *= d.y;
				(*tris)[i].vertex1.z *= d.z;
				(*tris)[i].vertex1.w *= d.y;

				(*tris)[i].vertex2.x *= d.x;
				(*tris)[i].vertex2.y *= d.y;
				(*tris)[i].vertex2.z *= d.z;
				(*tris)[i].vertex2.w *= d.z;
			}
		}

		void MoveToPlane(float height, vector<TriGPU>* tris) {
			float lowest = FLT_MAX;
			for (int i = index; i < index + triangleCount; i++) {
				if ((*tris)[i].vertex0.y < lowest) lowest = (*tris)[i].vertex0.y;
				if ((*tris)[i].vertex1.y < lowest) lowest = (*tris)[i].vertex1.y;
				if ((*tris)[i].vertex2.y < lowest) lowest = (*tris)[i].vertex2.y;
			}
			float diff = height - lowest;
			Translate(float3(0, diff + 0.001f, 0), tris);
		}
		int index; // Index of mesh triangles in the triangle list
		int triangleCount; // Amount of triangles in mesh
		int textureIndex; // Index of texture
		int width, height; // Texture width/height
		uint objIdx;
		mat4 M, invM;
		MaterialGPU* material;
	private:
	};

	class SceneGPU {
	public:
		SceneGPU() {
			//tris = new TriGPU[tri_count];
			//triExs = new TriExGPU[tri_count];
			//mats = new MaterialGPU[mat_count];

			MakeSimpleScene();
			//MakeComplexScene();

			cout << "Scene constructed with " << tris.size() << " triangles.\n";

			// Load skybox
			LoadSkyBox("assets/clarens_midday_4k.hdr");
		}

		void MakeComplexScene() {
			MaterialGPU default = {
				{1, 1, 1, 0},
				{0, 0, 0, 0},
				false,
				0,
			};
			mats.push_back(default);
			loadModel("assets/chessboard/chessboard.obj");
			meshPool[0]->MoveToPlane(-128, &tris);
		}

		void MakeSimpleScene() {
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
			MakeTriangle(float3(-3.f, -3.f, 3.f), float3(-3.f, -3.f, -3.f), float3(-3.f, 3.f, 3.f), 0);
			// Left wall 2
			MakeTriangle(float3(-3.f, 3.f, -3.f), float3(-3.f, 3.f, 3.f), float3(-3.f, -3.f, -3.f), 0);
			// Right wall 1
			MakeTriangle(float3(3.f, -3.f, 3.f), float3(3.f, 3.f, 3.f), float3(3.f, -3.f, -3.f), 1);
			// Right wall 2
			MakeTriangle(float3(3.f, 3.f, -3.f), float3(3.f, -3.f, -3.f), float3(3.f, 3.f, 3.f), 1);
			// Ceiling 1
			MakeTriangle(float3(3.f, 3.f, 3.f), float3(-3.f, 3.f, 3.f), float3(3.f, 3.f, -3.f), 2);
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
			float3 N = normalize(cross((v1 - v0), (v2 - v0)));
			float3 C = (v0 + v1 + v2) / 3.f;
			TriExGPU triEx = {
				{N.x, N.y, N.z, 0.f},
				{0, 0},
				{0, 0},
				{0, 0},
				mat,
				-1,
			};
			TriGPU tri = {
				{v0.x, v0.y, v0.z, C.x},
				{v1.x, v1.y, v1.z, C.y},
				{v2.x, v2.y, v2.z, C.z},
			};
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

		// Mesh loading using tutorial from learnopengl.com
		void loadModel(const char* file) {
			Assimp::Importer importer;
			const aiScene* scene = importer.ReadFile(file, aiProcess_Triangulate);
			if (!scene || scene->mFlags * AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
				cout << "Error loading Mesh: " << importer.GetErrorString() << endl;
				return;
			}
			processNode(scene->mRootNode, scene, file);
			cout << "Model loaded. Mesh pool size: " << meshPool.size() << endl;
			cout << "Texture pool size: " << textures.size() * 4 << " bytes.\n";
		}

		void processNode(aiNode* node, const aiScene* scene, const char* path) {
			// process parent node
			for (int i = 0; i < node->mNumMeshes; i++) {
				int id = meshPool.size(); //TODO: REMOVE THE START AT 2
				aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
				meshPool.push_back(makeMesh(mesh, scene, path, id, 0));
			}
			// go through children nodes
			for (uint i = 0; i < node->mNumChildren; i++) {
				processNode(node->mChildren[i], scene, path);
			}
		}

		MeshGPU* makeMesh(aiMesh* mesh, const aiScene* scene, string const& path, uint objId, int matId) {
			MeshGPU* m = new MeshGPU(objId);
			if (matId >= 0) m->material = &(mats[matId]);
			else m->material = NULL;
			m->index = tris.size();
			m->triangleCount = mesh->mNumFaces;
			m->textureIndex = -1;
			tris.reserve(tris.size() + mesh->mNumFaces);
			triExs.reserve(triExs.size() + mesh->mNumFaces);
			string directory = path.substr(0, path.find_last_of('/'));


			// Load texture
			aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
			// Diffuse texture
			if (material->GetTextureCount(aiTextureType_DIFFUSE) > 0) {
				aiString str;
				material->GetTexture(aiTextureType_DIFFUSE, 0, &str);
				string filename = string(str.C_Str());
				filename = directory + "/" + filename;
				Surface* textureMap = new Surface(filename.c_str());
				m->textureIndex = textureIndices.size();
				m->width = textureMap->width;
				m->height = textureMap->height;
				TextureData td;
				td.width = m->width;
				td.height = m->height;
				textureIndices.push_back(m->textureIndex);
				textureData.push_back(td);

				// Copy texture to vector
				int textureLength = textureMap->width * textureMap->height;
				textures.reserve(textures.size() + textureLength);
				std::copy(textureMap->pixels, textureMap->pixels + textureLength, back_inserter(textures));

				// Delete the allocated Surface (Look at me go, managing memory! I'm good at C++!)
				delete(textureMap);
			}

			// Load potential normal map
			if (material->GetTextureCount(aiTextureType_HEIGHT) > 0) {
				aiString str;
				material->GetTexture(aiTextureType_HEIGHT, 0, &str);
				string filename = string(str.C_Str());
				filename = directory + "/" + filename;
				Surface* normalMap = new Surface(filename.c_str());
				// NOT YET IMPLEMENTED
				//processing normal map here
				delete(normalMap);
			}

			// We load things face by face instead of vertex by vertex. Path tracers like accessing triangles compactly, I think.
			// Potentially can be changed into a list of vertices and a list of face indices, to save memory.
			for (uint i = 0; i < mesh->mNumFaces; i++) {
				TriGPU tri_i = TriGPU();
				TriExGPU triEx_i = TriExGPU();
				aiFace face = mesh->mFaces[i];
				aiVector3D c = (mesh->mVertices[face.mIndices[0]] + mesh->mVertices[face.mIndices[1]] + mesh->mVertices[face.mIndices[2]]) / 3.0f;
				// vertex0
				tri_i.vertex0.x = mesh->mVertices[face.mIndices[0]].x;
				tri_i.vertex0.y = mesh->mVertices[face.mIndices[0]].y;
				tri_i.vertex0.z = mesh->mVertices[face.mIndices[0]].z;
				tri_i.vertex0.w = c.x;
				// vertex1
				tri_i.vertex1.x = mesh->mVertices[face.mIndices[1]].x;
				tri_i.vertex1.y = mesh->mVertices[face.mIndices[1]].y;
				tri_i.vertex1.z = mesh->mVertices[face.mIndices[1]].z;
				tri_i.vertex1.w = c.y;
				// vertex2
				tri_i.vertex2.x = mesh->mVertices[face.mIndices[2]].x;
				tri_i.vertex2.y = mesh->mVertices[face.mIndices[2]].y;
				tri_i.vertex2.z = mesh->mVertices[face.mIndices[2]].z;
				tri_i.vertex2.w = c.z;

				// Load the normals
				// Since trying to use vertex normals is not really working, just calculate face normal
				float3 v0 = { tri_i.vertex0.x, tri_i.vertex0.y, tri_i.vertex0.z };
				float3 v1 = { tri_i.vertex1.x, tri_i.vertex1.y, tri_i.vertex1.z };
				float3 v2 = { tri_i.vertex2.x, tri_i.vertex2.y, tri_i.vertex2.z };
				float3 N = normalize(cross((v1 - v0), (v2 - v0)));
				triEx_i.N = { N.x, N.y, N.z, 0 };
				//// Load the normals
				//if (mesh->HasNormals()) {
				//}
				//else {
				//	float3 v0 = { tri_i.vertex0.x, tri_i.vertex0.y, tri_i.vertex0.z };
				//	float3 v1 = { tri_i.vertex1.x, tri_i.vertex1.y, tri_i.vertex1.z };
				//	float3 v2 = { tri_i.vertex2.x, tri_i.vertex2.y, tri_i.vertex2.z };
				//	float3 N = normalize(cross((v1 - v0), (v2 - v0)));
				//	triEx_i.N = { N.x, N.y, N.z, 0 };
				//}

				// Load Texture coordinates
				if (mesh->mTextureCoords[0]) {
					triEx_i.uv0 = { mesh->mTextureCoords[0][face.mIndices[0]].x, mesh->mTextureCoords[0][face.mIndices[0]].y };
					triEx_i.uv1 = { mesh->mTextureCoords[0][face.mIndices[1]].x, mesh->mTextureCoords[0][face.mIndices[1]].y };
					triEx_i.uv2 = { mesh->mTextureCoords[0][face.mIndices[2]].x, mesh->mTextureCoords[0][face.mIndices[2]].y };
				}
				else {
					triEx_i.uv0 = { 0, 0 };
					triEx_i.uv1 = { 0, 0 };
					triEx_i.uv2 = { 0, 0 };
				}
				triEx_i.matId = matId;
				triEx_i.textureId = m->textureIndex;
				//PrintTri(tri_i);
				//PrintTriEx(triEx_i);
				tris.push_back(tri_i);
				triExs.push_back(triEx_i);
			}
			return m;
		}

		void PrintTri(TriGPU tri) {
			cout << "Triangle Vertices:\n";
			cout << "v0\t{" << tri.vertex0.x << ", " << tri.vertex0.y << ", " << tri.vertex0.z << "}\n";
			cout << "v1\t{" << tri.vertex1.x << ", " << tri.vertex1.y << ", " << tri.vertex1.z << "}\n";
			cout << "v2\t{" << tri.vertex2.x << ", " << tri.vertex2.y << ", " << tri.vertex2.z << "}\n";
			cout << "Triangle Centroid:\t{" << tri.vertex0.w << ", " << tri.vertex1.w << ", " << tri.vertex2.w << "}\n";
		}

		void PrintTriEx(TriExGPU tri) {
			cout << "Triangle Normal:\t{" << tri.N.x << ", " << tri.N.y << ", " << tri.N.z << "}\n";
			cout << "Triangle uv:\n";
			cout << "uv0:\t{" << tri.uv0.x << ", " << tri.uv0.y << "}\n";
			cout << "uv1:\t{" << tri.uv1.x << ", " << tri.uv1.y << "}\n";
			cout << "uv2:\t{" << tri.uv2.x << ", " << tri.uv2.y << "}\n";
			cout << "Mat ID:\t" << tri.matId << "\n";
			cout << "Tex ID:\t" << tri.textureId << "\n";
		}

		vector<TriExGPU> triExs;
		vector<TriGPU> tris;
		vector<MaterialGPU> mats;
		vector<uint> textures;
		vector<int> textureIndices;
		vector<TextureData> textureData;
		vector<MeshGPU*> meshPool;
		//TriExGPU* triExs;
		//TriGPU* tris;
		//MaterialGPU* mats;
		//int tri_count = 18;
		//int mat_count = 5;

		//Skybox
		float4* skybox;
		int width, height, nrChannels;
		bool loadedSkybox;
	};

}