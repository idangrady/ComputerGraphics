#pragma once
#include <map>
#include <random>
#include <../lib/stb_image.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

// -----------------------------------------------------------
// scene.h
// Simple test scene for ray.I.tracing experiments. Goals:
// - Super-fast scene intersection
// - Easy interface: scene.FindNearest / IsOccluded
// - With normals and albedo: GetNormal / GetAlbedo
// - Area light source (animated), for light transport
// - Primitives can be hit from inside - for dielectrics
// - Can be extended with other primitives and/or a BVH
// - Optionally animated - for temporal experiments
// - Not everything is axis aligned - for cache experiments
// - Can be evaluated at arbitrary time - for motion blur
// - Has some high-frequency details - for filtering
// Some speed tricks that severely affect maintainability
// are enclosed in #ifdef SPEEDTRIX / #endif. Mind these
// if you plan to alter the scene in any way.
// -----------------------------------------------------------

#define SPEEDTRIX
#define PLANE_X(o,i) {if((t=-(ray.O.x+o)*ray.rD.x)<ray.I.t)ray.I.t=t,ray.objIdx=i;}
#define PLANE_Y(o,i) {if((t=-(ray.O.y+o)*ray.rD.y)<ray.I.t)ray.I.t=t,ray.objIdx=i;}
#define PLANE_Z(o,i) {if((t=-(ray.O.z+o)*ray.rD.z)<ray.I.t)ray.I.t=t,ray.objIdx=i;}

enum class Medium {
	Undefined = -1,
	Air = 0,
	Glass = 1,
};

struct Light {
	Light() = default;
	Light(float3 p, float i) : position(p), lumen(i) {}
	float3 position;
	float3 color;
	float lumen;
};

struct Material
{
	Material() = default;
	Material(float3 a, float s = 1, Medium m = Medium::Undefined) : albedo(a), specularity(s), mat_medium(m) {}
	float3 albedo = (0.9f, 0.9f, 0.9f); // color material
	float specularity;
	bool isLight = false;
	float3 absorption = (0, 0, 0);
	Medium mat_medium{ Medium::Undefined };
};

// intersection record, carefully tuned to be 16 bytes in size and equally carefully yoinked from jacco.ompf2.com
struct Intersection
{
	Intersection() = default;
	float t = 1e34f;		// intersection distance along ray
	float u, v;		// barycentric coordinates of the intersection
	uint instPrim = 0;	// instance index (12 bit) and primitive index (20 bit)
};

namespace Tmpl8 {

__declspec(align(64)) class Ray
{
public:
	Ray() = default;
	Ray(const Ray&) {}
	Ray(float3 origin, float3 direction, float distance = 1e34f, int depth = 0 )
	{
		O = origin, D = direction, I.t = distance;
		// calculate reciprocal ray direction for triangles and AABBs
		rD = float3( 1 / D.x, 1 / D.y, 1 / D.z );
		depthidx = depth;
	#ifdef SPEEDTRIX
		d0 = d1 = d2 = 0;
	#endif
	}
	float3 IntersectionPoint() { return O + I.t * D; }

	Ray Reflect(float3 I,float3 N) { 
		// create a secondary ray
		// intersection = origin, norm = norm from intersected object
		float3 dir = this->D;
		float3 reflected = normalize(dir - 2.0f * (dot(dir, N) * N));
		return Ray(I + (0.0002f * reflected), reflected, 1e34f, depthidx + 1);
		
	}
	// ray data
#ifndef SPEEDTRIX
	float3 O, D, rD;
#else
	union { struct { float3 O; float d0; }; __m128 O4; };
	union { struct { float3 D; float d1; }; __m128 D4; };
	union { struct { float3 rD; float d2; }; __m128 rD4; };
#endif
	Intersection I;
	int depthidx = 0;
	//float3 nearestcolor;
	//float4 nearestmat;
	//float3 Norm_surf;
};

// -----------------------------------------------------------
// Sphere primitive
// Basic sphere, with explicit support for rays that start
// inside it. Good candidate for a dielectric material.
// -----------------------------------------------------------
class Sphere
{
public:
	Sphere() = default;
	Sphere(uint idx, float3 p, float r ) : 
		pos(p), r2(r* r), invr(1 / r), objIdx(idx) {
	}
	Sphere(uint idx, float3 p, float r, Material mat) : Sphere(idx, p, r) {
		material = mat;
	}
	void Intersect( Ray& ray ) const
	{
		float3 oc = ray.O - this->pos;
		float b = dot( oc, ray.D );
		float c = dot( oc, oc ) - this->r2;
		float t, d = b * b - c;
		if (d <= 0) return;
		d = sqrtf( d ), t = -b - d;
		if (t < ray.I.t && t > 0)
		{
			ray.I.t = t, ray.I.instPrim = objIdx << 20;
			return;
		}
		t = d - b;
		if (t < ray.I.t && t > 0)
		{
			ray.I.t = t, ray.I.instPrim = objIdx << 20;
			return;
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		return (I - this->pos) * invr;
	}
	float3 GetAlbedo( const float3 I ) const
	{
		return float3( 0.93f );
	}
	float3 pos = 0;
	float r2 = 0, invr = 0;
	uint objIdx = 0;

	float3 color = float3(0, 0, 1);
	static int id; 
	Material material = {
		float3(0, 0, 1), //albedo
		0.0, //specularity
		Medium::Undefined, //medium
	};
};

// -----------------------------------------------------------
// Plane primitive
// Basic infinite plane, defined by a normal and a distance
// from the origin (in the direction of the normal).
// -----------------------------------------------------------
class Plane
{
public:
	Plane() = default;
	Plane(uint idx, float3 normal, float dist, float3 col = (1, 0.5, 0.5)) : N(normal), d(dist), objIdx(idx){ color = col; }
	void Intersect( Ray& ray ) const
	{
		float t = -(dot( ray.O, this->N ) + this->d) / (dot( ray.D, this->N ));
		if (t < ray.I.t && t > 0) ray.I.t = t, ray.I.instPrim = objIdx << 20;
	}
	float3 GetNormal( const float3 I ) const
	{
		return N;
	}
	float3 GetAlbedo( const float3 I ) const
	{
		if (N.y == 1)
		{
			// floor albedo: checkerboard
			int ix = (int)(I.x * 2 + 96.01f);
			int iz = (int)(I.z * 2 + 96.01f);
			// add deliberate aliasing to two tile
			if (ix == 98 && iz == 98) ix = (int)(I.x * 32.01f), iz = (int)(I.z * 32.01f);
			if (ix == 94 && iz == 98) ix = (int)(I.x * 64.01f), iz = (int)(I.z * 64.01f);
			return float3( ((ix + iz) & 1) ? 1 : 0.3f );
		}
		else if (N.z == -1)
		{
			// back wall: logo
			static Surface logo( "assets/logo.png" );
			int ix = (int)((I.x + 4) * (128.0f / 8));
			int iy = (int)((2 - I.y) * (64.0f / 3));
			uint p = logo.pixels[(ix & 127) + (iy & 63) * 128];
			uint3 i3( (p >> 16) & 255, (p >> 8) & 255, p & 255 );
			return float3( i3 ) * (1.0f / 255.0f);
		}
		return float3( 0.93f );
	}
	float3 N;
	float d;
	uint objIdx = 0;
	float3 color;
};

// -----------------------------------------------------------
// Cube primitive
// Oriented cube. Unsure if this will also work for rays that
// start inside it; maybe not the best candidate for testing
// dielectrics.
// -----------------------------------------------------------
class Cube
{
public:
	Cube() = default;
	Cube( uint idx, float3 pos, float3 size, mat4 transform = mat4::Identity() )
	{
		objIdx = idx;
		b[0] = pos - 0.5f * size, b[1] = pos + 0.5f * size; //-
		M = transform, invM = transform.FastInvertedTransformNoScale();
		
	}
	void Intersect( Ray& ray ) const
	{
		// 'rotate' the cube by transforming the ray into object space
		// using the inverse of the cube transform.
		float3 O = TransformPosition( ray.O, invM );
		float3 D = TransformVector( ray.D, invM );
		float rDx = 1 / D.x, rDy = 1 / D.y, rDz = 1 / D.z;
		int signx = D.x < 0, signy = D.y < 0, signz = D.z < 0;
		float tmin = (b[signx].x - O.x) * rDx;
		float tmax = (b[1 - signx].x - O.x) * rDx;
		float tymin = (b[signy].y - O.y) * rDy;
		float tymax = (b[1 - signy].y - O.y) * rDy;
		if (tmin > tymax || tymin > tmax) return;
		tmin = max( tmin, tymin ), tmax = min( tmax, tymax );
		float tzmin = (b[signz].z - O.z) * rDz;
		float tzmax = (b[1 - signz].z - O.z) * rDz;
		if (tmin > tzmax || tzmin > tmax) return;
		tmin = max( tmin, tzmin ), tmax = min( tmax, tzmax );
		if (tmin > 0)
		{
			if (tmin < ray.I.t) ray.I.t = tmin, ray.I.instPrim = objIdx << 20;
		}
		else if (tmax > 0)
		{
			if (tmax < ray.I.t) ray.I.t = tmax, ray.I.instPrim = objIdx << 20;
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		// transform intersection point to object space
		float3 objI = TransformPosition( I, invM );
		// determine normal in object space
		float3 N = float3( -1, 0, 0 );
		float d0 = fabs( objI.x - b[0].x ), d1 = fabs( objI.x - b[1].x );
		float d2 = fabs( objI.y - b[0].y ), d3 = fabs( objI.y - b[1].y );
		float d4 = fabs( objI.z - b[0].z ), d5 = fabs( objI.z - b[1].z );
		float minDist = d0;
		if (d1 < minDist) minDist = d1, N.x = 1;
		if (d2 < minDist) minDist = d2, N = float3( 0, -1, 0 );
		if (d3 < minDist) minDist = d3, N = float3( 0, 1, 0 );
		if (d4 < minDist) minDist = d4, N = float3( 0, 0, -1 );
		if (d5 < minDist) minDist = d5, N = float3( 0, 0, 1 );
		// return normal in world space
		return TransformVector( N, M );
	}
	float3 b[2];
	mat4 M, invM;
	uint objIdx = 0;
	Material material = {
		float3(0, 1, 0), //albedo
		0.2, //specularity
		Medium::Undefined, //medium
	};
};

// -----------------------------------------------------------
// Triangle primitive
// Just a triangle (with normal and centroid?).
// Should probably remove this when we start using meshes,
// with vertices and triangles just being defined as a list of indices
// -----------------------------------------------------------
class Triangle {
public:
	Triangle() = default;
	Triangle(float3 v0, float3 v1, float3 v2, uint id) : vertex0(v0), vertex1(v1), vertex2(v2) {
		centroid = (v0 + v1 + v2) / 3.0;
		normal = normalize(cross((v1 - v0), (v2 - v0)));
		objIdx = 10 + id;
	}
	Triangle(float3 v0, float3 v1, float3 v2, int id, Material mat) : Triangle(v0, v1, v2, id) {
		material = mat;
	}
	// Ripped straight from your BVH Tutorial, Jacco <3
	void Intersect(Ray& ray) const
	{
		const float3 edge1 = vertex1 - vertex0;
		const float3 edge2 = vertex2 - vertex0;
		const float3 h = cross(ray.D, edge2);
		const float a = dot(edge1, h);
		if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
		const float f = 1 / a;
		const float3 s = ray.O - vertex0;
		const float u = f * dot(s, h);
		if (u < 0 || u > 1) return;
		const float3 q = cross(s, edge1);
		const float v = f * dot(ray.D, q);
		if (v < 0 || u + v > 1) return;
		const float t = f * dot(edge2, q);
		if (t > 0.0001f && t < ray.I.t) {
			ray.I.t = t;	
			ray.I.instPrim = objIdx << 20;
		}
	}
	float3 GetNormal() const 
	{
		return normal;
	}
	union {
		struct { float3 vertex0, vertex1, vertex2; };
		float3 cell[3];
	};
	uint objIdx;
	float3& operator [] (const bool n) { return cell[n]; }
	float3 centroid;
	float3 normal;
	Material material = {
		float3(1, 0, 0), //albedo
		0.0, //specularity
		Medium::Undefined, //medium
	};
};

// -----------------------------------------------------------
// Smaller Triangle struct for use in Mesh
// Perhaps this is stupid.
// -----------------------------------------------------------
__declspec(align(64)) struct Tri {
	float3 vertex0, vertex1, vertex2;
	float3 centroid;
};


// -----------------------------------------------------------
// Holds Tri data for texturing and shading
// This is 
// -----------------------------------------------------------
__declspec(align(64)) struct TriEx {
	float2 uv0, uv1, uv2;
	float3 N0, N1, N2;
};

// -----------------------------------------------------------
// Mesh class, partially yoinked from Jacco BVH tutorial
//
// -----------------------------------------------------------
class RTXMesh
{
public:
	RTXMesh() = default;
	RTXMesh(uint objId, mat4 transform = mat4::Identity()) {
		objIdx = objId;
		M = transform;
		invM = transform.FastInvertedTransformNoScale();
	};
	void Intersect(Ray& ray) {
		for (uint i = 0; i < triCount; i++) {
			const float3 edge1 = tri[i].vertex1 - tri[i].vertex0;
			const float3 edge2 = tri[i].vertex2 - tri[i].vertex0;
			const float3 h = cross(ray.D, edge2);
			const float a = dot(edge1, h);
			if (a > -0.0001f && a < 0.0001f) continue; // ray parallel to triangle
			const float f = 1 / a;
			const float3 s = ray.O - tri[i].vertex0;
			const float u = f * dot(s, h);
			if (u < 0 || u > 1) continue;
			const float3 q = cross(s, edge1);
			const float v = f * dot(ray.D, q);
			if (v < 0 || u + v > 1) continue;
			const float t = f * dot(edge2, q);
			if (t > 0.0001f && t < ray.I.t) {
				ray.I.t = t;
				ray.I.u = u;
				ray.I.v = v;
				ray.I.instPrim = (objIdx << 20) + i;
			}
		}
	}
	float3 GetColor(Intersection& I) {
		const uint mask = ~(~0 << 20); // Mask for last 20 bits
		uint id = I.instPrim & mask;
		if (textureLoaded) {
			float2 uv = I.u * triEx[id].uv1 + I.v * triEx[id].uv2 + (1.0f - (I.u + I.v)) * triEx[id].uv0;
			int iu = (int)(uv.x * texture->width) % texture->width;
			int iv = (int)(uv.y * texture->height) % texture->height;
			uint texel = texture->pixels[iu + iv * texture->width];
			bitset<32> texbit(texel);
			//cout << "Texture pixel found: " << texbit << endl;
			unsigned char x = (texel >> 16);
			unsigned char y = (texel >> 8);
			unsigned char z = texel;
			float3 returnval(x / 255.f, y / 255.f, z / 255.f);
			//cout << "Returning:\t" << returnval.x << "|" << returnval.y << "|" << returnval.z << endl;
			return returnval;
		}
		else {
			return material.albedo;
		}
	}
	float3 GetNormal(Intersection& I) {
		const uint mask = ~(~0 << 20); // Mask for last 20 bits
		uint id = I.instPrim & mask;
		if (normalMapLoaded) {
			float2 uv = I.u * triEx[id].uv1 + I.v * triEx[id].uv2 + (1.0f - (I.u + I.v)) * triEx[id].uv0;
			int iu = (int)(uv.x * normalMap->width) % normalMap->width;
			int iv = (int)(uv.y * normalMap->height) % normalMap->height;
			uint texel = normalMap->pixels[iu + iv * normalMap->width];
			unsigned char x = (texel >> 16);
			unsigned char y = (texel >> 8);
			unsigned char z = texel;
			return float3(x / 255.f, y / 255.f, z / 255.f);
		}
		else {
			float3 N = I.u * triEx[id].N1 + I.v * triEx[id].N2 + (1.0f - (I.u + I.v)) * triEx[id].N0;
			return normalize(N);
		}
	}
	static RTXMesh* makeMesh(aiMesh* mesh, const aiScene* scene, string const& path, uint objId) {
		RTXMesh* m = new RTXMesh(objId);
		m->tri = new Tri[mesh->mNumFaces];
		m->triEx = new TriEx[mesh->mNumFaces];
		string directory = path.substr(0, path.find_last_of('/'));

		// We load things face by face instead of vertex by vertex. ray.I.tracers like accessing triangles compactly, I think.
		// Potentially can be changed into a list of vertices and a list of face indices, to save memory.
		for (uint i = 0; i < mesh->mNumFaces; i++) {
			aiFace face = mesh->mFaces[i];
			// vertex0
			m->tri[i].vertex0.x = mesh->mVertices[face.mIndices[0]].x;
			m->tri[i].vertex0.y = mesh->mVertices[face.mIndices[0]].y;
			m->tri[i].vertex0.z = mesh->mVertices[face.mIndices[0]].z;
			// vertex1
			m->tri[i].vertex1.x = mesh->mVertices[face.mIndices[1]].x;
			m->tri[i].vertex1.y = mesh->mVertices[face.mIndices[1]].y;
			m->tri[i].vertex1.z = mesh->mVertices[face.mIndices[1]].z;
			// vertex2
			m->tri[i].vertex2.x = mesh->mVertices[face.mIndices[2]].x;
			m->tri[i].vertex2.y = mesh->mVertices[face.mIndices[2]].y;
			m->tri[i].vertex2.z = mesh->mVertices[face.mIndices[2]].z;
			// Calculate centroid
			m->tri[i].centroid = (m->tri[i].vertex0 + m->tri[i].vertex1 + m->tri[i].vertex2) / 3.0f;

			// Load the normals
			if (mesh->HasNormals()) {
				// Vertex 0 Normal
				m->triEx[i].N0.x = mesh->mNormals[face.mIndices[0]].x;
				m->triEx[i].N0.y = mesh->mNormals[face.mIndices[0]].y;
				m->triEx[i].N0.z = mesh->mNormals[face.mIndices[0]].z;
				// Vertex 1 Normal
				m->triEx[i].N1.x = mesh->mNormals[face.mIndices[1]].x;
				m->triEx[i].N1.y = mesh->mNormals[face.mIndices[1]].y;
				m->triEx[i].N1.z = mesh->mNormals[face.mIndices[1]].z;
				// Vertex 2 Normal
				m->triEx[i].N2.x = mesh->mNormals[face.mIndices[2]].x;
				m->triEx[i].N2.y = mesh->mNormals[face.mIndices[2]].y;
				m->triEx[i].N2.z = mesh->mNormals[face.mIndices[2]].z;
				// Normalize them, because for some reason they are not.
				m->triEx[i].N0 = normalize(m->triEx[i].N0);
				m->triEx[i].N1 = normalize(m->triEx[i].N1);
				m->triEx[i].N2 = normalize(m->triEx[i].N2);
			}

			// Load the texture coordinates
			if (mesh->mTextureCoords[0]) {
				// Texture coordinates for vertex0. I think.
				m->triEx[i].uv0.x = mesh->mTextureCoords[0][face.mIndices[0]].x;
				m->triEx[i].uv0.y = mesh->mTextureCoords[0][face.mIndices[0]].y;
				// Texture coordinates for vertex1. I hope.
				m->triEx[i].uv1.x = mesh->mTextureCoords[0][face.mIndices[1]].x;
				m->triEx[i].uv1.y = mesh->mTextureCoords[0][face.mIndices[1]].y;
				// Texture coordinates for vertex2. I pray.
				m->triEx[i].uv2.x = mesh->mTextureCoords[0][face.mIndices[2]].x;
				m->triEx[i].uv2.y = mesh->mTextureCoords[0][face.mIndices[2]].y;
			}
			else {
				m->material.albedo = float3(1, 0, 0); //default red color
			}
		}
		m->triCount = mesh->mNumFaces;
		// Time to get the textures (and normal map)
		aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
		// Diffuse (aka just texture?)
		if (material->GetTextureCount(aiTextureType_DIFFUSE) > 0) {
			aiString str;
			material->GetTexture(aiTextureType_DIFFUSE, 0, &str);
			string filename = string(str.C_Str());
			filename = directory + "/" + filename;
			m->texture = new Surface(filename.c_str());
			m->textureLoaded = true;
		}
		// Bump map!
		if (material->GetTextureCount(aiTextureType_HEIGHT) > 0) {
			cout << "Loaded normal map??" << endl;
			aiString str;
			material->GetTexture(aiTextureType_HEIGHT, 0, &str);
			string filename = string(str.C_Str());
			filename = directory + "/" + filename;
			m->normalMap = new Surface(filename.c_str());
			m->normalMapLoaded = true;
		}
		return m;
	}
	Tri* tri;
	TriEx* triEx;
	int triCount;
	uint objIdx;
	mat4 M, invM;
	Material material;
	Surface* texture;
	bool textureLoaded = false;
	Surface* normalMap;
	bool normalMapLoaded = false;
};

// -----------------------------------------------------------
// Quad primitive
// Oriented quad, intended to be used as a light source.
// -----------------------------------------------------------
class Quad
{
public:
	Quad() = default;
	Quad( uint idx, float s, mat4 transform = mat4::Identity() )
	{
		objIdx = idx;
		size =  s * 0.5f;
		T = transform, invT = transform.FastInvertedTransformNoScale();
		
	}
	void Intersect( Ray& ray ) const
	{
		const float3 O = TransformPosition( ray.O, invT );
		const float3 D = TransformVector( ray.D, invT );
		const float t = O.y / -D.y;
		if (t < ray.I.t && t > 0)
		{
			float3 I = O + t * D;
			if (I.x > -size && I.x < size && I.z > -size && I.z < size)
				ray.I.t = t, ray.I.instPrim = objIdx << 20;
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		// TransformVector( float3( 0, -1, 0 ), T ) 
		return float3( -T.cell[1], -T.cell[5], -T.cell[9] );
	}
	float3 GetAlbedo( const float3 I ) const
	{
		return float3( 10 );
	}
	Material material = {
		float3(1, 1, 0), //albedo
		0.2, //specularity
		Medium::Undefined, //medium
	};
	float size;
	mat4 T, invT;
	uint objIdx = 0;
};

// -----------------------------------------------------------
// Scene class
// We intersect this. The query is internally forwarded to the
// list of primitives, so that the nearest hit can be returned.
// For this hit (distance, obj id), we can query the normal and
// albedo.
// -----------------------------------------------------------
class Scene
{
public:
	Scene()
	{
		// Load cat
		//loadModel("assets/wolf/Wolf.obj");
		loadModel("assets/chessboard/chessboard.obj");
		cout << meshPool.size() << endl;
		if(meshPool[0]->textureLoaded) cout << meshPool[0]->texture->width << "||" << meshPool[0]->texture->height << endl;
		for (int i = 0; i < meshPool[0]->triCount; i++) {
			float3 N0 = meshPool[0]->triEx[i].N0;
			float3 N1 = meshPool[0]->triEx[i].N1;
			float3 N2 = meshPool[0]->triEx[i].N2;
			cout << i << ":N0\t" << N0.x << "||" << N0.y << "||" << N0.z << endl;
			cout << i << ":N1\t" << N1.x << "||" << N1.y << "||" << N1.z << endl;
			cout << i << ":N2\t" << N2.x << "||" << N2.y << "||" << N2.z << endl;
		}
		/*for (int i = 0; i < meshPool[0]->triCount; i++) {
			Intersection I;
			I.instPrim = i;
			I.u = 1.f / 3.f;
			I.v = 1.f / 3.f;
			float3 N = meshPool[0]->GetNormal(I);
			cout << i << ":\t" << N.x << "|" << N.y << "|" << N.z << endl;
			cout << "Length:\t" << length(N) << endl;
		}*/
		cout << "Texture:" << endl;
		bitset<32> black(meshPool[0]->texture->pixels[0]);
		bitset<32> white(meshPool[0]->texture->pixels[255]);
		cout << "Black: " << black << endl;
		cout << "White: " << white << endl;
		cout << "Triangles: " << meshPool[0]->triCount << endl;


		// Skybox
		unsigned char* data = stbi_load("assets/clarens_midday_4k.png", &width, &height, &nrChannels, 0);
		if (data) {
			skybox = (unsigned char*)MALLOC64(width * height * 4 * sizeof(char));
			loadedSkybox = true;
			memcpy(skybox, data, width * height * 4);
			stbi_image_free(data);
		}
		else {
			cout << stbi_failure_reason() << endl;
			throw exception("Failed to load Skybox.");
		}

		lights[0] = Light(float3(0, 80, 0.5), 240);
		// Precalc refractive mappings
		refractiveIndex[Medium::Air] = 1.0;
		refractiveIndex[Medium::Glass] = 1.52;
		for (int i = 0; i < refractiveIndex.size(); i++) {
			Medium a = static_cast<Medium>(i);
			for (int j = 0; j < refractiveIndex.size(); j++) {
				Medium b = static_cast<Medium>(j);
				if (j == i) refractiveTransmissions[a][b] = 1.0;
				else refractiveTransmissions[a][b] = refractiveIndex[a] / refractiveIndex[b];
			}
		}
		// Precalc some materials
		Material mirror = Material(float3(0, 1.0, 0), 1.0, Medium::Undefined);
		materialMap["Mirror"] = mirror;


		// we store all primitives in one continuous buffer
		//quad = Quad( 0, 1 );									// 1: light source
		lightSphere = Sphere(1, lights[0].position, 5);
		lightSphere.material.albedo = (1, 1, 1);
		sphere = Sphere( 2, float3( 0 ), 10.0f);					// 2: bouncing ball   0.5f
		sphere2 = Sphere( 3, float3( 0, -1.5f, -1.05f ), 0.5f );	// 3: glass ball
		sphere2.material.albedo = float3(0.2, 0.8, 0.2);
		sphere2.material.mat_medium = Medium::Glass;
		sphere2.material.absorption = float3(0.2f, 2.0f, 4.0f);
		cube = Cube(4, float3(0), float3(1.15f));				// 4: cube 		cube = Cube( 3, float3( 0 ), float3( 1.15f ) );		
		cube2 = Cube(5, float3(0, -1.0f, -2.0f), float3(0.8f)); // 5 Glass cube
		cube2.material.mat_medium = Medium::Glass;
		cube2.material.albedo = float3(0.2f, 0.2f, 0.2f);
		cube2.material.absorption = float3(0);
		
		Triangle* triangle = new Triangle(float3(-0.9f, 0, -1), float3(0.2f, 0, -1), float3(0, 1.0f, -0.5), 0); // 10: triangle
		// Left wall
		Triangle* wall_l0 = new Triangle(float3(-3, -3, 3), float3(-3, -3, -3), float3(-3, 3, 3), 1);
		Triangle* wall_l1 = new Triangle(float3(-3, 3, -3), float3(-3, 3, 3), float3(-3, -3, -3), 2);
		// Right wall mirro
		Triangle* wall_r0 = new Triangle(float3(2.99, -2.5, 2.5), float3(2.99, 2.5, 2.5), float3(2.99, -2.5, -2.5), 3);
		wall_r0->material = materialMap["Mirror"];
		Triangle* wall_r1 = new Triangle(float3(2.99, 2.5, -2.5), float3(2.99, -2.5, -2.5), float3(2.99, 2.5, 2.5), 4);
		wall_r1->material = materialMap["Mirror"];
		// Right wall backdrop
		Triangle* wall_r_backdrop0 = new Triangle(float3(3.0, -3.0, 3), float3(3, 3, 3), float3(3, -3, -3), 5);
		Triangle* wall_r_backdrop1 = new Triangle(float3(3.0, 3.0, -3), float3(3, -3, -3), float3(3, 3, 3), 6);
		wall_r_backdrop0->material.albedo = float3(0.5, 0.5, 0);
		wall_r_backdrop1->material.albedo = float3(0.5, 0.5, 0);
		// Ceiling
		Triangle* ceiling_0 = new Triangle(float3(3, 3, 3), float3(-3, 3, 3), float3(3, 3, -3), 7);
		Triangle* ceiling_1 = new Triangle(float3(-3, 3, -3), float3(3, 3, -3), float3(-3, 3, 3), 8);
		ceiling_0->material.albedo = float3(0.0, 0.5, 0.5);
		ceiling_1->material.albedo = float3(0.0, 0.5, 0.5);
		// Back wall
		Triangle* wall_b0 = new Triangle(float3(-3, -3, 3), float3(-3, 3, 3), float3(3, -3, 3), 9);
		Triangle* wall_b1 = new Triangle(float3(3, 3, 3), float3(3, -3, 3), float3(-3, 3, 3), 10);
		wall_b0->material.albedo = float3(0.5, 0, 0.5);
		wall_b1->material.albedo = float3(0.5, 0, 0.5);
		// floor
		Triangle* floor_0= new Triangle(float3(3, -3, 3), float3(3, -3, -3), float3(-3, -3, 3), 11);
		Triangle* floor_1 = new Triangle(float3(-3, -3, -3), float3(-3, -3, 3), float3(3, -3, -3), 12);
		floor_0->material.albedo = float3(0.2, 1, 0);
		floor_1->material.albedo = float3(0.2, 1, 0);
		trianglePool.push_back(triangle);
		trianglePool.push_back(wall_l0);
		trianglePool.push_back(wall_l1);
		trianglePool.push_back(wall_r0);
		trianglePool.push_back(wall_r1);
		trianglePool.push_back(wall_r_backdrop0);
		trianglePool.push_back(wall_r_backdrop1);
		trianglePool.push_back(ceiling_0);
		trianglePool.push_back(ceiling_1);
		trianglePool.push_back(wall_b0);
		trianglePool.push_back(wall_b1);
		trianglePool.push_back(floor_0);
		trianglePool.push_back(floor_1);
		SetTime( 0 );
		// Note: once we have triangle support we should get rid of the class
		// hierarchy: virtuals reduce performance somewhat.
	}

	Scene::~Scene() {
		if (loadedSkybox) FREE64(skybox);
	}

	void SetTime( float t )
	{
		// default time for the scene is simply 0. Updating/ the time per frame 
		// enables animation. Updating it per ray can be used for motion blur.
		animTime = t;
		// light source animation: swing
		mat4 M1base = mat4::Translate( float3( 0, 2.6f, 2 ) );
		mat4 M1 = M1base * mat4::RotateZ( sinf( animTime * 0.6f ) * 0.1f ) * mat4::Translate( float3( 0, -0.9, 0 ) );
		//quad.T = M1, quad.invT = M1.FastInvertedTransformNoScale();
		// cube animation: spin
		mat4 M2base = mat4::RotateX( PI / 4 ) * mat4::RotateZ( PI / 4 );
		mat4 M2 = mat4::Translate( float3( 1.4f, 0, 2 ) ) * mat4::RotateY( animTime * 0.5f ) * M2base;
		cube.M = M2, cube.invM = M2.FastInvertedTransformNoScale();
		// sphere animation: bounce
		float tm = 1 - sqrf( fmodf( animTime, 2.0f ) - 1 );
		sphere.pos = float3( -1.4f, -0.5f + tm, 2 );
	}

	float3 getSkyBox(float3 dir) const {
		float phi = atan2f(dir.x, dir.z);
		float theta = atan2f(hypotf(dir.x, dir.z), dir.y);
		if (phi < 0) phi += 2.0f * PI;
		if (theta < 0) throw exception("Negative theta?");
		float x = phi / (2.0f * PI);
		float y = theta / PI;
		int x_ = (int)(x * width);
		int y_ = (int)(y * height);
		unsigned char* pixel = skybox + (y_ * width + x_) * nrChannels;
		unsigned char r = pixel[0];
		unsigned char g = pixel[1];
		unsigned char b = pixel[2];
		return float3(r / 255.f, g / 255.f, b / 255.f);
	}

	float3 GetLightPos() const
	{
		// light point position is the middle of the swinging quad
		//float3 corner1 = TransformPosition( float3( -0.5f, 0, -0.5f ), quad.T  );
		//float3 corner2 = TransformPosition( float3( 0.5f, 0, 0.5f ) , quad.T );
		//return (corner1 + corner2) * 0.5f - float3( 0, 0.01f, 1 );
		return lights[0].position;
	}

	Material getMaterial(uint objIdx) const
	{
		uint objId = objIdx >> 20;
		if (objId == 0) throw exception("There's no material for nothing"); // or perhaps we should just crash
		//if (objIdx == 0) return quad.material;
		if (objId == 1) return lightSphere.material; //Jax:  chancged idx to 2
		else if (objId == 2) return sphere.material;
		else if (objId == 3) return sphere2.material;
		else if (objId == 4) return cube.material;
		else if (objId == 5) return cube2.material;
		else if (objId >= 10) return trianglePool[objId - 10]->material;
		throw exception("ID not known");
	}

	float3 GetLightColor() const
	{
		return lights[0].color; 
	}

	float3 GetDiffuseRefelectDir(float3 N) {
	// steps: create a random cube around the hit point
	// find a random point in the cube 
	// check if greater than one -> if so, recalculate. 
	// if not, normlize. 
	// check if dot(Point, Normal) <0 -> p = -p
	// create a cobe -1 1 around the point
		bool gen = false;
		float3 random_vec;
		while (!gen) {
			float x = (2*RandomFloat())-1;
			float y = (2 * RandomFloat()) - 1;
			float z = (2 * RandomFloat()) - 1;

			random_vec = float3( x,  y,  z);
			if (length(random_vec) < 1) {
				random_vec = normalize(random_vec); // puting on the unit circle
				if (dot(N, random_vec) < 0) {
					random_vec = -random_vec;
				}
				gen = true;
			}
		}
		return random_vec;
	}

	float3 directIllumination(uint objIdx, float3 intersection, float3 norm, float3 albedo) {
		uint objId = objIdx >> 20;
		//cout << "Direct Illumination " << endl;
		if (objId == 1) return lightSphere.material.albedo;
		float3 color = (0, 0, 0);
		for (int i = 0; i < numLightSouces; i++) { 
			float3 dir_light = normalize(lights[i].position - intersection);
			float dot_p = dot(dir_light, norm);
			//cout << "dotP: " << dot_p << endl;
			float intensity = lights[i].lumen / pow(length(lights[i].position - intersection), 2);

			// -----------------------------------------------------------
			// regular
			// -----------------------------------------------------------
			color += intensity * albedo * maxFloat( dot_p, 0.0f); //
			//cout << "Color in directIlum: " << color.x << "|" << color.y << "|" << color.z << endl;
		
			// -----------------------------------------------------------
			// For specular highlights
			// color += normalize(GetLightColor())* maxFloat(dot_p* dot_p, 0.0f); //
			// -----------------------------------------------------------
		}
		
		return color;
	}
	void FindNearest( Ray& ray ) const
	{
		// room walls - ugly shortcut for more speed
		float t;
		lightSphere.Intersect(ray);
		//sphere.Intersect( ray );
		//sphere2.Intersect( ray );
		//cube.Intersect( ray );
		//cube2.Intersect(ray);
		//for (int i = 0; i < (int)trianglePool.size(); i++) {
		//	trianglePool[i]->Intersect(ray);
		//}
		for (RTXMesh* mesh : meshPool) {
			mesh->Intersect(ray);
		}
		//if (ray.objIdx >= 10) cout << "Triangle hit: " << ray.objIdx << endl;
	}
	bool IsOccluded( Ray& ray ) const
	{
		float rayLength = ray.I.t;
		// skip planes: it is not possible for the walls to occlude anything
		//quad.Intersect( ray );
		//sphere.Intersect( ray );
		//sphere2.Intersect( ray );
		//cube.Intersect( ray );
		//cube2.Intersect(ray);
		//for (int i = 0; i < (int)trianglePool.size(); i++) {
		//	trianglePool[i]->Intersect(ray);
		//}
		for (RTXMesh* mesh : meshPool) mesh->Intersect(ray);
		return ray.I.t < rayLength;
		// technically this is wasteful: 
		// - we potentially search beyond rayLength
		// - we store objIdx and t when we just need a yes/no
		// - we don't 'early out' after the first occlusion
	}
	float3 GetNormal(Intersection& Inters, float3 I, float3 wo, bool& hit_back) const
	{
		// we get the normal after finding the nearest intersection:
		// this way we prevent calculating it multiple times.
		uint objId = Inters.instPrim >> 20;
		if (objId == 0) return float3( 0 ); // or perhaps we should just crash
		float3 N;
		//if (objIdx == 0) N = quad.GetNormal(I);
		if (objId == 1) N = lightSphere.GetNormal(I); // check objIDX to 2 from 0
		//if (objId == 2) N = sphere.GetNormal(I);
		//else if (objId == 3) N = sphere2.GetNormal(I);
		//else if (objId == 4) N = cube.GetNormal(I);
		//else if (objId == 5) N = cube2.GetNormal(I);
		//else if (objId >= 10) N = trianglePool[objId - 10]->GetNormal();
		else if (objId > 1) meshPool[objId - 2]->GetNormal(Inters);
		//else 
		//{
		//	// faster to handle the 6 planes without a call to GetNormal
		//	N = float3( 0 );
		//	N[(objIdx - 4) / 2] = 1 - 2 * (float)(objIdx & 1);
		//}
		if (dot(N, wo) > 0) { N = -N; hit_back = true; } // hit backside / inside
		return N;
	}

	float3 GetColor(Intersection& I) {
		uint objId = I.instPrim >> 20;
		if (objId == 1) return lightSphere.material.albedo;
		else if (objId > 1) 
		{
			return meshPool[objId - 2]->GetColor(I); 
		}
	}

	// Mesh loading using tutorial from learnopengl.com
	void loadModel(const char* file) {
		Assimp::Importer importer;
		const aiScene* scene = importer.ReadFile(file, aiProcess_Triangulate | aiProcess_FlipUVs);
		if (!scene || scene->mFlags * AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
			cout << "Error loading Mesh: " << importer.GetErrorString() << endl;
			return;
		}

		processNode(scene->mRootNode, scene, file);
	}

	void processNode(aiNode* node, const aiScene* scene, const char* path) {
		// process parent node
		for (int i = 0; i < node->mNumMeshes; i++) {
			int id = meshPool.size() + 2; //TODO: REMOVE THE START AT 2
			aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
			meshPool.push_back(RTXMesh::makeMesh(mesh, scene, path, id));
		}
		// go through children nodes
		for (uint i = 0; i < node->mNumChildren; i++) {
			processNode(node->mChildren[i], scene, path);
		}
	}

	__declspec(align(64)) // start a new cacheline here
	float animTime = 0;
	const float ambient = 0.005;
	//Quad quad;
	Light lights[1];
	Sphere sphere;
	Sphere lightSphere;
	Sphere sphere2;
	Cube cube;
	Cube cube2;
	//Plane plane[6];
	static inline vector<Triangle*> trianglePool;
	static inline vector<Sphere*> spherePool;
	static inline vector<Cube*> cubePool;
	static inline vector<RTXMesh*> meshPool;
	static inline map<string, Material> materialMap;
	unsigned char* skybox;
	int width, height, nrChannels;
	bool loadedSkybox = false;
	map<Medium, float> refractiveIndex;
	map<Medium, map<Medium, float>> refractiveTransmissions;
	int numLightSouces=1; // should be initilize better in the future
};

}