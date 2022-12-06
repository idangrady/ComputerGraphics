#pragma once
#include <map>
#include <random>
#include <../lib/stb_image.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <config.h>
// -----------------------------------------------------------
// scene.h
// Simple test scene for ray tracing experiments. Goals:
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

namespace Tmpl8 {
	const uint sphereID = 0;
	const uint cubeID = 1;
	const uint triangleID = 2;
	const uint meshID = 3;
	const uint planeID = 4;
	const uint lightID = 5;
	const uint skyBoxID = 6;

	static inline uint MakeID(uint type, uint id, uint tri) {
		return (type << 29) + (id << 20) + (tri);
	}

	static inline uint GetObjectType(uint idx) {
		const uint mask = ~0 << 29;
		return (idx & mask) >> 29;
	}
	static inline uint GetObjectIndex(uint idx) {
		const uint mask = ((~0 << 20) & ~(7 << 29)); // This mask is giving me a headache;
		return (idx & mask) >> 20;
	}
	static inline uint GetTriangleIndex(uint idx) {
		const uint mask = ~(~0 << 20);
		return idx & mask;
	}
	enum class Medium {
		Undefined = -1,
		Air = 0,
		Glass = 1,
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
		uint instPrim = MakeID(skyBoxID, 0, 0);	// Type indedx (3 bit), instance index (9 bit) and primitive index (20 bit)
	};

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
			ray.I.t = t, ray.I.instPrim = MakeID(sphereID, objIdx, 0);
			return;
		}
		t = d - b;
		if (t < ray.I.t && t > 0)
		{
			ray.I.t = t, ray.I.instPrim = MakeID(sphereID, objIdx, 0);
			return;
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		return (I - this->pos) * invr;
	}
	float3 GetAlbedo() const
	{
		return material.albedo;
	}
	float3 pos = 0;
	float r2 = 0, invr = 0;
	uint objIdx = 0;
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
		if (t < ray.I.t && t > 0) ray.I.t = t, ray.I.instPrim = MakeID(planeID, objIdx, 0);
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
			if (tmin < ray.I.t) ray.I.t = tmin, ray.I.instPrim = Tmpl8::MakeID(cubeID, objIdx, 0);;
		}
		else if (tmax > 0)
		{
			if (tmax < ray.I.t) ray.I.t = tmax, ray.I.instPrim = Tmpl8::MakeID(cubeID, objIdx, 0);;
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
		objIdx = id;
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
			ray.I.instPrim = Tmpl8::MakeID(triangleID, objIdx, 0);;
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


class light {

public:
	light() = default;
	light(float i, uint8_t id) : lumen(i), idx(id) {}
	float3 color{ 1.0f, 1.0f, 1.0f };
	float lumen;
	uint idx;
	Material material;

	void setLightColor(float3 c) { color = c; };
	uint8_t getidx() { return idx; }
	Material getMaterial() const { return material; }
	float3 getLightColor() const { return lumen* color; };
};
//class Triangle; class Ray; class Sphere;
class areaLight : public light
{
public:
	areaLight() = default;
	areaLight(float i, uint id, float3 p1_, float3 p2_, float3 p3_, float3 p4_, float3 p5_, float3 p6_) :light(i, id) {
		p1 = p1_; p2 = p2_; p3 = p3_;  p4 = p4_;  p5 = p5_;  p6 = p6_;
		triangle_p1 = new Triangle(p1_, p2_, p3_, id);
		triangle_p2 = new Triangle(p4_, p5_, p6, id);
		triangle_p1->material.isLight = true;
		triangle_p2->material.isLight = true;
		//triangle_p1->material.;
		//triangle_p2->material.isLight = true;


	};

	void Intersect(Ray& ray) const  {
		triangle_p1->Intersect(ray); triangle_p2->Intersect(ray);
	}
	float3 getcolor() const { return getLightColor(); };
	int getIdx() const { return idx; };
	float getLumen() const { return lumen; };

	float3 p1; float3 p2; float3 p3; float3 p4; float3 p5; float3 p6;
	Triangle* triangle_p1;
	Triangle* triangle_p2;
};

class pointLight : public light
{
public:
	pointLight() = default;
	pointLight(float lumen, float3 p1_, uint id) :light(lumen, id) {
		position = p1_;
		sphereLight = new Sphere(id, position, 0.05);
	};
	float3 position;
	Sphere* sphereLight;
	float3 getcolor() const { return getLightColor(); };
	uint8_t getIdx() const { return idx; };
	float getLumen() const { return lumen; };

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

	void Translate(float3 d) {
		for (int i = 0; i < tri.size(); i++) {
			tri[i].vertex0 += d;
			tri[i].vertex1 += d;
			tri[i].vertex2 += d;
			tri[i].centroid += d;
		}
	}

	void Scale(float3 d) {
		for (int i = 0; i < tri.size(); i++) {
			tri[i].vertex0 *= d;
			tri[i].vertex1 *= d;
			tri[i].vertex2 *= d;
			tri[i].centroid *= d;
		}
	}

	void MoveToPlane(float height) {
		float lowest = FLT_MAX;
		for (Tri triangle : tri) {
			if (triangle.vertex0.y < lowest) lowest = triangle.vertex0.y;
			if (triangle.vertex1.y < lowest) lowest = triangle.vertex1.y;
			if (triangle.vertex2.y < lowest) lowest = triangle.vertex2.y;
		}
		float diff = height - lowest;
		Translate(float3(0, diff + 0.001f, 0));
	}

	void Intersect(Ray& ray) {
		for (uint i = 0; i < tri.size(); i++) {
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
				ray.I.instPrim = Tmpl8::MakeID(meshID, objIdx, i);
			}
		}
	}
	float3 GetColor(Intersection& I) {
		uint id = GetTriangleIndex(I.instPrim);
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
			return returnval;
		}
		else {
			return material.albedo;
		}
	}
	float3 GetNormal(Intersection& I) {
		uint id = GetTriangleIndex(I.instPrim);
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
			//cout << triEx[id].N0.x << "|" << triEx[id].N0.y << "|" << triEx[id].N0.z << endl;
			//cout << id << endl;
			//float3 N = I.u * triEx[id].N1 + I.v * triEx[id].N2 + (1.0f - (I.u + I.v)) * triEx[id].N0;
			return triEx[id].N0; // Since per vertex normals dont work properly we just do face normal
			//cout << "NORMAL?: " << N.x << "|" << N.y << "|" << N.z << endl;
			//return normalize(N);
		}
		//float3 N =  -normalize(cross((tri[id].vertex1 - tri[id].vertex0), (tri[id].vertex2 - tri[id].vertex0)));
		//float3 N = float3(0, 1, 0);
		//cout << "NORMAL?: " << N.x << "|" << N.y << "|" << N.z << endl;
		//return N;
	}
	static RTXMesh* makeMesh(aiMesh* mesh, const aiScene* scene, string const& path, uint objId) {
		RTXMesh* m = new RTXMesh(objId);
		m->tri.reserve(mesh->mNumFaces);// = new Tri[mesh->mNumFaces];
		m->triEx.reserve(mesh->mNumFaces);// = new TriEx[mesh->mNumFaces];
		string directory = path.substr(0, path.find_last_of('/'));

		// We load things face by face instead of vertex by vertex. ray.I.tracers like accessing triangles compactly, I think.
		// Potentially can be changed into a list of vertices and a list of face indices, to save memory.
		for (uint i = 0; i < mesh->mNumFaces; i++) {
			Tri tri_i = Tri();
			TriEx triEx_i = TriEx();
			aiFace face = mesh->mFaces[i];
			// vertex0
			tri_i.vertex0.x = mesh->mVertices[face.mIndices[0]].x;
			tri_i.vertex0.y = mesh->mVertices[face.mIndices[0]].y;
			tri_i.vertex0.z = mesh->mVertices[face.mIndices[0]].z;
			// vertex1
			tri_i.vertex1.x = mesh->mVertices[face.mIndices[1]].x;
			tri_i.vertex1.y = mesh->mVertices[face.mIndices[1]].y;
			tri_i.vertex1.z = mesh->mVertices[face.mIndices[1]].z;
			// vertex2
			tri_i.vertex2.x = mesh->mVertices[face.mIndices[2]].x;
			tri_i.vertex2.y = mesh->mVertices[face.mIndices[2]].y;
			tri_i.vertex2.z = mesh->mVertices[face.mIndices[2]].z;
			// Calculate centroid
			tri_i.centroid = (tri_i.vertex0 + tri_i.vertex1 + tri_i.vertex2) / 3.0f;

			// Load the normals
			if (mesh->HasNormals()) {
				//// Vertex 0 Normal
				//triEx_i.N0.x = mesh->mNormals[face.mIndices[0]].x;
				//triEx_i.N0.y = mesh->mNormals[face.mIndices[0]].y;
				//triEx_i.N0.z = mesh->mNormals[face.mIndices[0]].z;
				//// Vertex 1 Normal
				//triEx_i.N1.x = mesh->mNormals[face.mIndices[1]].x;
				//triEx_i.N1.y = mesh->mNormals[face.mIndices[1]].y;
				//triEx_i.N1.z = mesh->mNormals[face.mIndices[1]].z;
				//// Vertex 2 Normal
				//triEx_i.N2.x = mesh->mNormals[face.mIndices[2]].x;
				//triEx_i.N2.y = mesh->mNormals[face.mIndices[2]].y;
				//triEx_i.N2.z = mesh->mNormals[face.mIndices[2]].z;
				//// Normalize them, because for some reason they are not.
				//triEx_i.N0 = normalize(triEx_i.N0);
				//triEx_i.N1 = normalize(triEx_i.N1);
				//triEx_i.N2 = normalize(triEx_i.N2);
				float3 N = normalize(cross((tri_i.vertex1 - tri_i.vertex0), (tri_i.vertex2 - tri_i.vertex0)));

				triEx_i.N0 = N;
				triEx_i.N1 = N;
				triEx_i.N2 = N;
			}
			else {
				float3 N = normalize(cross((tri_i.vertex1 - tri_i.vertex0), (tri_i.vertex2 - tri_i.vertex0)));
				
				triEx_i.N0 = N;
				triEx_i.N1 = N;
				triEx_i.N2 = N;
			}

			// Load the texture coordinates
			if (mesh->mTextureCoords[0]) {
				// Texture coordinates for vertex0. I think.
				triEx_i.uv0.x = mesh->mTextureCoords[0][face.mIndices[0]].x;
				triEx_i.uv0.y = mesh->mTextureCoords[0][face.mIndices[0]].y;
				// Texture coordinates for vertex1. I hope.
				triEx_i.uv1.x = mesh->mTextureCoords[0][face.mIndices[1]].x;
				triEx_i.uv1.y = mesh->mTextureCoords[0][face.mIndices[1]].y;
				// Texture coordinates for vertex2. I pray.
				triEx_i.uv2.x = mesh->mTextureCoords[0][face.mIndices[2]].x;
				triEx_i.uv2.y = mesh->mTextureCoords[0][face.mIndices[2]].y;
			}
			else {
				m->material.albedo = float3(1, 0, 0); //default red color
			}
			m->tri.push_back(tri_i);
			m->triEx.push_back(triEx_i);
		}
		//m->triCount = mesh->mNumFaces;
		m->triCount = m->tri.size();
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
			aiString str;
			material->GetTexture(aiTextureType_HEIGHT, 0, &str);
			string filename = string(str.C_Str());
			filename = directory + "/" + filename;
			m->normalMap = new Surface(filename.c_str());
			m->normalMapLoaded = true;
		}
		cout << "Loaded model with " << m->tri.size() << " triangles." << endl;
		return m;
	}
	vector<Tri> tri;
	vector<TriEx> triEx;
	int triCount;
	uint objIdx;
	mat4 M, invM;
	Material material = {
		float3(1, 0, 0), //albedo
		0.2, //specularity
		Medium::Undefined, //medium
	};
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


#if PRETTY
		// Load cat
		//loadModel("assets/wolf/Wolf.obj");
		loadModel("assets/chessboard/chessboard.obj");
		loadModel("assets/wolf/Wolf.obj");
		//meshPool[0]->material.specularity = 0.2f;
		meshPool[1]->material.mat_medium = Medium::Glass;
		meshPool[1]->Scale(float3(0.8f, 0.8f, 0.8f));
		meshPool[1]->MoveToPlane(32.f);
		meshPool[1]->material.absorption = float3(0.f, 1.0f, 1.0f);

		// Reserve object IDs for the lights
		area_lights[0] = areaLight(100, areaID, float3(32, 256, 32), float3(-32, 256, 32), float3(32, 256, -32), float3(-32, 256, -32), float3(32, 256, -32), float3(-32, 256, 32));
		spot_lights[0] = pointLight(100, float3(0, 256, 0.5), spotID);
#else
		// we store all primitives in one continuous buffer
		Sphere* sphere = new Sphere(0, float3(0), 0.5f);					// 0: bouncing ball   0.5f
		spherePool.push_back(sphere);
		Sphere* sphere2 = new Sphere(1, float3(0, -1.5f, -1.05f), 0.5f);	// 1: glass ball
		sphere2->material.albedo = float3(0.2, 0.8, 0.2);
		sphere2->material.mat_medium = Medium::Glass;
		sphere2->material.absorption = float3(0.2f, 2.0f, 4.0f);
		spherePool.push_back(sphere2);
		Cube* cube = new Cube(0, float3(0), float3(1.15f));				// 0: spinning cube
		Cube* cube2 = new Cube(1, float3(0, -1.0f, -2.0f), float3(0.8f)); // 1 Glass cube
		cube2->material.mat_medium = Medium::Glass;
		cube2->material.albedo = float3(0.2f, 0.2f, 0.2f);
		cube2->material.absorption = float3(0);
		cubePool.push_back(cube);
		cubePool.push_back(cube2);
		Triangle* triangle = new Triangle(float3(-0.9f, 0, -1), float3(0.2f, 0, -1), float3(0, 1.0f, -0.5), 0); // 0: triangle
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
		Triangle* floor_0 = new Triangle(float3(3, -3, 3), float3(3, -3, -3), float3(-3, -3, 3), 11);
		Triangle* floor_1 = new Triangle(float3(-3, -3, -3), float3(-3, -3, 3), float3(3, -3, -3), 12);
		floor_0->material.albedo = float3(0.2, 1, 0);
		floor_1->material.albedo = float3(0.2, 1, 0);

		//// area light
		//Triangle* light_1 = new Triangle(float3(1.5, 2.9, 1.5), float3(-1.5, 2.9, 1.5), float3(1.5, 2.9, -1.5), 13);
		//Triangle* light_2 = new Triangle(float3(-1.5, 2.9, -1.5), float3(1.5, 2.9, -1.5), float3(-1.5, 2.9, 1.5), 14);
		//light_1->material.albedo = float3(1.0, 1.0, 1.0);
		//light_2->material.albedo = float3(1.0, 1.0, 1.0);


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

		// Reserve object IDs for the lights
		area_lights[0] = areaLight(10, areaID, float3(1.5, 2.9, 1.5), float3(-1.5, 2.9, 1.5), float3(1.5, 2.9, -1.5), float3(-1.5, 2.9, -1.5), float3(1.5, 2.9, -1.5), float3(-1.5, 2.9, 1.5));
		spot_lights[0] = pointLight(10, float3(0, 1.5, 0.5), spotID);
#endif


		

		//lightSphere = L[0];



		//trianglePool.push_back(light_1);
		//trianglePool.push_back(light_2);

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
		if(cubePool.size() > 0) cubePool[0]->M = M2, cubePool[0]->invM = M2.FastInvertedTransformNoScale();
		// sphere animation: bounce
		float tm = 1 - sqrf( fmodf( animTime, 2.0f ) - 1 );
		if(spherePool.size() > 0 )spherePool[0]->pos = float3(-1.4f, -0.5f + tm, 2);
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

	float3 GetLightPos(int i = 0) const
	{
		// light point position is the middle of the swinging quad
		//float3 corner1 = TransformPosition( float3( -0.5f, 0, -0.5f ), quad.T  );
		//float3 corner2 = TransformPosition( float3( 0.5f, 0, 0.5f ) , quad.T );
		//return (corner1 + corner2) * 0.5f - float3( 0, 0.01f, 1 );
		return spot_lights[i].position;
	}

	Material getMaterial(uint objIdx) const
	{
		uint Id = GetObjectIndex(objIdx);
		if (Id == spotID) return spot_lights[0].getMaterial();					// if we add multiple light souce this would need to change
		else if (objIdx == areaID) return area_lights[0].getMaterial();	// if we add multiple light souce this would need to change
		uint type = GetObjectType(objIdx);
		if (type == sphereID) return spherePool[Id]->material;
		if (type == cubeID) return cubePool[Id]->material;
		if (type == triangleID) return trianglePool[Id]->material;
		if (type == meshID) return meshPool[Id]->material;
		//uint objId = objIdx >> 20;
		//if (objId == 0) throw exception("There's no material for nothing"); // or perhaps we should just crash
		////if (objIdx == 0) return quad.material;
		////if (objIdx == 0) return lightSphere.position; 
		//else if (objIdx == 1) return sphere.material;
		//else if (objIdx == 2) return sphere2.material;
		//else if (objIdx == 3) return cube.material;
		//else if (objIdx == 4) return cube2.material;
		//else if (objIdx == 5) return spot_lights[0].getcolor();					// if we add multiple light souce this would need to change
		//else if (objIdx == 6) return area_lights[0].getcolor();	// if we add multiple light souce this would need to change

		//else if (objIdx >= 10) return trianglePool[objIdx - 10]->material;
		//cout << objIdx << endl;
		throw exception("ID not known");
	}

	float3 GetLightColor(int i =0) const
	{
		if (isSpotLight) {
			return spot_lights[i].getcolor();}
		else {
			return area_lights[i].getcolor();}
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

	float3 directIllumination(Intersection& I, float3 intersection, float3 norm) {
		uint objType = GetObjectType(I.instPrim);
		uint objIdx = GetObjectIndex(I.instPrim);
		float3 obj_Color = GetColor(I);
		/*if (objIdx == 1000) { 
			return spot_lights[0].getcolor(); }*/
		
		float3 color = (0, 0, 0);
		for (int i = 0; i < size(spot_lights); i++) {
			float3 dir_light = normalize(spot_lights[i].position - intersection);
			float dot_p = dot(dir_light, norm);
			float intensity = 1.0f / pow(length(spot_lights[i].position - intersection), 2);

			// -----------------------------------------------------------
			// regular
			// -----------------------------------------------------------
			color += intensity * spot_lights[i].getcolor() * obj_Color * maxFloat( dot_p, 0.0f); //
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
		if (isSpotLight) {
			spot_lights[0].sphereLight->Intersect(ray);

		}else{
			area_lights[0].Intersect(ray);
		}
		for (int i = 0; i < spherePool.size(); i++)	spherePool[i]->Intersect(ray);
		for (int i = 0; i < cubePool.size(); i++) cubePool[i]->Intersect(ray);
		for (Triangle* tri : trianglePool) tri->Intersect(ray);
		for (RTXMesh* mesh : meshPool) 	mesh->Intersect(ray);
		//if (ray.objIdx >= 10) cout << "Triangle hit: " << ray.objIdx << endl;
	}
	bool IsOccluded( Ray& ray ) const
	{
		float rayLength = ray.I.t;
		// skip planes: it is not possible for the walls to occlude anything
		//quad.Intersect( ray );
		for (int i = 0; i < spherePool.size(); i++)	spherePool[i]->Intersect(ray);
		for (int i = 0; i < cubePool.size(); i++) cubePool[i]->Intersect(ray);
		for (Triangle* tri : trianglePool) tri->Intersect(ray);
		for (RTXMesh* mesh : meshPool) 	mesh->Intersect(ray);
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
		uint typeId = GetObjectType(Inters.instPrim);
		uint objId = GetObjectIndex(Inters.instPrim);
		if (objId == areaID || objId == spotID) throw exception("Normal of light not implemented.");
		if (typeId == 6) throw exception("NOT ALLOWED");//return float3( 0 ); // or perhaps we should just crash
		float3 N;
		//if (objIdx == 0) N = quad.GetNormal(I);
		//if (objId == 1) N = lightSphere.GetNormal(I); // check objIDX to 2 from 0
		//if (objId == 2) N = sphere.GetNormal(I);
		//else if (objId == 3) N = sphere2.GetNormal(I);
		//else if (objId == 4) N = cube.GetNormal(I);
		//else if (objId == 5) N = cube2.GetNormal(I);
		if (typeId == sphereID) N = spherePool[objId]->GetNormal(I);
		if (typeId == cubeID) N = cubePool[objId]->GetNormal(I);
		if (typeId == triangleID) N = trianglePool[objId]->GetNormal();
		if (typeId == meshID) N = meshPool[objId]->GetNormal(Inters);
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
		uint typeId = GetObjectType(I.instPrim);
		uint Id = GetObjectIndex(I.instPrim);
		if (Id == areaID) return area_lights[0].getcolor();
		if (Id == spotID) return spot_lights[0].getcolor();
		if (typeId == sphereID) return spherePool[Id]->material.albedo;
		if (typeId == cubeID) return cubePool[Id]->material.albedo;
		if (typeId == triangleID) return trianglePool[Id]->material.albedo;
		if (typeId == meshID) return meshPool[Id]->GetColor(I);
		if (typeId == skyBoxID) {
			cout << "Did not expect to reach this code." << endl;
			return getSkyBox(float3(0, 1, 0));
		}
		throw exception("Unidentified type id in scene.GetColor()");
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
	}

	void processNode(aiNode* node, const aiScene* scene, const char* path) {
		// process parent node
		for (int i = 0; i < node->mNumMeshes; i++) {
			int id = meshPool.size(); //TODO: REMOVE THE START AT 2
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
	const int spotID = 500;
	const int areaID = 501;
	//Quad quad;
	//Light lights[1];
	pointLight spot_lights[1];
	areaLight area_lights[1];

	pointLight lightSphere;

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
	//int numLightSouces=1;								// should be initilize better in the future


	bool isSpotLight= sendWhittedCONFIG;								// we will need to find one variable that connect both is spot to whitted
};



}