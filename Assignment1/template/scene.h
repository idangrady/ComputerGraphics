#pragma once
#include <map>
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
#define PLANE_X(o,i) {if((t=-(ray.O.x+o)*ray.rD.x)<ray.t)ray.t=t,ray.objIdx=i;}
#define PLANE_Y(o,i) {if((t=-(ray.O.y+o)*ray.rD.y)<ray.t)ray.t=t,ray.objIdx=i;}
#define PLANE_Z(o,i) {if((t=-(ray.O.z+o)*ray.rD.z)<ray.t)ray.t=t,ray.objIdx=i;}

enum class Medium {
	Undefined = -1,
	Air = 0,
	Glass = 1,
};

struct Light {
	Light() = default;
	Light(float3 p, float i) : position(p), intensity(i) {}
	float3 position;
	float3 color;
	float intensity;
};

struct Material
{
	Material() = default;
	Material(float3 a, float d, float s, float dif, Medium m) : albedo(a), diffuse(d), specular(s), diffractive(dif), mat_medium(m) {}
	float3 albedo = (0.9, 0.9, 0.9); // color material
	float diffuse = 1; // material 
	float specular = 0;
	float diffractive = 0;
	Medium mat_medium{ Medium::Undefined };
};

namespace Tmpl8 {

__declspec(align(64)) class Ray
{
public:
	Ray() = default;
	Ray(const Ray&) {}
	Ray(float3 origin, float3 direction, float distance = 1e34f, int depth = 0 )
	{
		O = origin, D = direction, t = distance;
		// calculate reciprocal ray direction for triangles and AABBs
		rD = float3( 1 / D.x, 1 / D.y, 1 / D.z );
		depthidx = depth;
	#ifdef SPEEDTRIX
		d0 = d1 = d2 = 0;
	#endif
	}
	float3 IntersectionPoint() { return O + t * D; }
	Ray Reflect(float3 intersec,float3 norm, int idx) { 
		// create a secondary ray
		// intersection = origin, norm = norm from intersected object
		float3 dir = this->D;
		Ray ray = Ray(intersec + (0.0002 * norm), dir - 2 * (dot(dir, norm) * norm), 1e34f, idx + 1);
		//ray.objIdx = this->objIdx;
		//ray.inside = this->inside;
		return ray;
	}
	// ray data
#ifndef SPEEDTRIX
	float3 O, D, rD;
#else
	union { struct { float3 O; float d0; }; __m128 O4; };
	union { struct { float3 D; float d1; }; __m128 D4; };
	union { struct { float3 rD; float d2; }; __m128 rD4; };
#endif
	float t = 1e34f;
	int objIdx = -1;
	bool inside = false; // true when in medium
	int depthidx;

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
	Sphere( int idx, float3 p, float r ) : 
		pos(p), r2(r* r), invr(1 / r), objIdx(idx) {
	}
	Sphere(int idx, float3 p, float r, Material mat) : Sphere(idx, p, r) {
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
		if (t < ray.t && t > 0)
		{
			ray.t = t, ray.objIdx = objIdx;
			return;
		}
		t = d - b;
		if (t < ray.t && t > 0)
		{
			ray.t = t, ray.objIdx = objIdx;
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
	int objIdx = -1;

	float3 color = float3(0, 0, 1);
	static int id; 
	Material material = {
		float3(0, 0, 1), //albedo
		1.0, //diffuse
		0.0, //specular
		0.0, //diffraction
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
	Plane(int idx, float3 normal, float dist, float3 col = (1, 0.5, 0.5)) : N(normal), d(dist), objIdx(idx){ color = col; }
	void Intersect( Ray& ray ) const
	{
		float t = -(dot( ray.O, this->N ) + this->d) / (dot( ray.D, this->N ));
		if (t < ray.t && t > 0) ray.t = t, ray.objIdx = objIdx;
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
	int objIdx = -1;
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
	Cube( int idx, float3 pos, float3 size, mat4 transform = mat4::Identity() )
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
			if (tmin < ray.t) ray.t = tmin, ray.objIdx = objIdx; 
		}
		else if (tmax > 0)
		{
			if (tmax < ray.t) ray.t = tmax, ray.objIdx = objIdx;
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
	int objIdx = -1;
	Material material = {
		float3(0, 1, 0), //albedo
		0.8, //diffuse
		0.2, //specular
		0.0, //diffraction
		Medium::Undefined, //medium
	};
	float3 color = float3(0, 1, 0);;
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
	Triangle(float3 v0, float3 v1, float3 v2, int id) : vertex0(v0), vertex1(v1), vertex2(v2) {
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
		if (t > 0.0001f && t < ray.t) {
			ray.t = t;	
			ray.objIdx = objIdx;
		}
	}
	float3 GetNormal() const 
	{
		return normal;
		//return float3(0, 0, -1);
	}
	union {
		struct { float3 vertex0, vertex1, vertex2; };
		float3 cell[3];
	};
	int objIdx;
	float3& operator [] (const bool n) { return cell[n]; }
	float3 centroid;
	float3 normal;
	Material material = {
		float3(1, 0, 0), //albedo
		1.0, //diffuse
		0.0, //specular
		0.0, //diffraction
		Medium::Undefined, //medium
	};
};

// -----------------------------------------------------------
// Quad primitive
// Oriented quad, intended to be used as a light source.
// -----------------------------------------------------------
class Quad
{
public:
	Quad() = default;
	Quad( int idx, float s, mat4 transform = mat4::Identity() )
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
		if (t < ray.t && t > 0)
		{
			float3 I = O + t * D;
			if (I.x > -size && I.x < size && I.z > -size && I.z < size)
				ray.t = t, ray.objIdx = objIdx;
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
		1.0, //diffuse
		0.0, //specular
		0.0, //diffraction
		Medium::Undefined, //medium
	};
	float size;
	mat4 T, invT;
	int objIdx = -1;
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
		lights[0] = Light(float3(0, 1.5, 0.5), 24);
		// Precalc indices
		refractiveIndices[Medium::Glass][Medium::Air] = 1.0 / 1.52;
		refractiveIndices[Medium::Air][Medium::Glass] = 1.52;
		// Precalc some materials
		Material mirror = Material(float3(0, 0.2, 0), 0.0, 1.0, 0.0, Medium::Undefined);
		materialMap["Mirror"] = mirror;


		// we store all primitives in one continuous buffer
		//quad = Quad( 0, 1 );									// 0: light source
		lightSphere = Sphere(0, lights[0].position, 0.05);
		lightSphere.material.albedo = (24, 24, 24);
		sphere = Sphere( 1, float3( 0 ), 0.5f);					// 1: bouncing ball   0.5f
		//sphere2 = Sphere( 2, float3( 0, 2.5f, -3.07f ), 8 );	// 2: rounded corners
		cube = Cube(3, float3(0), float3(1.15f));				// 3: cube 		cube = Cube( 3, float3( 0 ), float3( 1.15f ) );		
		//plane[0] = Plane( 4, float3( 1, 0, 0 ), 3 , float3(1,0.5,1));			// 4: left wall
		//plane[1] = Plane( 5, float3( -1, 0, 0 ), 2.99f, float3(0, 0.23, 0.23));		// 5: right wall
		//plane[2] = Plane( 6, float3( 0, 1, 0 ), 1, float3(1, 1, 0.23));			// 6: floor
		//plane[3] = Plane( 7, float3( 0, -1, 0 ), 2, float3(0.23, 0.23, 1));			// 7: ceiling
		//plane[4] = Plane( 8, float3( 0, 0, 1 ), 3, float3(0, 0.23, 1));			// 8: front wall
		//plane[5] = Plane( 9, float3( 0, 0, -1 ), 3.99f, float3(0.23, 0, 1));		// 9: back wall
		
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
		//cout << "Triangle 0: " << triangle->objIdx << endl;
		trianglePool.push_back(triangle);
		//cout << "Triangle wall L 0: " << wall_l0->objIdx << endl;
		trianglePool.push_back(wall_l0);
		//cout << "Triangle wall L 1: " << wall_l1->objIdx << endl;
		trianglePool.push_back(wall_l1);
		//cout << "Triangle wall R 0: " << wall_r0->objIdx << endl;
		trianglePool.push_back(wall_r0);
		//cout << "Triangle wall R 0: " << wall_r1->objIdx << endl;
		trianglePool.push_back(wall_r1);
		trianglePool.push_back(wall_r_backdrop0);
		trianglePool.push_back(wall_r_backdrop1);
		trianglePool.push_back(ceiling_0);
		trianglePool.push_back(ceiling_1);
		trianglePool.push_back(wall_b0);
		trianglePool.push_back(wall_b1);
		trianglePool.push_back(floor_0);
		trianglePool.push_back(floor_1);
		//cout << "Sanity checks: " << endl;
		//cout << "10 == " << trianglePool[0]->objIdx << endl;
		//cout << "11 == " << trianglePool[1]->objIdx << endl;
		//cout << "12 == " << trianglePool[2]->objIdx << endl;
		//cout << "13 == " << trianglePool[3]->objIdx << endl;
		//cout << "14 == " << trianglePool[4]->objIdx << endl;
		SetTime( 0 );
		// Note: once we have triangle support we should get rid of the class
		// hierarchy: virtuals reduce performance somewhat.
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
	float3 GetLightPos() const
	{
		// light point position is the middle of the swinging quad
		//float3 corner1 = TransformPosition( float3( -0.5f, 0, -0.5f ), quad.T  );
		//float3 corner2 = TransformPosition( float3( 0.5f, 0, 0.5f ) , quad.T );
		//return (corner1 + corner2) * 0.5f - float3( 0, 0.01f, 1 );
		return lights[0].position;
	}

	Material getMaterial(int objIdx) const
	{
		if (objIdx == -1) throw exception("There's no material for nothing"); // or perhaps we should just crash
		//if (objIdx == 0) return quad.material;
		if (objIdx == 0) return lightSphere.material;
		else if (objIdx == 1) return sphere.material;
		//else if (objIdx == 2) return sphere2.material;
		else if (objIdx == 3) return cube.material;
		else if (objIdx >= 10) return trianglePool[objIdx - 10]->material;
		throw exception("ID not known");
	}

	float3 GetLightColor() const
	{
		return lights[0].color; //return float3( 24, 24, 22 );
	}
	float3 directIllumination(int objIdx, float3 intersection, float3 norm, float3 albedo) {
		if (objIdx == 0) return lightSphere.material.albedo;
		float3 color = (0, 0, 0);
		for (int i = 0; i < numLightSouces; i++) { // would change once we add more lights
			float3 dir_light = normalize(lights[i].position - intersection);
			float dot_p = dot(dir_light, norm);
			float intensity = lights[i].intensity / pow(length(lights[i].position - intersection), 2);

			// -----------------------------------------------------------
			// regular
			// -----------------------------------------------------------
			color += intensity * albedo * maxFloat( dot_p, 0.0f); //
		
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
		//if (ray.D.x < 0) PLANE_X( 3, 4 ) else PLANE_X( -2.99f, 5 );
		//if (ray.D.y < 0) PLANE_Y( 1, 6 ) else PLANE_Y( -2, 7 );
		//if (ray.D.z < 0) PLANE_Z( 3, 8 ) else PLANE_Z( -3.99f, 9 );
		//quad.Intersect( ray );
		lightSphere.Intersect(ray);
		sphere.Intersect( ray );
		//sphere2.Intersect( ray );
		cube.Intersect( ray );
		for (int i = 0; i < (int)trianglePool.size(); i++) {
			trianglePool[i]->Intersect(ray);
		}
		//if (ray.objIdx >= 10) cout << "Triangle hit: " << ray.objIdx << endl;
	}
	bool IsOccluded( Ray& ray ) const
	{
		float rayLength = ray.t;
		// skip planes: it is not possible for the walls to occlude anything
		//quad.Intersect( ray );
		sphere.Intersect( ray );
		//sphere2.Intersect( ray );
		cube.Intersect( ray );
		for (int i = 0; i < (int)trianglePool.size(); i++) {
			trianglePool[i]->Intersect(ray);
		}
		return ray.t < rayLength;
		// technically this is wasteful: 
		// - we potentially search beyond rayLength
		// - we store objIdx and t when we just need a yes/no
		// - we don't 'early out' after the first occlusion
	}
	float3 GetNormal( int objIdx, float3 I, float3 wo ) const
	{
		// we get the normal after finding the nearest intersection:
		// this way we prevent calculating it multiple times.

		if (objIdx == -1) return float3( 0 ); // or perhaps we should just crash
		float3 N;
		//if (objIdx == 0) N = quad.GetNormal(I);
		if (objIdx == 0) N = lightSphere.GetNormal(I);
		if (objIdx == 1) N = sphere.GetNormal(I);
		//else if (objIdx == 2) N = sphere2.GetNormal(I);
		else if (objIdx == 3) N = cube.GetNormal(I);
		else if (objIdx >= 10) N = trianglePool[objIdx - 10]->GetNormal();
		//else 
		//{
		//	// faster to handle the 6 planes without a call to GetNormal
		//	N = float3( 0 );
		//	N[(objIdx - 4) / 2] = 1 - 2 * (float)(objIdx & 1);
		//}
		if (dot( N, wo ) > 0) N = -N; // hit backside / inside
		return N;
	}
	__declspec(align(64)) // start a new cacheline here
	float animTime = 0;
	const float ambient = 0.005;
	//Quad quad;
	Light lights[1];
	Sphere sphere;
	Sphere lightSphere;
	//Sphere sphere2;
	Cube cube;
	//Plane plane[6];
	static inline vector<Triangle*> trianglePool;
	static inline map<string, Material> materialMap;
	map<Medium, map<Medium, float>> refractiveIndices;
	int numLightSouces=1; // should be initilize better in the future
};

}