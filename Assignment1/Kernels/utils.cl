
#define SCRWIDTH2	1280
#define SCRHEIGHT2	720  


const uint sphereID = 0;
const uint cubeID = 1;
const uint triangleID = 2;
const uint meshID = 3;
const uint planeID = 4;
const uint lightID = 5;
const uint skyBoxID = 6;



struct Intersection
{
	Intersection() = default;
	float t = 1e34f;								// intersection distance along ray
	float u, v;										// barycentric coordinates of the intersection
	uint instPrim = MakeID(skyBoxID, 0, 0);			// Type indedx (3 bit), instance index (9 bit) and primitive index (20 bit)
};


typedef struct Ray {
 float3 O;
  float3 D;
  Intersection I;
  int depthidx =0;
} Ray;


static inline uint MakeID(uint type, uint id, uint tri) {
	return (type << 29) + (id << 20) + (tri);
}

struct primitive
{
int object;										// this is to 	
vector<float3> points;							// This will be used to initilize which object it is -> if size =3-> rechtangle. If 2 cube.. etc
int r = 2;										// we should make this better 

void intersect(Ray ray)
{
	if(points.size()==3) // triangle
	{
		const float3 edge1 = points[1] - points[0];
		const float3 edge2 = points[2] - points[1];
		const float3 h = cross(ray.D, edge2);
		const float a = dot(edge1, h);
		if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
		const float f = 1 / a;
		const float3 s = ray.O - points[0];
		const float u = f * dot(s, h);
		if (u < 0 || u > 1) return;
		const float3 q = cross(s, edge1);
		const float v = f * dot(ray.D, q);
		if (v < 0 || u + v > 1) return;
		const float t = f * dot(edge2, q);
		if (t > 0.0001f && t < ray.I.t) 
		{
			ray.I.t = t;
			ray.I.u = u;
			ray.I.v = v;
			ray.I.instPrim = Tmpl8::MakeID(triangleID, objIdx, 0);;
		}
	}
	if(points.size()==2) // cube
		{
		// currently not supported
		}

	if(points.size()==1) // sphere we init all the sphere currently with the same size
		{
			float3 pos = points[0];
			float3 oc = ray.O - points[0];
			float b = dot( oc, ray.D );
			float c = dot( oc, oc ) - r;
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
	}

	float3 GetNormal( const float3 I ) const
	{
		if (points.size()==1) return (I - points[0]) * (1/r);												// sphere
		if(points.size()==3) return normalize(cross((points[1] - points[0]), (points[2]- points[0])));		// trig
	}
	float3 getCentroid()
	{
		if (points.size()==1) return (points[0]);															// sphere
		if(points.size()==3) return (points[1]+points[0]+points[2])/3;										// trig
	}
}

struct BVHNode
{
	float3 aabbMin, aabbMax;
	uint leftFirst, triCount;
};