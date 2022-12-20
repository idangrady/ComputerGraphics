

#define SCRWIDTH2	1280
#define SCRHEIGHT2	720


typedef struct Ray {
 float3 O;
  float3 D;
  Intersection I;
  int depthidx =0;
} Ray;


__attribute__((always_inline)) Ray GetPrimaryRay( const int x, const int y, __global float3* camPos )
	{
		// calculate pixel position on virtual screen plane
		const float u = (float)x * (1.0f / SCRWIDTH);
		const float v = (float)y * (1.0f / SCRHEIGHT);
		const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
		Ray ray;
		ray.O = *camPos;
		ray.D = normalize( P - ray.O );
		return Ray;
	}

	// make primary rays
__kernel void makePrimaryRays(write_only MakePrimaryRays,__global float3* camPos ) 
{
	int threadIdx = get_global_id(0);
	int x = threadIdx % SCRWIDTH;
	int y = threadIdx / SCRWIDTH;
	Ray ray = GetPrimaryRay(x, y, camPos); 
	MakePrimaryRays[threadIdx] =ray;
}


