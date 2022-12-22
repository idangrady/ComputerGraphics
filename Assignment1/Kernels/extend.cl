#include "template/common.h"
#include "Kernels/utils.cl"

inline uint makeId(uint type, uint id, uint tri){
    return (type << 30) + (id << 20) + (tri);
}

void IntersectTri(Ray* ray, __constant Triangle* tri) 
{
    const float3 edge1 = tri->vertex1.xyz - tri->vertex0.xyz;
    const float3 edge2 = tri->vertex2.xyz - tri->vertex0.xyz;
    const float3 h = cross(ray->D.xyz, edge2);
    const float a = dot(edge1, h);
    if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
    const float f = 1 / a;
    const float3 s = ray->O.xyz - tri->vertex0.xyz;
    const float u = f * dot(s, h);
    if (u < 0 || u > 1) return;
    const float3 q = cross(s, edge1);
    const float v = f * dot(ray->D.xyz, q);
    if (v < 0 || u + v > 1) return;
    const float t = f * dot(edge2, q);
    if (t > 0.0001f && t < ray->t){
        ray->t = t;
        ray->uv = (float2)(u, v);
        ray->primIdx = makeId(1, 0, tri->id);
    }
}

void IntersectSphere(Ray* ray, __constant Triangle* tri, int r) 
{
		{
			float3 pos = tri->vertex0.xyz;
			float3 oc = ray->O.xyz - tri->vertex0.xyz;
			float b = dot( oc, ray->D.xyz );
			float c = dot( oc, oc ) - 5;
			float t, d = b * b - c;
			if (d <= 0) return;
			d = sqrt( d ), t = -b - d;
			if (t < ray->t && t > 0)
			{
				ray->t = t; ray->primIdx = makeId(1, 0, tri->id);
				return;
			}
			t = d - b;
			if (t < ray->t && t > 0)
			{
               ray->t = t;
               ray->primIdx = makeId(1, 0, tri->id);
			   return;
			}
		}
}


__kernel void extend(__global Ray* rays, __constant Triangle* triangles,__constant BVHNode* bvhnodes,__constant int * primitives_idx, int triangleCount)    
{
	int threadIdx = get_global_id(0);
    Ray* ray = &rays[threadIdx];
	int stack_size = 0;

	__local BVHNode stack[10];									// stack for traversing the tree
	stack[stack_size++] = bvhnodes[0];							// init with the first bvhnode
																	
	while(stack_size>0 && stack_size<10)											// iterate over all values -> would run at most SCRWIDTH2*SCRHEIGHT2
	{
		BVHNode node = stack[--stack_size];						// check if init is currect

		float tx1 = (node.aabbMin.x - ray->O.x) / ray->D.x, tx2 = (node.aabbMax.x - ray->O.x) / ray->D.x;
		float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
		float ty1 = (node.aabbMin.y - ray->O.y) / ray->D.y, ty2 = (node.aabbMax.y - ray->O.y) / ray->D.y;
		tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
		float tz1 = (node.aabbMin.z - ray->O.z) / ray->D.z, tz2 = (node.aabbMax.z - ray->O.z) / ray->D.z;
		tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));


		if (tmax >= tmin && tmin < ray->t && tmax > 0) {}
		else{
		int i =0;

		if (node.triCount>0)										// if root
			{
			int i=0;
			if(triangles[primitives_idx[node.leftFirst + i]].objId==0) IntersectTri(ray, &triangles[primitives_idx[node.leftFirst + i]]);
			else{IntersectSphere(ray, &triangles[primitives_idx[node.leftFirst + i]], 2);}

				for (uint i = 0; i < node.triCount; i++) 
				{ 
					if(triangles[primitives_idx[node.leftFirst + i]].objId==0) IntersectTri(ray, &triangles[primitives_idx[node.leftFirst + i]]);
					else{IntersectSphere(ray, &triangles[primitives_idx[node.leftFirst + i]], 2);}
				}
			}
			else												// add to the stack 
				{
					stack[stack_size++]= bvhnodes[primitives_idx[node.leftFirst]];
					stack[stack_size++]= bvhnodes[primitives_idx[node.leftFirst+1]] ;
				}
			}
	}

}




