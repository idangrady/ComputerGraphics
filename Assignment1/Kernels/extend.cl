#include "Kernels/utils.cl"

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
    if (t > 0.0001f && t < ray->rD_t.w){
        ray->rD_t.w = t;
        ray->uv = (float2)(u, v);
        ray->primIdx = makeId(1, 0, tri->id);
    }
}

__kernel void extend(__global Ray* rays, __constant Triangle* triangles,__constant int*arrPrimitivesIdx , __constant BVHNode* bvhnodes, int triangleCount) 
 {
	int threadIdx = get_global_id(0);
    Ray* ray = &rays[threadIdx];
 
 
 	int stack_size = 0;

	__local BVHNode stack[PrimitiveCount*2-1];									// stack for traversing the tree
	stack[stack_size++] = bvhnodes[0];	
	BVHNode node = stack[0];						// check if init is currect


    while(stack_size>0  )											// iterate over all values -> would run at most SCRWIDTH2*SCRHEIGHT2
	{
		BVHNode *node = &stack[--stack_size];						// check if init is currect

		float tx1 = (node->aabbMin.x - ray->O.x) / ray->D.x, tx2 = (node->aabbMax.x - ray->O.x) / ray->D.x;
		float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
		float ty1 = (node->aabbMin.y - ray->O.y) / ray->D.y, ty2 = (node->aabbMax.y - ray->O.y) / ray->D.y;
		tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
		float tz1 = (node->aabbMin.z - ray->O.z) / ray->D.z, tz2 = (node->aabbMax.z - ray->O.z) / ray->D.z;
		tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
		
		if (tmax >= tmin && tmin < ray->rD_t.w && tmax > 0) {} 
		else{

		if (node->triCount>0)										// if root
			{
			    for(int i = 0; i < triangleCount; i++){
				IntersectTri(ray, &triangles[i]);
				}

			}
			else													// add to the stack 
				{
					stack[stack_size++]= bvhnodes[arrPrimitivesIdx[node->leftFirst]];
					stack[stack_size++]= bvhnodes[arrPrimitivesIdx[node->leftFirst+1]] ;
				}

		}
	}

}
