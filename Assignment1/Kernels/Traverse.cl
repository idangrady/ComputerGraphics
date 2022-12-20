#include <utils.cl>



__attribute__((always_inline)) void traverse_tree( buffer B, int max_depth)
	{
	int root = 0;

	BVHNode& node = *B[nodeIdx];
	int i =0;
	while(!node.is_leaf() || i>max_depth){
	
	if (!IntersectAABB(ray, node.aabbMin, node.aabbMax)) return;
	if (node.triCount>0)
	{
		for (uint i = 0; i < node.triCount; i++) { 
			arrPrimitive[arrPrimitiveIdx[node.leftFirst + i]]->Intersect(ray);
		}
	}
	else
	{
		IntersectBVH(ray, node.leftFirst);
		IntersectBVH(ray, node.leftFirst + 1); 
	}
	

	i++;
	}
	}


	// make primary rays
__kernel void makePrimaryRays(read_only MakePrimaryRays) 
{

for each buffered ray r
{
cur_ray = MakePrimaryRays[i]
}
}



O,D,dist,primIdx= MakePrimaryRays[i]
I = IntersectionPoint( O, D, dist)
N = PrimNormal( primIdx, I )
if (NEE) {
si= atomicInc( shadowRayIdx)
shadowBuffer[si] = ShadowRay( … )
}
if (bounce) {
ei= atomicInc( extensionRayIdx)
newRayBuffer[ei] = ExtensionRay( … )
