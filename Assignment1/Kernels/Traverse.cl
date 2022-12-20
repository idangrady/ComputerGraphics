#include <utils.cl>



__attribute__((always_inline)) void traverse_tree(global buffer &nodes,global buffer &nodes_idx int max_depth)
	{

	  int stack_size = 0;
	  __local BVHNode stack[SCRWIDTH2*SCRHEIGHT2];				// stack for traversing the tree
	    stack[stack_size++] = bvhnodes[0];						// init with the first bvhnode
																	
		while(stack_size>0)										// iterate over all values 
		{
			// pop the head bvhnode
			BVHNode node = [--stack_size];						// check if init is currect

			// if intersect and the count >0
			if (!IntersectAABB(ray, node.aabbMin, node.aabbMax)) return;
			if (node.triCount>0)
			{
				for (uint i = 0; i < node.triCount; i++) 
				{ 
					nodes[nodes_idx[node.leftFirst + i]]->Intersect(ray);
				}
			}
			else
			{
				stack[stack_size++]= nodes[nodes_idx[node.leftFirst]];
				stack[stack_size++]= nodes[nodes_idx[node.leftFirst+1]] ;
			}
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

