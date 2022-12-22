
				if (false) // sphere we init all the sphere currently with the same size
			{
				float3 pos = *points[0];
				float3 oc = ray.O - *points[0];
				float b = dot(oc, ray.D);
				float c = dot(oc, oc) - r;
				float t, d = b * b - c;
				if (d <= 0) return;
				d = sqrtf(d), t = -b - d;
				if (t < ray.I.t && t > 0)
				{
					ray.I.t = t, ray.I.instPrim = MakeID(sphereID, 1, 0);
					return;
				}
				t = d - b;
				if (t < ray.I.t && t > 0)
				{
					ray.I.t = t, ray.I.instPrim = MakeID(sphereID, 1, 0);
					return;
				}
			}
		}

		float3 GetNormal(const float3 I) const
		{
			if (object == 1) return (I - *points[0]) * (1 / r);													// sphere
			if (object == 3) return normalize(cross((*points[1] - *points[0]), (*points[2] - *points[0])));		// trig
		}
		float3 getCentroid()
		{
			if (object == 1) return (*points[0]);																// sphere
			if (object == 3) return (*points[1] + *points[0] + *points[2]) / 3;									// trig
		}






		void traverse_tree(__global BVHNode *nodes,__global int *nodes_idx)
	{
	  int stack_size = 0;
	  __local BVHNode stack[SCRWIDTH*SCRHEIGHT];			// stack for traversing the tree
	    stack[stack_size++] = nodes[0];						// init with the first bvhnode
																	
		while(stack_size>0)										// iterate over all values -> would run at most SCRWIDTH2*SCRHEIGHT2
		{
			// pop the head bvhnode
			BVHNode node = [--stack_size];						// check if init is currect

			// if intersect and the count >0
			//if (!IntersectAABB(ray, node->aabbMin, node->aabbMax)) return;
			if (node->triCount>0)								// if root
			{
				for (uint i = 0; i < node->triCount; i++) 
				{ 
					nodes[nodes_idx[node->leftFirst + i]]->IntersectTri(ray);
				}
			}
			else												// add to the stack 
			{
				stack[stack_size++]= nodes[nodes_idx[node->leftFirst]];
				stack[stack_size++]= nodes[nodes_idx[node->leftFirst+1]] ;
			}
		}

	}



	
__kernel void extend(__global Ray* rays, __constant Triangle* triangles,__constant BVHNode* bvhnodes,__constant int * primitives_idx, int triangleCount)    
{
	int threadIdx = get_global_id(0);
    Ray* ray = &rays[threadIdx];
	int stack_size = 0;

	__local BVHNode stack[200];									// stack for traversing the tree
	stack[stack_size++] = bvhnodes[0];							// init with the first bvhnode
																	
	while(stack_size>0)											// iterate over all values -> would run at most SCRWIDTH2*SCRHEIGHT2
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
			if (node.triCount>0)										// if root
			{
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












