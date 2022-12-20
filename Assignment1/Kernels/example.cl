

#define SCRWIDTH	1280
#define SCRHEIGHT	720

// Generate a rainbow hue from HSV
 __attribute__((always_inline)) float3 getColor(int i)
{
	int c = i % (SCRWIDTH + 1);
	float hue = (1.0f / SCRWIDTH) * (float)c * 360.0f;
	float X = 1.0f - fabs(fmod(hue / 60.0f, 2.0f) - 1.0f);
	if(hue <= 60.0f) return (float3)(1, X, 0);
	if(hue <= 120.0f) return (float3)(X, 1, 0);
	if(hue <= 180.0f) return (float3)(0, 1, X);
	if(hue <= 240.0f) return (float3)(0, X, 1);
	if(hue <= 300.0f) return (float3)(X, 0, 1);
	return (float3)(1, 0, X);
}


// Make a nice moving rainbow hue
__kernel void renderToScreen(write_only image2d_t target, __constant int* frame)
{
	int threadIdx = get_global_id(0);
	int x = threadIdx % SCRWIDTH;
	int y = threadIdx / SCRWIDTH;
	int num = (x + frame[0]);
    write_imagef( target, (int2)(x, y), (float4)(getColor(num), 1 ) );
}
