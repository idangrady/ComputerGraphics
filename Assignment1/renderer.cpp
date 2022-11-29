#include "precomp.h"
// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------

int i = 1;
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);

	POINT cursorPosition;
	GetCursorPos(&cursorPosition);
	mousePos = int2(cursorPosition.x, cursorPosition.y);
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace(Ray& ray)
{


	scene.FindNearest(ray);
	if (ray.objIdx == -1) return float3(0, 0, 0.2); // or a fancy sky color

	float3 color(0, 0, 0);
	float3 I = ray.O + ray.t * ray.D;
	float3 N = scene.GetNormal(ray.objIdx, I, ray.D);
	Material& m = scene.getMaterial(ray.objIdx);

	float s = m.specular;
	float d = m.diffuse; //or: float d = 1 - s; 

	//Ray secondary_ray = ray.reflect(I, N, ray.depthidx);
	float3 offset_O = I + (0.0002 * N);
	float3 distance_to_light = scene.lights[0].position - offset_O;
	float3 dirtolight = normalize(distance_to_light);
	Ray occlusion_ray = Ray(offset_O, dirtolight, length(distance_to_light));
	//color += scene.ambient * m.albedo;
	if (ray.depthidx > max_depth) {
		// at maximum depth we try to return the last object hit's color, or just darkness for "eternal reflection"
		if (scene.IsOccluded(occlusion_ray))
		{
			return color;
		}
		if (d > 0.0) color += d * scene.directIllumination(ray.objIdx, I, N, m.albedo); // If diffuse
		return color;
	}
	if (!scene.IsOccluded(occlusion_ray)) {
		// -----------------------------------------------------------
		// less efficient
		// -----------------------------------------------------------
		//return  ray.dist.color * (d * directIllum + s * Trace(ray.reflect(I, N, ray.t + 1)));

		// -----------------------------------------------------------
		// more efficient
		// -----------------------------------------------------------

		if (d > 0.0) color += d * scene.directIllumination(ray.objIdx, I, N, m.albedo); // If diffuse
	}
	if (s > 0.0) {
		// Only if we hit front of the material
		color += s * Trace(ray.Reflect(I, N, ray.depthidx)); // If specular
	}
	return color;

	/* visualize normal */// return (N + 1) * 0.5f;//* directIllum;
	/*return*/  //col_;//(N + 1) * directIllum;
	/* visualize distance */   //return 0.1f * float3( ray.t, ray.t, ray.t );
	/* visualize albedo */  //return albedo ;
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
	// animation
	static float animTime = 0;
	scene.SetTime(animTime += deltaTime * 0.002f);

	//Camera
	const float speed = 0.2f;
	camera.move(mov, speed);
	camera.fovUpdate(fovc,speed);
	// pixel loop
	Timer t;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		if (y != 0) {
			int s = 0;
		}
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++) {
			accumulator[x + y * SCRWIDTH] =
				float4(Trace(camera.GetPrimaryRay(x, y)), 0); // Note for myself to remember: x+y*SCRWIDTH simply make sure we have a consecative num seq.

		}

		// translate accumulator contents to rgb32 pixels
		for (int dest = y * SCRWIDTH, x = 0; x < SCRWIDTH; x++)
			screen->pixels[dest + x] =
			RGBF32_to_RGB8(&accumulator[x + y * SCRWIDTH]);
	}
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
	//printf( "%5.2fms (%.1fps) - %.1fMrays/s\n", avg, fps, rps / 1000000 );
}

void Tmpl8::Renderer::MouseMove(int x, int y)
{
	int2 mouseDiff = mousePos - int2(x, y);
	const float angular_speed = 0.003f;
	camera.rotate(mouseDiff, angular_speed);
	mousePos = int2(x, y);
}

void Tmpl8::Renderer::KeyUp(int key)
{
	if (key == 0x57) mov -= float3(0, 0, 1);  // W
	if (key == 0x41) mov -= float3(-1, 0, 0);  //A
	if (key == 0x53) mov -= float3(0, 0, -1); //S
	if (key == 0x44) mov -= float3(1, 0, 0);; //D
	if (key == VK_SPACE) mov -= float3(0, 1, 0); // Space bar
	if (key == 0x43) mov -= float3(0, -1, 0); // C
	if (key == 0x45) fovc -= float3(-1,0,0); // E
	if (key == 0x51) fovc -= float3(1, 0, 0); // Q

}

void Tmpl8::Renderer::KeyDown(int key)
{
	if (key == 0x57) mov += float3(0, 0, 1);  // W
	if (key == 0x41) mov += float3(-1, 0, 0);  //A
	if (key == 0x53) mov += float3(0, 0, -1); //S
	if (key == 0x44) mov += float3(1, 0, 0);; //D
	if (key == VK_SPACE) mov += float3(0, 1, 0); // Space bar
	if (key == 0x43) mov += float3(0, -1, 0); // C
	if (key == 0x45) fovc += float3(-1, 0, 0);// E
	if (key == 0x51) fovc += float3(1, 0, 0);// Q


}
