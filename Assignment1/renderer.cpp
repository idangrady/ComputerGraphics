#include "precomp.h"
// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------

int i = 1;
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	accumulator_visit = (int*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);

	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
	memset(accumulator_visit, 0, SCRWIDTH * SCRHEIGHT * 16);

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
	if (ray.objIdx == -1) return float3(0.5f, 0.5f, 0.5f); // or a fancy sky color

	float3 I = ray.O + ray.t * ray.D;
	float3 N = scene.GetNormal(ray.objIdx, I, ray.D);
	Material& m = scene.getMaterial(ray.objIdx);

	if (sendWhitted) { return Whitted(I, N, ray, m); }
	else return RE(I, N, ray, m);

	// earlier work
	// 
	//return color;
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
	if(!static_scene) scene.SetTime(animTime += deltaTime * 0.002f);

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
			if (static_scene) {
				accumulator_visit[x + y * SCRWIDTH] += 1;
				accumulator[x + y * SCRWIDTH] +=							
					float4(Trace(camera.GetPrimaryRay(x, y)), 0)/ accumulator_visit[x + y * SCRWIDTH];
			}
			else {
				if (num_antiAlias > 1) {																		// anti aliasing over num_antiAlias values
					accumulator[x + y * SCRWIDTH] = Antialiasing(x,y);}
				else accumulator[x + y * SCRWIDTH] = float4(Trace(camera.GetPrimaryRay(x, y)), 0);
				}
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

float3 Tmpl8::Renderer::Whitted(float3 I, float3 N, Ray& ray, Material& m)
{
	float d = m.diffuse;
	float s = m.specular;
	float3 color = (0, 0, 0);

	float3 distance_to_light = scene.lights[0].position - I;
	float3 dirtolight = normalize(distance_to_light);
	Ray occlusion_ray = Ray(I + 0.0002f * dirtolight, dirtolight, length(distance_to_light), ray.objIdx + 1);

	if (ray.depthidx <= max_depth) {
		if (s > 0.0)  color += s * Trace(ray.Reflect(I, N)); // If specular
		if (d > 0.0 && !scene.IsOccluded(occlusion_ray)) color += d * scene.directIllumination(ray.objIdx, I, N, m.albedo);	// If diffuse
	}
	else {
		// at maximum depth we try to return the last object hit's color, or just darkness for "eternal reflection"
		if (d > 0.0 && !scene.IsOccluded(occlusion_ray)) color += (d + s) * scene.directIllumination(ray.objIdx, I, N, m.albedo);
	}
	return color;
}


float3 Tmpl8::Renderer::RE(float3 I, float3 N, Ray& ray, Material& m)
{	
	float d = m.diffuse;
	float s = m.specular;
	// path tracing 
	if (ray.objIdx == 2) {
		//hit light
		return m.albedo; // I think this should be normlized other wise it will cause some issue 
	}
	if (ray.depthidx > max_depth) {
		return float3(0.9f, 0.9f, 0.9f);
	}
	else {
		float random = RandomFloat();
		if (random > d) // Randomly reflect
		{ // Mirror
			return Trace(ray.Reflect(I, N));
		}
		else { // Randomly diffuse
			float3 BRDF_m = m.albedo; // I deleted the PI because it was cancalled in the return * PI 
			float3 random_dir = scene.GetDiffuseRefelectDir(N);
			Ray newRay(I + 0.0002f * random_dir, random_dir, 1e34f, ray.depthidx + 1); //+ 0.0002f * random_dir
			float3 EI = Trace(newRay) * dot(N, random_dir);
			return 2.0f * BRDF_m * EI;
		}
	}
}

float4 Tmpl8::Renderer::Antialiasing(int x, int y)
{
	float4 color(0.0f, 0.0f, 0.0f, 0.0f);
	for (uint8_t i = num_antiAlias; i > 0; i--) {
		float x_ = (float)x + (2 * RandomFloat() - 1) +0.5;						// + 0.5 to bring it to the middle of the pixel of the x axis
		float y_ = (float)y + (2 * RandomFloat() - 1) + 0.5;					// + 0.5 to bring it to the middle of the pixel of the y axis
		color += float4(Trace(camera.GetPrimaryRayRandomized(x_, y_)), 0);		// I checked, if I send float x, float y to GetPrimaryRay, it convert it to int. 
	}
	return color / (float)num_antiAlias;

}
