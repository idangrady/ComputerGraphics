#include "precomp.h"
// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------

int i = 1;
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	//accumulator_visit = (int*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	frame = 0;

	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
	//memset(accumulator_visit, 0, SCRWIDTH * SCRHEIGHT * 16);

	POINT cursorPosition;
	GetCursorPos(&cursorPosition);
#if PRETTY
	camera.move(float3(180, 150, 180), 1.0f);
	mat4 m = mat4::LookAt(camera.camPos, float3(0, 128, 0));
	camera.matRotate(m);
#else
	camera.move(float3(-3, 0, -4), 1.0f);
	mat4 m = mat4::LookAt(camera.camPos, float3(0, -1, 0));
	camera.matRotate(m);

	//scene.bvh->BuildBVH();
#endif
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace(Ray& ray)
{
	scene.FindNearest(ray);
	uint typeId = GetObjectType(ray.I.instPrim);
	uint id = GetObjectIndex(ray.I.instPrim);
	//if (ray.objIdx == -1) return float3(0.5f, 0.5f, 0.5f); // or a fancy sky color
	//if (typeId == skyBoxID) return float3(0, 0, 0); // The void
	if (typeId == skyBoxID) return scene.getSkyBox(ray.D); // fancy sky color
	if (id == scene.areaID) return scene.area_lights[0].getcolor(); // I think this should be normlized other wise it will cause some issue 
	if (id == scene.spotID) return scene.area_lights[0].getcolor(); // I think this should be normlized other wise it will cause some issue 


	float3 I = ray.O + ray.I.t * ray.D;
	bool hit_back = false;
	float3 N = scene.GetNormal(ray.I, I, ray.D, hit_back);
	Material& m = scene.getMaterial(ray.I.instPrim);
	//if (m.isLight) return scene.GetColor(I);

	if (sendWhitted) { return Whitted(N, ray, m, hit_back); }
	else return RE(N, ray, m, hit_back);

}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
	// animation
	static float animTime = 0;
#if !STATIC
	scene.SetTime(animTime += deltaTime * 0.002f);
#endif
	//Camera
	const float speed = 0.2f;
	camera.move(mult * mov, speed);
	camera.fovUpdate(fovc,speed);
	// pixel loop
	Timer t;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
	frame += 1.f;
#if STATIC
	cout << "Frame: " << frame << endl;
#endif
	float invFrame = 1.f / frame;
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		if (y != 0) {
			int s = 0;
		}
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++) {
#if STATIC
			if (num_antiAlias > 1) {																		// anti aliasing over num_antiAlias values
				accumulator[x + y * SCRWIDTH] += Antialiasing(x, y);
			}
			accumulator[x + y * SCRWIDTH] += float4(Trace(camera.GetPrimaryRay(x, y)), 0);
#else
			if (num_antiAlias > 1) {																		// anti aliasing over num_antiAlias values
				accumulator[x + y * SCRWIDTH] = Antialiasing(x, y);
			}
			else accumulator[x + y * SCRWIDTH] = float4(Trace(camera.GetPrimaryRay(x, y)), 0);
#endif
		}
		// translate accumulator contents to rgb32 pixels
		for (int dest = y * SCRWIDTH, x = 0; x < SCRWIDTH; x++)
#if STATIC
		{
			float4 color = accumulator[x + y * SCRWIDTH];
			color *= invFrame;
			screen->pixels[dest + x] = RGBF32_to_RGB8(&color);
		}
#else
			screen->pixels[dest + x] = RGBF32_to_RGB8(&accumulator[x + y * SCRWIDTH]);
#endif
		}
	
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
}

#if !STATIC
void Tmpl8::Renderer::MouseMove(int x, int y)
{
	int2 mouseDiff = int2(SCRWIDTH / 2, SCRHEIGHT / 2) - int2(x, y);
	const float angular_speed = 0.003f;
	camera.rotate(mouseDiff, angular_speed);
}

void Tmpl8::Renderer::KeyUp(int key)
{
	if (key == 0x58) mult = 1.0f;
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
	if (key == 0x58) mult = 100.0f;
	if (key == 0x57) mov += float3(0, 0, 1);  // W
	if (key == 0x41) mov += float3(-1, 0, 0);  //A
	if (key == 0x53) mov += float3(0, 0, -1); //S
	if (key == 0x44) mov += float3(1, 0, 0);; //D
	if (key == VK_SPACE) mov += float3(0, 1, 0); // Space bar
	if (key == 0x43) mov += float3(0, -1, 0); // C
	if (key == 0x45) fovc += float3(-1, 0, 0);// E
	if (key == 0x51) fovc += float3(1, 0, 0);// Q
}
#else
void Tmpl8::Renderer::MouseMove(int x, int y){}
void Tmpl8::Renderer::KeyUp(int key) {}
void Tmpl8::Renderer::KeyDown(int key) {}
#endif

float Tmpl8::Renderer::FresnelReflection(float cos_theta_i, float n1, float n2, float refr, float3 N, float3 D)
{
	float sin_theta_i = length(cross(N, -D));
	float refr_sin = (refr * sin_theta_i);
	float cos_theta_t = sqrtf(1.0f - (refr_sin * refr_sin));
	float first_term = ((n1 * cos_theta_i) - (n2 * cos_theta_t)) / ((n1 * cos_theta_i) + (n2 * cos_theta_t));
	float second_term = ((n1 * cos_theta_t) - (n2 * cos_theta_i)) / ((n1 * cos_theta_t) + (n2 * cos_theta_i));
	return 0.5f * ((first_term * first_term) + (second_term * second_term));
}

// Use Schlick to calculate if it fully reflects
float Tmpl8::Renderer::GetSnell(float refr, float cos_theta_i)
{
	float k = 1.0f - (refr * refr) * (1.0f - (cos_theta_i * cos_theta_i));
	return k;
}

float3 Tmpl8::Renderer::Whitted(float3 N, Ray& ray, Material& m, bool hit_back)
{
	float3 color = (0, 0, 0);
	float3 I_loc = ray.O + ray.I.t * ray.D;

	float s = m.specularity;
	float d = 1.0f - s;

	//Ray secondary_ray = ray.reflect(I, N, ray.depthidx);
	float3 dirtolight = scene.spot_lights[0].position - I_loc;
	float distance_to_light = length(dirtolight);
	dirtolight /= distance_to_light;
	float3 offset_O = I_loc + (0.0002f * dirtolight);
	Ray occlusion_ray = Ray(offset_O, dirtolight, distance_to_light - 0.0002f);
	//color += scene.ambient * m.albedo;
	if (m.mat_medium == Medium::Glass) {
		float refr, n1, n2;
		if (hit_back) { // From Glass to Air
			refr = scene.refractiveTransmissions[Medium::Glass][Medium::Air];
			n1 = scene.refractiveIndex[Medium::Glass];
			n2 = scene.refractiveIndex[Medium::Air];
		}
		else { // From Air to Glass
			refr = scene.refractiveTransmissions[Medium::Air][Medium::Glass];
			n1 = scene.refractiveIndex[Medium::Air];
			n2 = scene.refractiveIndex[Medium::Glass];
		}
		float R;
		float T;
		float traveled = ray.I.t;
		float3 interim_color = float3(0, 0, 0);
		float cos_theta_i = dot(N, -ray.D);
		float k = GetSnell(refr, cos_theta_i);
		if (k > 0.00001f) { // Inner scope cuz lots of terms. Also we use k > 0.0001f to account for stupid float inaccuracies.
			R = FresnelReflection(cos_theta_i, n1, n2, refr, N, ray.D);
			T = 1.0f - R;
			// Double check for correctness later
			float3 t_dir = (refr * ray.D) + (N * (refr * cos_theta_i - sqrtf(k)));
			t_dir /= length(t_dir);
			interim_color += T * Trace(Ray(I_loc + (0.0002f * t_dir), t_dir, 1e34f, ray.depthidx + 1));
		}
		else
		{
			R = 1.0f;
			T = 0.0f;
		}
		if (R > 0.0f && ray.depthidx <= max_depth) {
			interim_color += R * Trace(ray.Reflect(I_loc, N));
		}

		if (hit_back) { // If we go from glass to air, we have to absorb some of the light we found (because we traverse in reverse order!)
			interim_color.x *= exp(-m.absorption.x * traveled);
			interim_color.y *= exp(-m.absorption.y * traveled);
			interim_color.z *= exp(-m.absorption.z * traveled);
		}
		color += interim_color;
	}
	else {
		if (s > 0.0f && ray.depthidx <= max_depth)
		{
			color += s * Trace(ray.Reflect(I_loc, N)); // If specular
		}
		if (!scene.IsOccluded(occlusion_ray))
		{
			if (d > 0.0f) 
			{
				color += d * scene.directIllumination(ray.I, I_loc, N); // If diffuse
			}
		}
	}
	return color;
}


float3 Tmpl8::Renderer::RE(float3 N, Ray& ray, Material& m, bool hit_back)
{	
	float3 I_loc = ray.O + ray.I.t * ray.D;
	float s = m.specularity;
	float d = 1.0f - s;
	// path tracing 

	if (m.mat_medium == Medium::Glass) {
		float refr, n1, n2;
		if (hit_back) { // From Glass to Air
			refr = scene.refractiveTransmissions[Medium::Glass][Medium::Air];
			n1 = scene.refractiveIndex[Medium::Glass];
			n2 = scene.refractiveIndex[Medium::Air];
		}
		else { // From Air to Glass
			refr = scene.refractiveTransmissions[Medium::Air][Medium::Glass];
			n1 = scene.refractiveIndex[Medium::Air];
			n2 = scene.refractiveIndex[Medium::Glass];
		}
		float R;
		float traveled = ray.I.t;
		float3 interim_color = float3(0, 0, 0);
		float cos_theta_i = dot(N, -ray.D);
		float k = GetSnell(refr, cos_theta_i);
		if (k > 0.00001f) {
			R = FresnelReflection(cos_theta_i, n1, n2, refr, N, ray.D);
			float random = RandomFloat();
			if (random > R) { //Randomly refracted
				float3 t_dir = (refr * ray.D) + (N * (refr * cos_theta_i - sqrtf(k)));
				t_dir /= length(t_dir);
				interim_color = Trace(Ray(I_loc + (0.0002f * t_dir), t_dir, 1e34f, ray.depthidx + 1));
			}
			else if(ray.depthidx <= max_depth) { // Randomly reflected
				interim_color = Trace(ray.Reflect(I_loc, N));
			}
		}
		else if(ray.depthidx <= max_depth) {
			interim_color = Trace(ray.Reflect(I_loc, N));
		}
		if (hit_back) { // If we go from glass to air, or stay inside glass, we have to absorb some of the light we found (because we traverse in reverse order!)
			interim_color.x *= exp(-m.absorption.x * traveled);
			interim_color.y *= exp(-m.absorption.y * traveled);
			interim_color.z *= exp(-m.absorption.z * traveled);
		}
		return interim_color;
		throw exception("Reached unreachable part?");
	}
	else {
		float random = RandomFloat();
		if (random > d) // Randomly reflect
		{ // Mirror
			return Trace(ray.Reflect(I_loc, N));
		}
		else if(ray.depthidx <= max_depth){ // Randomly diffuse
			float3 BRDF_m = scene.GetColor(ray.I, I_loc) * PI; // I deleted the PI because it was cancalled in the return * PI 
			float3 random_dir = scene.GetDiffuseRefelectDir(N);
			Ray newRay(I_loc + 0.0002f * random_dir, random_dir, 1e34f, ray.depthidx + 1); //+ 0.0002f * random_dir
			float3 EI = Trace(newRay) * dot(N, random_dir);
			float3 r_val = BRDF_m * EI * INVPI;
			if (length(EI) < length(r_val)) cout << "MORE OUTGOING THAN INCOMING??" << endl;
			return r_val;
		}
		else {
			return float3(0, 0, 0);
		}
	}
	throw exception("Reached unreachable part 2?");
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
