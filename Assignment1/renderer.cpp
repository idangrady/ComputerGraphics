#include "precomp.h"
#include "userInput.h"
// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64( SCRWIDTH * SCRHEIGHT * 16 );
	memset( accumulator, 0, SCRWIDTH * SCRHEIGHT * 16 );
	

	mousePos = User.get_mouse();
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace( Ray& ray )
{
	scene.FindNearest( ray );
	if (ray.objIdx == -1) return 0; // or a fancy sky color
	float3 I = ray.O + ray.t * ray.D;
	float3 N = scene.GetNormal( ray.objIdx, I, ray.D );
	float3 albedo = scene.GetAlbedo( ray.objIdx, I );
	/* visualize normal */ return (N + 1) * 0.5f;
	/* visualize distance */ // return 0.1f * float3( ray.t, ray.t, ray.t );
	/* visualize albedo */ // return albedo;
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick( float deltaTime )
{

	updateKey();
	// animation
	static float animTime = 0;
	scene.SetTime( animTime += deltaTime * 0.002f );
	// pixel loop
	Timer t;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
			accumulator[x + y * SCRWIDTH] =
				float4( Trace( camera.GetPrimaryRay( x, y ) ), 0 );
		// translate accumulator contents to rgb32 pixels
		for (int dest = y * SCRWIDTH, x = 0; x < SCRWIDTH; x++)
			screen->pixels[dest + x] = 
				RGBF32_to_RGB8( &accumulator[x + y * SCRWIDTH] );
	}
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
	printf( "%5.2fms (%.1fps) - %.1fMrays/s\n", avg, fps, rps / 1000000 );
}


void Renderer::updateKey() {
	float3 norm = normalize(camera.topLeft - camera.topRight);
	int action = User.get_input();


	float angular = PI / 2 * 0.003f;
	float speed = 0.1f;
	float3 direction = camera.getdirection();
	float3 forward = direction * float3(0, 0, 1);
	float3 side = normalize(camera.topRight - camera.topLeft);


	int2 curMousePose = User.get_mouse();

	if (curMousePose != mousePos) {
		camera.rotate_cam(normalize(float3(0, (mousePos-curMousePose ).x, 0)) * speed*0.3);
		MouseMove(curMousePose.x, curMousePose.y);
	}


	switch (action)
	{
	case 1:
		camera.move(speed, side);
		break;
	case 2:
		camera.move(-speed, side);
		break;
	case 3:
		camera.move(speed, forward);
		break;
	case 4:
		camera.move(-speed, forward);
		break;
	case 5:
		camera.rotate_cam(normalize(float3(0, -angular, 0))*speed);
		break;
	case 6:
		camera.rotate_cam(normalize(float3(0, angular, 0))*speed);
		break;
	case 7:
		camera.rotate_cam(normalize(float3(1, 0, 0))*speed);
		break;
	case 8:
		camera.rotate_cam(normalize(float3(1, 0, 0))*-speed);
		break;
	}

};
