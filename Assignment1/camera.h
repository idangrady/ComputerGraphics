#pragma once

// default screen resolution
#define SCRWIDTH	1280
#define SCRHEIGHT	720
// #define FULLSCREEN
// #define DOUBLESIZE

namespace Tmpl8 {

class Camera
{
public:
	Camera()
	{

		double vfov = 90;
		float theta = degrees_to_radians(vfov);
		float h = tan(theta / 2.0f);
		viewWidth = aspect * h;



		// setup a basic view frustum
		camPos_start = float3(0, 80, -128);
		startDir = float3(0, 0, 1);
		camPos = camPos_start;
		topLeft_start = float3( -viewWidth, camPos_start.y + 1, camPos_start.z + 2 );
		topRight_start = float3(viewWidth, camPos_start.y + 1, camPos_start.z + 2); // changed
		bottomLeft_start = float3( -viewWidth, -camPos_start.y - 1, camPos_start.z + 2);
		topLeft = topLeft_start;
		topRight = topRight_start;
		bottomLeft = bottomLeft_start;

		auto dist_to_focus = length(camPos - ((topRight- bottomLeft) / 2)); // TODO: This should be the depth.

	}

	void rotate(int2 mouse_movement, float speed) {
		yaw += mouse_movement.x * speed;
		pitch += mouse_movement.y * speed;
		while (yaw > PI) yaw -= 2 *PI;
		while (yaw < -PI) yaw += 2 * PI;
		const float min_pitch = -(PI / 2) + FLT_MIN;
		const float max_pitch = (PI / 2) - FLT_MIN;
		pitch = clamp(pitch, min_pitch, max_pitch);
		topLeft =  mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * float3(-viewWidth, 1, 2));
		topRight = mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * float3(viewWidth, 1, 2));
		bottomLeft = mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * float3(-viewWidth, -1, 2));
		topRight += camPos;
		bottomLeft += camPos;
		topLeft += camPos;
	}

	void move(float3 dir, float speed) {
		// Movement in camera space
		float3 move = mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * dir);
		camPos += move * speed;
		topLeft += move * speed;
		topRight += move * speed;
		bottomLeft += move * speed;
	}
	void fovUpdate(float3 dir,  float speed) { // interactive FOV
		topLeft -= dir * speed;
		topRight += dir * speed;
		bottomLeft -= dir * speed;
	};

	float3 getdirection() {
		return normalize((topRight + bottomLeft)/2 - camPos );}


	Ray GetPrimaryRay( const int x, const int y )
	{
		// calculate pixel position on virtual screen plane
		const float u = (float)x * (1.0f / SCRWIDTH);
		const float v = (float)y * (1.0f / SCRHEIGHT);
		const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
		return Ray( camPos, normalize( P - camPos ) );
	}
	Ray GetRandomPrimaryRay(const int x, const int y) 
	{
		float x_ = (float)x + (2 * RandomFloat() - 1) + 0.5;					// + 0.5 to bring it to the middle of the pixel of the x axis
		float y_ = (float)y + (2 * RandomFloat() - 1) + 0.5;					// + 0.5 to bring it to the middle of the pixel of the y axis
	}

	Ray GetPrimaryRayRandomized(const float x, const float  y)								// I am still wondering whether to leave it like
																							// that (twoo similar func) , or to try to write it differently
	{
		// calculate pixel position on virtual screen plane for randomized 
		const float u = x * (1.0f / SCRWIDTH);
		const float v = y * (1.0f / SCRHEIGHT);
		const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
		return Ray(camPos, normalize(P - camPos));
	}


	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float ZoomLevel = 1.0f;
	float3 camPos, camPos_start, startDir;
	float yaw = 0, pitch = 0;
	mat4 cam_matrix;
	float3 topLeft_start, topRight_start, bottomLeft_start;
	float3 topLeft, topRight, bottomLeft;

	
	float vfov;
	float theta;
	float h;
	float viewWidth;
	float viewport_height;
	float viewport_width;

};

}