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
		// setup a basic view frustum
		camPos_start = float3(0, 0, -2);
		startDir = float3(0, 0, 1);
		camPos = camPos_start;
		topLeft_start = float3( -aspect, 1, 0 );
		topRight_start = float3( aspect, 1, 0 );
		bottomLeft_start = float3( -aspect, -1, 0 );
		topLeft = topLeft_start;
		topRight = topRight_start;
		bottomLeft = bottomLeft_start;

		auto dist_to_focus = length(camPos - ((topRight- bottomLeft) / 2)); // TODO: This should be the depth.

	}

	void rotate(int2 mouse_movement, float speed) {
		yaw += mouse_movement.x * speed;
		pitch += mouse_movement.y * speed;
		while (yaw > PI) yaw -= PI;
		while (yaw < -PI) yaw += PI;
		const float min_pitch = -(PI / 2) + FLT_MIN;
		const float max_pitch = (PI / 2) - FLT_MIN;
		pitch = clamp(pitch, min_pitch, max_pitch);
		topLeft =  mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * (topLeft_start - camPos_start));
		topRight = mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * (topRight_start - camPos_start));
		bottomLeft = mat4::RotateY(-yaw) * (mat4::RotateX(-pitch) * (bottomLeft_start - camPos_start));
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

	float3 getdirection() {
		return normalize((topRight + bottomLeft)/2 - camPos );
	}

	//void zoom(float speed) {
	//	float3 cam_dir = getdirection();
	//	topLeft += speed * cam_dir; topRight += speed * cam_dir; bottomLeft += speed * cam_dir;
	//}


	Ray GetPrimaryRay( const int x, const int y )
	{
		// calculate pixel position on virtual screen plane
		const float u = (float)x * (1.0f / SCRWIDTH);
		const float v = (float)y * (1.0f / SCRHEIGHT);
		const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
		return Ray( camPos, normalize( P - camPos ) );
	}
	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float ZoomLevel = 1.0f;
	float3 camPos, camPos_start, startDir;
	float yaw = 0, pitch = 0;
	mat4 cam_matrix;
	float3 topLeft_start, topRight_start, bottomLeft_start;
	float3 topLeft, topRight, bottomLeft;

	


	//mat4 mat = mat4();

	//auto theta = degrees_to_radians(vfov);

	//auto h = tan(theta / 2);
	//auto viewport_height = 2.0 * h;
	//auto viewport_width = aspect * viewport_height;

};

}