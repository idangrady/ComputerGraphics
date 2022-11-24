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
		camPos = float3( 0, 0, -2 );
		topLeft = float3( -aspect, 1, 0 );
		topRight = float3( aspect, 1, 0 );
		bottomLeft = float3( -aspect, -1, 0 );

		auto dist_to_focus = length(camPos - ((topRight- bottomLeft) / 2)); // TODO: This should be the depth.

	}

	void move(float speed, float3 dircm_dir) {
		//float3 dircm_dir = getdirection()* side;
		camPos += speed * dircm_dir ;
		topRight += speed * dircm_dir ;
		bottomLeft += speed * dircm_dir ;
		topLeft += speed * dircm_dir ;
	}
	float3 getdirection() {
		return normalize((topRight - bottomLeft)/2 - camPos );
	}

	void rotate_camAxis(mat4 TMatrix)
	{

			auto ss = topLeft - camPos;
			topLeft =  ((float3)(TMatrix * (topLeft - camPos)) + camPos);
			topRight =  ((float3)(TMatrix * (topRight - camPos)) + camPos);
			bottomLeft = ((float3)(TMatrix * (bottomLeft - camPos)) + camPos);
		
	}

	void rotate_cam(float3 rot)
	{
		if (rot.x + rot.y + rot.z != 0)
		{
			mat4 TMatrix = mat4::Rotate(normalize(rot), length(rot));
			auto ss = topLeft - camPos;
			topLeft = (float3)(TMatrix * (topLeft - camPos)) + camPos;
			topRight = (float3)(TMatrix * (topRight - camPos)) + camPos;
			bottomLeft = (float3)(TMatrix * (bottomLeft - camPos)) + camPos;
		}
	}

	inline double degrees_to_radians(double degrees) {
		return degrees * PI / 180.0;
	}

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
	float3 camPos;
	float3 topLeft, topRight, bottomLeft;
	


	//mat4 mat = mat4();

	//auto theta = degrees_to_radians(vfov);

	//auto h = tan(theta / 2);
	//auto viewport_height = 2.0 * h;
	//auto viewport_width = aspect * viewport_height;

};

}