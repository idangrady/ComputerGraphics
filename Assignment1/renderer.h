#pragma once

namespace Tmpl8
{

	class Renderer : public TheApp
	{
	public:
		// game flow methods
		void Init();
		float3 Trace(Ray& ray);
		float3 Whitted(float3 I, float3 N, Ray& ray, Material& m);
		float3 RE(float3 I, float3 N, Ray& ray, Material& m);
		float4 Antialiasing(int  x, int y);
		void Tick(float deltaTime);
		void Shutdown() { /* implement if you want to do something on exit */ }

		// input handling
		void MouseUp(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseDown(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseMove(int x, int y);
		void MouseWheel(float y) { /* implement if you want to handle the mouse wheel */ }
		void KeyUp(int key);
		void KeyDown(int key);

		// data members
		int2 mousePos;
		float4* accumulator;
		int* accumulator_visit;

		bool sendWhitted = true;
		bool static_scene = false;

		uint8_t num_antiAlias = 5;

		float3 mov;
		float3 fovc; // interactive FOV
		Scene scene;
		Camera camera;

		int max_depth = 20;
	};

} // namespace Tmpl8
