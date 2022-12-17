#pragma once
#include<config.h>
namespace Tmpl8
{
	class GPURenderer : public TheApp {
		void Init();
		void Tick(float deltaTime);
		void Shutdown() { /* implement if you want to do something on exit */ }

		// input handling
		void MouseUp(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseDown(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseWheel(float y) { /* implement if you want to handle the mouse wheel */ }
		void MouseMove(int x, int y);
		void KeyUp(int key);
		void KeyDown(int key);

		static inline int* frame;
		float mult = 1.0f;
		float3 mov;
		float3 fovc; // interactive FOV

		Scene scene;
		Camera camera;

		static inline Kernel* screenKernel;

		static inline Buffer* accumulatorBuffer;
		static inline Buffer* screenBuffer;
		static inline Buffer* frameCountBuffer; //We really need a buffer for this?

		// Config options
		int max_depth = MAX_DEPTH;
		bool sendWhitted = sendWhittedCONFIG;
		uint8_t num_antiAlias = AA_COUNT;
	};
}