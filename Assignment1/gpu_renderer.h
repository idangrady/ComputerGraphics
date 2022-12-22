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
		int* framesSinceLastMoved;
		float mult = 1.0f;
		float3 mov;
		float3 fovc; // interactive FOV
		uint* seedDepth;

		SceneGPU scene;

		int* counters; // Counters for new rays
		static inline Buffer* counterBuffer;
		static inline Buffer* newRayBuffer;
		static inline Buffer* seedBuffer;
		static inline Buffer* movedBuffer;

		Camera camera;
		static inline Buffer* cameraBuffer;
		static inline Buffer* rayBuffer;

		static inline Kernel* generateKernel;
		static inline Kernel* extendKernel;
		static inline Kernel* shadeKernel;
		//static inline Kernel* connectKernel; // Don't need this I think

		static inline Buffer* triBuffer;
		static inline Buffer* matBuffer;
		static inline Buffer* triExBuffer;

		static inline Kernel* screenKernel;
		static inline Buffer* screenBuffer;

		static inline Buffer* accumulatorBuffer;
		static inline Buffer* intermediateBuffer;
		static inline Buffer* frameCountBuffer; //We really need a buffer for this?

		// Kernel that copies new rays to rays
		static inline Kernel* copyKernel;
		// Kernel that clears accumulator
		static inline Kernel* clearKernel;
		static inline Kernel* resetKernel;

		// Config options
		int max_depth = MAX_DEPTH;
		bool sendWhitted = sendWhittedCONFIG;
		uint8_t num_antiAlias = AA_COUNT;
	};
}