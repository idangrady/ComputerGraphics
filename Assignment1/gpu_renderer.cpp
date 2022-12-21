#include "precomp.h"
#include <chrono>;
#pragma OPENCL EXTENSION cl_nvidia_printf : enable

const uint32_t RAY_SIZE = 64;

void GPURenderer::Init() {
	frame = new int[1];
	frame[0] = 0;
	seedDepth = new uint[2];
	seedDepth[0] = 0x1234578;
	seedDepth[1] = 0;
	justMoved = new bool;

	counters = new int[1];
	counters[0] = 0; // Index 0 for normal rays
	counterBuffer = new Buffer(sizeof(int) * 1, counters, 0);
	newRayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * RAY_SIZE, 0, 0);
	seedBuffer = new Buffer(sizeof(uint) * 2, seedDepth, CL_MEM_READ_WRITE);

	// We will explicity not use the (CPU located) screen as intermediate buffer, 
	// instead we will use GPU to draw directly to the renderTarget, skipping an
	// expensive copy from GPU to CPU.
	screen = 0;
	// Bind render target to OpenCL buffer
	screenBuffer = new Buffer(GetRenderTarget()->ID, 0, Buffer::TARGET);

	//frameCountBuffer = new Buffer(sizeof(uint), frame, CL_MEM_READ_ONLY);
	cameraBuffer = new Buffer(sizeof(float4) * 4, camera.cameraFloats, CL_MEM_READ_ONLY);
	rayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * RAY_SIZE, 0, 0);

	// Set the accumulator buffer. We write all the ray results to this 
	// buffer, and then run a Kernel that copies it to screenBuffer
	accumulatorBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), 0, 0);

	// Generate Kernel
	generateKernel = new Kernel("Kernels/generate.cl", "generate");
	// Extend Kernel
	extendKernel = new Kernel("Kernels/extend.cl", "extend");
	// Shade Kernel
	shadeKernel = new Kernel("Kernels/shade.cl", "shade");

	// Make some dummy triangles
	triBuffer = new Buffer(scene.tri_count * sizeof(TriGPU), scene.tris, CL_MEM_READ_ONLY);
	matBuffer = new Buffer(scene.tri_count * sizeof(MaterialGPU), scene.mats, CL_MEM_READ_ONLY);
	//triColorBuffer = new Buffer(2 * sizeof(cl_float4), scene.triColors, CL_MEM_READ_ONLY);

	// Generate Kernel arguments
	generateKernel->SetArguments(rayBuffer, cameraBuffer);
	// Extend Kernel Arguments
	extendKernel->SetArguments(rayBuffer, triBuffer, scene.tri_count);
	// Shade Kernel Arguments
	shadeKernel->SetArguments(rayBuffer, triBuffer, matBuffer, accumulatorBuffer, counterBuffer, newRayBuffer, seedBuffer);

	// Screen kernel
	screenKernel = new Kernel("Kernels/screen.cl", "renderToScreen");
	screenKernel->SetArguments(accumulatorBuffer, screenBuffer);

	// Clear kernel. Clears the accumulator.
	clearKernel = new Kernel("Kernels/clear.cl", "clear");
	clearKernel->SetArguments(accumulatorBuffer);
	// Copy Kernel for copying rays to ray buffer
	copyKernel = new Kernel("Kernels/enqueue.cl", "copy");
	copyKernel->SetArguments(rayBuffer, newRayBuffer);

	cameraBuffer->CopyToDevice();
	triBuffer->CopyToDevice();
	matBuffer->CopyToDevice();
	counterBuffer->CopyToDevice();
	seedBuffer->CopyToDevice();

	//triColorBuffer->CopyToDevice();
}

void Tmpl8::GPURenderer::Tick(float deltaTime)
{
	//int info = clGetCommandQueueInfo(Kernel::GetQueue(), CL_QUEUE_REFERENCE_COUNT, 4, &c_size, &ac_csize);
	//if (info == CL_SUCCESS) cout << "Size of queue:" << c_size << endl;
	//else cout << "Error getting command queue info.\n";
	if (length(mov) > 0) *justMoved = true;
	const float speed = 0.02f;
	camera.move(mult * mov, speed);
	cameraBuffer->CopyToDevice();
	Timer t;
	frame[0] += 1;
	cl_event k_events[3];
	// Clear accumulator
	clearKernel->Run(SCRWIDTH * SCRHEIGHT, 0, 0, &k_events[0]);
	// Run generate Kernel. Creates SCRWIDTH * SCRHEIGHT primary rays
	generateKernel->Run(SCRWIDTH * SCRHEIGHT, 0, &k_events[0], &k_events[1]);
	clFinish(Kernel::GetQueue());
	int val = SCRWIDTH * SCRHEIGHT;
	int it = 0;
	seedDepth[1] = 0;
	while (true) {
		seedDepth[0] = RandomUInt(seedDepth[0]);
		seedDepth[1] += 1;
		counters[0] = 0;
		seedBuffer->CopyToDevice();
		counterBuffer->CopyToDevice();
		cl_event c_events[3];
		// Run Extend Kernel
		extendKernel->Run(val, 0, &k_events[1], &c_events[0]);
		// Run Shade Kernel
		shadeKernel->Run(val, 0, &c_events[0], 0);
		clFinish(Kernel::GetQueue());
		counterBuffer->CopyFromDevice(true);
		if (counters[0] == 0) break;
		val = counters[0];
		copyKernel->Run(val, 0, 0, &k_events[1]);
	}
	// Draw to Screen
	screenKernel->Run(SCRWIDTH * SCRHEIGHT, 0, &k_events[1], 0);
	// Wait for finish
	clFinish(Kernel::GetQueue());
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
	printf("%5.2fms (%.1fps) - %.1fMrays/s\n", avg, fps, rps / 1000000);
	//*justMoved = false;
}

void Tmpl8::GPURenderer::MouseMove(int x, int y)
{
	int2 mouseDiff = int2(SCRWIDTH / 2, SCRHEIGHT / 2) - int2(x, y);
	const float angular_speed = 0.003f;
	camera.rotate(mouseDiff, angular_speed);
	*justMoved = true;
}

void Tmpl8::GPURenderer::KeyUp(int key)
{
	if (key == 0x58) mult = 1.0f;
	if (key == 0x57) mov -= float3(0, 0, 1);  // W
	if (key == 0x41) mov -= float3(-1, 0, 0);  //A
	if (key == 0x53) mov -= float3(0, 0, -1); //S
	if (key == 0x44) mov -= float3(1, 0, 0);; //D
	if (key == VK_SPACE) mov -= float3(0, 1, 0); // Space bar
	if (key == 0x43) mov -= float3(0, -1, 0); // C
	if (key == 0x45) fovc -= float3(-1, 0, 0); // E
	if (key == 0x51) fovc -= float3(1, 0, 0); // Q
}

void Tmpl8::GPURenderer::KeyDown(int key)
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