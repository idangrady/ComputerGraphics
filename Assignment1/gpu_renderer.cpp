#include "precomp.h"

void GPURenderer::Init() {
	frame = new int[1];
	frame[0] = 0;

	// We will explicity not use the (CPU located) screen as intermediate buffer, 
	// instead we will use GPU to draw directly to the renderTarget, skipping an
	// expensive copy from GPU to CPU.
	screen = 0;
	// Bind render target to OpenCL buffer
	screenBuffer = new Buffer(GetRenderTarget()->ID, 0, Buffer::TARGET);

	frameCountBuffer = new Buffer(sizeof(uint), frame, CL_MEM_READ_ONLY);

	// Screen kernel
	screenKernel = new Kernel("Kernels/example.cl", "renderToScreen");
	screenKernel->SetArguments(screenBuffer, frameCountBuffer);


	MakePrimaryRays = new Kernel("Kernels/MakePrimaryRays.cl", "renderToScreen");
	primaryRays = new Buffer(GetRenderTarget()->ID, 0, Buffer::TARGET);
	MakePrimaryRays->SetArguments(MakePrimaryRays, &camera.camPos);




	// Set the accumulator buffer. Currently unused, but we should probably write all the ray results 
	// to this buffer, and then run a Kernel that copies it to screenBuffer
	accumulatorBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * sizeof(float4), 0, 0);
}

void Tmpl8::GPURenderer::Tick(float deltaTime)
{
	Timer t;
	frame[0] += 1;
	// Copy frameCountBuffer to GPU
	frameCountBuffer->CopyToDevice(true);

	screenKernel->Run(SCRWIDTH * SCRHEIGHT);
	// Wait for finish
	clFinish(Kernel::GetQueue());
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
	printf("%5.2fms (%.1fps) - %.1fMrays/s\n", avg, fps, rps / 1000000);
}

void Tmpl8::GPURenderer::MouseMove(int x, int y)
{
	int2 mouseDiff = int2(SCRWIDTH / 2, SCRHEIGHT / 2) - int2(x, y);
	const float angular_speed = 0.003f;
	camera.rotate(mouseDiff, angular_speed);
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