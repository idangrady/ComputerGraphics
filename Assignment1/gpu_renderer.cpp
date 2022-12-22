#include "precomp.h"

#pragma OPENCL EXTENSION cl_nvidia_printf : enable

void GPURenderer::Init() {
	frame = new int[1];
	frame[0] = 0;

	// We will explicity not use the (CPU located) screen as intermediate buffer, 
	// instead we will use GPU to draw directly to the renderTarget, skipping an
	// expensive copy from GPU to CPU.
	screen = 0;
	// Bind render target to OpenCL buffer
	screenBuffer = new Buffer(GetRenderTarget()->ID, 0, Buffer::TARGET);

	//frameCountBuffer = new Buffer(sizeof(uint), frame, CL_MEM_READ_ONLY);
	cameraBuffer = new Buffer(sizeof(float4) * 4, camera.cameraFloats, CL_MEM_READ_ONLY);
	rayBuffer = new Buffer(SCRWIDTH * SCRHEIGHT * 64, 0, 0);

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
	triBuffer = new Buffer(3 * sizeof(Primitive_GPU), scene.arrPrimitive, CL_MEM_READ_ONLY);
	triColorBuffer = new Buffer(3 * sizeof(cl_float4), scene.triColors, CL_MEM_READ_ONLY);

	// Make some dummy triangles
	bvhNodeBuffer = new Buffer(3 * sizeof(BVHNode), scene.bvhNode, CL_MEM_READ_ONLY);					// the 3 has to change!
	arrPrimitiveBuffer = new Buffer(3 * sizeof(Primitive_GPU), scene.arrPrimitive, CL_MEM_READ_ONLY);		// the 3 has to change!
	arrPrimitiveIdxBuffer = new Buffer(3 * sizeof(cl_uint), scene.arrPrimitiveIdx, CL_MEM_READ_ONLY);		// the 3 has to change!


	// Generate Kernel arguments
	generateKernel->SetArguments(rayBuffer, cameraBuffer);
	// Extend Kernel Arguments
	extendKernel->SetArguments(rayBuffer,  triBuffer, bvhNodeBuffer, arrPrimitiveIdxBuffer, 3); //, bvhNodeBuffer,arrPrimitiveBuffer, arrPrimitiveIdxBuffer
	// Shade Kernel Arguments
	shadeKernel->SetArguments(rayBuffer, triBuffer, triColorBuffer, accumulatorBuffer);

	// Screen kernel
	screenKernel = new Kernel("Kernels/screen.cl", "renderToScreen");
	screenKernel->SetArguments(accumulatorBuffer, screenBuffer);

	// Clear kernel. Clears the accumulator.
	clearKernel = new Kernel("Kernels/clear.cl", "clear");
	clearKernel->SetArguments(accumulatorBuffer);

	cameraBuffer->CopyToDevice();
	triBuffer->CopyToDevice();
	triColorBuffer->CopyToDevice();
}

void Tmpl8::GPURenderer::Tick(float deltaTime)
{
	const float speed = 0.02f;
	camera.move(mult * mov, speed);
	cameraBuffer->CopyToDevice();
	Timer t;
	frame[0] += 1;
	cl_event k_events[5];
	// Clear accumulator
	clearKernel->Run(SCRWIDTH * SCRHEIGHT, 0, 0, &k_events[0]);
	// Run generate Kernel. Creates SCRWIDTH * SCRHEIGHT primary rays
	generateKernel->Run(SCRWIDTH * SCRHEIGHT, 0, &k_events[0], &k_events[1]);
	// Run Extend Kernel
	extendKernel->Run(SCRWIDTH * SCRHEIGHT, 0, &k_events[1], &k_events[2]);
	// Run Shade Kernel
	shadeKernel->Run(SCRWIDTH * SCRHEIGHT, 0, &k_events[2], &k_events[3]);
	// Draw to Screen
	screenKernel->Run(SCRWIDTH * SCRHEIGHT, 0, &k_events[3], &k_events[4]);
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