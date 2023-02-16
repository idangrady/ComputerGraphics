#include "Kernels/enqueue.cl"
#include "Kernels/extend.cl"
#include "Kernels/screen.cl"
#include "Kernels/shade.cl"
#include "Kernels/utils.cl"

__kernel void
schedule(__global Ray *rays, __constant Triangle *triangles,
         __constant TriEx *triExes, __constant Material *materials,
         __global float4 *intermediate, volatile __global int *counter,
         __global Ray *newRays, __global ulong *seeds, __global uint *depth,
         __constant float4 *skybox, __private int width, __private int height,
         __private int triangleCount, __global float4 *accumulator,
         write_only image2d_t target, __constant int *framesSinceMoved) {
  int val = atomic_xchg(counter, 0);
  if (val == 0) {
    ndrange_t n = ndrange_1D(SCRWIDTH * SCRHEIGHT);
    int suc =
        enqueue_kernel(get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, n, ^{
          renderToScreen(intermediate, accumulator, target, framesSinceMoved);
        });
    // printf("Succes: %i\n", suc);
  } else {
    ndrange_t n = ndrange_1D(val);
    ndrange_t r = ndrange_1D(1);
    int old = atomic_add(depth, 1);
    // depth[0] += 1;
    clk_event_t extendDone;
    clk_event_t shadeDone;
    clk_event_t copyDone;
    int suc;
    suc = enqueue_kernel(get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, n, 0,
                         NULL, &extendDone, ^{
                           extend(rays, triangles, triangleCount);
                         });
    release_event(extendDone);
    // printf("Succes: %i\n", suc);
    suc = enqueue_kernel(get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, n, 1,
                         &extendDone, &shadeDone, ^{
                           shade(rays, triangles, triExes, materials,
                                 intermediate, counter, newRays, seeds, depth,
                                 skybox, width, height);
                         });
    // printf("Succes: %i\n", suc);
    suc = enqueue_kernel(get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, n, 1,
                         &shadeDone, &copyDone, ^{
                           copy(rays, newRays);
                         });
    release_event(shadeDone);
    // printf("Succes: %i\n", suc);
    suc = enqueue_kernel(
        get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, r, 1, &copyDone, 0, ^{
          schedule(rays, triangles, triExes, materials, intermediate, counter,
                   newRays, seeds, depth, skybox, width, height, triangleCount,
                   accumulator, target, framesSinceMoved);
        });
    release_event(copyDone);
    // printf("Succes: %i\n", suc);
    // printf("iteration %i\n", old + 1);
  }
}