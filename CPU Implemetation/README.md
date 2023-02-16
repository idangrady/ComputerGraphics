
# Project Description

<img src="https://user-images.githubusercontent.com/63195210/219384232-d9ae721d-e61c-48fc-a08e-df947853c9f0.png" width="400">


This is assignment 2 for Advance Computer Graphic course at Utrecht University. In this project, we have built upon the ray tracing algorithm that we developed for our first assignment. In particular, we have implemented a BVH traversal on the CPU to enhance the speed and efficiency of our ray tracing. Our BVH traversal algorithm allows us to rapidly traverse the scene and quickly identify which objects the rays intersect with, resulting in faster and more accurate rendering.

## Outputs
<div style="display:flex;">
  <img src="https://user-images.githubusercontent.com/63195210/219384688-d1bdc671-1adb-4d24-8730-b9618ca07dfe.png" width="300">
  <img src="https://user-images.githubusercontent.com/63195210/219384312-578b400e-7edc-47eb-a116-9565b818cc68.png" width="300">
  <img src="https://user-images.githubusercontent.com/63195210/219384651-1b7f7487-58c7-4a9f-a8a2-9e23bffd832a.png" width="300">
</div>

## Point For Improvement
We built the simple version of the BVH (article 1). Yet, as we experimented, we noticed a strange phenomenon. A scene with 18 triangles had only 3 BVH nodes, where the right and left children of the root nodes had the same AABB coordinates. A specific scene construction was the main cause of this phenomenon. Due to a lack of time, this was not fixed for this submission, but will be for the following. 
