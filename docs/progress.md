# Progress report 6/5
We are trying to implement progressive photon mapping to better render dielectrics and caustic effects.
## stuff done so far
 - fix bugs in path tracing and multiple importance sampling.
 - refactor code for readability. Added more files to split up code and also split rendering by strategies.
 - implement MIS dielectrics
 - create a few scenes to help demonstrate caustics, render them
 - reading and start figuring out photon mapping algorithm and structures needed
 - Setup testing machine

## initial results
Here are some images we have rendered using the base renderer from the homework


The water with large light renders pretty good with MIS already. However when the area light is small it is very noisy and the image is not very clear even with 10000 samples per pixel.

## What we need to still do
 - integrate nanoflann(kd tree)
   - Transforming Vector3 and other Torrey defined types to fit into nanoflann
   - Working with the parallel support provided by both libraries
 - start implementation with Diffuse, Mirror, Dielectric materials first.
 - lighting, start with point lights, then go on to area lights.
 - Implement the phases of the algorithm
   1. tracing rays and getting "visible points"
   2. Disperse photons into the scene, store these in kd-tree
   3. update the progressive radiance estimates
   4. evaluate radiance
 - Make some more scenes to display caustic effects