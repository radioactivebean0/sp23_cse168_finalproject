# Final Project Proposal: Photon Mapping
I want to implement photon mapping to better model caustics with global illumination. Caustics create natural and beautiful patterns but are unable to be easily rendered using standard global illumination path tracing since the patterns are most easily observed with small light sources passing through glass which is not modeled well due to zero probability of path traced rays going through the glass and successfully hitting a light source.

## Rough outline of tasks
 - Find scenes/objects: I need to find or make some scenes and objects with glass that are good at showing the effects of caustics and improvements using photon mapping versus the standard path tracing methods.
 - implement glass materials: The current path tracing algorithm I have implemented does not model dielectric materials. There would need to be changes done to both the path tracer and also the parsers to support new material types.
 - photon mapping: the main goal of the project. I'll start with reading the PBR chaper on photon mapping. https://www.pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/Stochastic_Progressive_Photon_Mapping

 ## What I have currently
 The current starting point of the project is the accumulation of work done on the homework assignments.
 - path tracing for global illumination
 - diffuse, mirror, phong, blinn-phong, blinn-phong microfacet BRDFs for materials
 - one sample multiple importance sampling
 - anti-aliasing
 - compacted bounding volume heirarchies built with some surface area heuristic with some SIMD and using iterative stack based methods and some SIMD

## Additional improvements I'd like to implement if time permits
 - better optimize bvh with SIMD and improve the surface area heuristic implementation
 - implement iterative path tracing instead of recursive path tracing
 - if there is a lot of time then maybe even implementing combination of irradiance caching and photon mapping