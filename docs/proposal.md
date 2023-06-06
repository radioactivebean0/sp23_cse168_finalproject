## Final Project Proposal: Photon Mapping

Ryan and I want to implement photon mapping to better model caustics with global illumination.
Caustics create natural and beautiful patterns but are unable to be easily rendered using standard
global illumination path tracing since the patterns are most easily observed with small light
sources passing through glass which is not modeled well due to zero probability of path traced rays
going through the glass and successfully hitting a light source.

## Rough outline of tasks

- Find scenes/objects: we need to find or make some scenes and objects with glass that are good
  at showing the effects of caustics and improvements using photon mapping versus the standard path
  tracing methods.
- Implement glass materials: The current path tracing algorithm Ryan has implemented does not model
  dielectric materials. There would need to be changes done to both the path tracer and also the
  parsers to support new material types.
- Photon mapping: the main goal of the project. We'll start by reading the
  [PBR chapter](https://www.pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/Stochastic_Progressive_Photon_Mapping) on photon mapping.

## What We have currently

The current starting point of the project is the accumulation of work done on the homework
assignments.

- Path tracing for global illumination
- Diffuse, mirror, Phong, Blinn-Phong, Blinn-Phong microfacet BRDFs for materials
- One sample multiple importance sampling
- Anti-aliasing
- Compacted bounding volume heirarchies built with some surface area heuristic with some SIMD and
  using iterative stack based methods and some SIMD

## Additional improvements we'd like to implement if time permits

- Better optimized BVH with SIMD and improve the surface area heuristic implementation
- Implement iterative path tracing instead of recursive path tracing
- If there is a lot of time, maybe even implementing combination of irradiance caching and photon
  mapping
