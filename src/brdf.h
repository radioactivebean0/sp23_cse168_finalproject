#include "scene.h"
// todo, move this out
Vector3 brdf_bounce(const Scene &scene, Shape *hs, const Vector3 &omeganot, Vector3 &hit_pt, const Vector2 &uv, pcg32_state &pcg_state);
Vector3 brdf_color(const Scene &scene, Shape *hs, const Vector3 &omeganot, const Vector3 &hit_pt, const Vector2 &uv);
