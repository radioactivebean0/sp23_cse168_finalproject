#pragma once
#include "scene.h"


Vector3 get_color_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, pcg32_state &pcg_state);
Vector3 radiance_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, pcg32_state &pcg_state);
