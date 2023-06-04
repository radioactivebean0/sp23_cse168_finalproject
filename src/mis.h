#pragma once
#include "scene.h"

Vector3 mis_path_trace(const Scene &scene, const Vector3 &ray, const Vector3 &ray_origin, pcg32_state &pcg_state, int max_depth);