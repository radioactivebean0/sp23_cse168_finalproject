#pragma once
#include "vector.h"
#include "parse_scene.h"
#include "3rdparty/pcg.h"

struct Camera {
    Vector3 lookfrom;
    Vector3 lookat;
    Vector3 up;
    Real vfov;
    Vector3 llc, horizontal, vertical;
    int width, height;
};

Camera from_parsed_camera(const ParsedCamera &camera);
Vector3 get_ray(Camera &camera, int x, int y, pcg32_state &pcg_state);
