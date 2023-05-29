#pragma once
#include "vector.h"
#include <variant>

struct PointLight {
    Vector3 position;
    Vector3 intensity;
};


struct AreaLight {
    int shape_idx;
    Vector3 radiance;
};

using Light = std::variant<PointLight, AreaLight>;