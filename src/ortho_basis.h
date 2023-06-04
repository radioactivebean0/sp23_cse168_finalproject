#pragma once
#include "vector.h"


inline Vector3 ortho_basis(const Vector3 &sampled, const Vector3 &z_axis){
    const Vector3 a = (abs(z_axis.x) > 0.9) ? Vector3{0.0,1.0,0.0} : Vector3{1.0,0.0,0.0};
    const Vector3 t = normalize(cross(a,z_axis));
    const Vector3 s = cross(t,z_axis);
    return sampled.x*s + sampled.y*t + sampled.z*z_axis;
}