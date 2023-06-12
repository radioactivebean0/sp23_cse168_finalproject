#pragma once
#include "3rdparty/pcg.h"
#include "torrey.h"
#include "vector.h"


inline Vector3 rand_phong_cos(pcg32_state &pcg_state, Real alpha){
    Real u1 = next_pcg32_real<Real>(pcg_state);
    Real u2 = next_pcg32_real<Real>(pcg_state);
    Real theta = acos(pow(1.0 - u1, Real(1)/(alpha+1)));
    Real phi = c_TWOPI*u2;

    Real x = sin(theta)*cos(phi);
    Real y = sin(theta)*sin(phi);


    return normalize(Vector3{x, y, cos(theta)});
}

inline Vector3 rand_cos(pcg32_state &pcg_state){
    Real r1 = next_pcg32_real<Real>(pcg_state);
    Real r2 = next_pcg32_real<Real>(pcg_state);
    Real z = sqrt(1.0-r2);

    Real phi = c_TWOPI*r1;
    Real x = cos(phi)*sqrt(r2);
    Real y = sin(phi)*sqrt(r2);


    return Vector3{x, y, z};
}

inline Vector3 rand_uniform_sphere(pcg32_state &pcg_state){
    Real r1 = next_pcg32_real<Real>(pcg_state);
    Real r2 = next_pcg32_real<Real>(pcg_state);
    return Vector3 {cos(c_TWOPI*r1)*2.0*sqrt(r2*(1.0-r2)),
                    sin(c_TWOPI*r1)*2.0*sqrt(r2*(1.0-r2)),
                    1.0 - 2.0*r2 };
}