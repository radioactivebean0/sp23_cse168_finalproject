#pragma once

// CMake insert NDEBUG when building with RelWithDebInfo
// This is an ugly hack to undo that...
#undef NDEBUG

#include "scene.h"

// stores a hitpoint, not final maybe need more or less fields
struct PPMHitPoint {
    Vector3 position;
    Vector3 normal;
    Vector3 omeganot;
    Vector3 emission;
    Vector3 beta; // throughput of the color multiplied over the bounces
    material_e mat;
    Real r; // photon radius
    Real n; // photon count accumulated
    Vector3 tau; // accumulated flux
};

// grid of vectors which are lists of the hit points sampled from a given pixel
// useful for parallelization
struct PPMGrid {
    PPMGrid(int w, int h, int spp) : 
            width(w),
            height(h),
            spp(spp){
        ppm_grid.resize(w * h);
        for (auto hit_point_vec: ppm_grid) {
            hit_point_vec = std::vector<PPMHitPoint>(); // initialize PPM hit points
        }
    }

    std::vector<PPMHitPoint> &operator()(int x, int y) {
        return ppm_grid[y * width + x];
    }

    int width;
    int height;
    int spp;

    std::vector< std::vector<PPMHitPoint> > ppm_grid;
};

// TODO: trace a ray from the camera to the scene to get visible point
PPMHitPoint generate_visible_point(
    const Scene &scene,
    const Vector3 &ray,
    pcg32_state &pcg_state,
    int max_depth,
    const Real default_radius
);