#pragma once

#include "shapes.h"
#include "hw3_scene.h"



inline bool wide_aabb_x_compare (const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b) {
    return a->a.x < b->a.x;
}
inline bool wide_aabb_y_compare (const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b) {
    return a->a.y < b->a.y;
}
inline bool wide_aabb_z_compare (const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b) {
    return a->a.z < b->a.z;
}

inline bool wide_aabb_x_compare_max (const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b) {
    return a->b.x < b->b.x;
}
inline bool wide_aabb_y_compare_max (const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b) {
    return a->b.y < b->b.y;
}
inline bool wide_aabb_z_compare_max (const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b) {
    return a->b.z < b->b.z;
}

std::shared_ptr<WideAABB> build_wide_bvh(std::vector<Shape> &shapes);

bool hit_bvh(const WideAABB &bvh, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, Shape **hit_shape);