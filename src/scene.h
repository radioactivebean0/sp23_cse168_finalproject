#pragma once
#include "vector.h"
#include "parse_scene.h"
#include "shapes.h"
#include "materials.h"
#include "lights.h"
#include "camera.h"

struct Scene {
    Scene(const ParsedScene &scene);
    Scene(const ParsedScene &scene, bool with_cbvh);

    Camera camera;
    int width, height;
    std::vector<Shape> shapes;
    std::vector<Material> materials;
    std::vector<Light> lights;
    Vector3 background_color;
    int samples_per_pixel;
    AABB bvh;
    //std::shared_ptr<WideAABB> wbvh;
    // For the Triangle in the shapes to reference to.
    std::vector<TriangleMesh> meshes;
    compact_AABB *cbvh;
};

Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv);
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd);
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms);
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms, const bool fresnel);
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms, const bool fresnel, pcg32_state &pcg_state);
Vector3 get_color_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, pcg32_state &pcg_state);

bool closest_hit(const std::vector<Shape> &shapes, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, int &obj_id, Real &t_val);

Vector3 get_color_nobvh(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt, const Real eps, const int obj_id);

Vector3 radiance(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms, const bool fresnel, pcg32_state &pcg_state);

Vector3 radiance_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, pcg32_state &pcg_state);

Vector3 path_trace(const Scene &scene, const Vector3 &ray, const Vector3 &ray_origin, pcg32_state &pcg_state, int max_depth);
Vector3 mis_path_trace(const Scene &scene, const Vector3 &ray, const Vector3 &ray_origin, pcg32_state &pcg_state, int max_depth);