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
