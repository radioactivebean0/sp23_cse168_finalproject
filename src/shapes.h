#pragma once
#include "vector.h"
#include <vector>
#include <variant>
#include "3rdparty/pcg.h"

struct HW2TriangleMesh {
    int material_id = -1;
    std::vector<Vector3> positions;
    std::vector<Vector3i> indices;
};

struct TriangleMesh {
    int material_id = -1;
    int area_light_id = -1;
    std::vector<Vector3> positions;
    std::vector<Vector3i> indices;
    std::vector<Vector3> normals;
    std::vector<Vector2> uvs;
};

struct Sphere {
    Vector3 center;
    Real radius;
    int material_id;
    int area_light_id;
};

struct Triangle{
    int face_index;
    const TriangleMesh *mesh;
};
struct HW2Triangle{
    int face_index;
    const HW2TriangleMesh *mesh;
};
using Shape = std::variant<Sphere, Triangle>;

struct AABB {
    Vector3 a, b;
    int axis;
    std::vector<Shape*> shapes;
    std::shared_ptr<AABB> left, right;
};

struct compact_AABB {
    Vector3 a, b;
    uint8_t axis;
    int next;
    Shape* shape;
};

struct WideAABB {
    WideAABB();
    WideAABB(std::shared_ptr<WideAABB> aabb);
    Vector3 a, b;
    int axis;
    std::vector<Shape*> shapes;
    std::shared_ptr<WideAABB> c1, c2, c3, c4;
};

Vector2 sphere_uv(const Sphere* sphere, const Vector3 &pt);
Vector2 triangle_uv(const Triangle* triangle, const Vector2 &uv);
Vector3 shading_norm(const Triangle* tri, const Vector2 &uv);
Vector3 sample_shape_point(const Shape* shape, pcg32_state &pcg_state, Vector3 &norm, Real &pdf);
Vector3 shape_shade_norm(const Shape* shape, const Vector3 &pt, const Vector2 &uv);
Vector3 shape_geo_norm(const Shape* shape, const Vector3 &pt, const Vector2 &uv);
int shape_matid(const Shape* shape);

void triangle_points(const Triangle* tri, Vector3 &p0, Vector3 &p1, Vector3 &p2);

double hit_sphere(const Sphere & sphere, const Vector3 &ray, const Vector3 &origin);
double get_tri_intersect(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &ray, const Vector3 &lookfrom, const Real eps, Vector2 &uv);
double hit_shape(const Shape* shape, const Vector3 &ray, const Vector3 &origin, const Real eps, Vector2 &uv);
double hit_shape(const Shape* shape, const Vector3 &ray, const Vector3 &origin, const Real eps);
