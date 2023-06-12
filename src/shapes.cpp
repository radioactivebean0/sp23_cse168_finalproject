#include "shapes.h"

double hit_sphere(const Sphere &sphere, const Vector3 &ray, const Vector3 &origin){
    Real a = dot(ray, ray);
    Real b = Real(2.0) * dot(origin-sphere.center, ray);
    Real c = dot(origin-sphere.center, origin-sphere.center) - sphere.radius*sphere.radius;
    Real disc = b*b - 4*a*c;
    if (disc < 0){
        return -1.0;
    } else {
        Real min_disc = (-b - sqrt(disc)) / (Real(2.0)*a);
        if (min_disc < 0.0){
            return (-b + sqrt(disc)) / (Real(2.0)*a);
        } else {
            return min_disc;
        }
    }
}

double get_tri_intersect(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &ray, const Vector3 &lookfrom, const Real eps, Vector2 &uv) {
        Vector3 vec1 = p1 - p0;
        Vector3 vec2 = p2 - p0;
        Vector3 pvec = cross(ray, vec2);
        Real det = dot(vec1, pvec);

        if (abs(det) < eps) {// ray is parallel
            return -1.0;
        }
        Vector3 tvec = lookfrom-p0;
        Real m = dot(tvec, pvec) * (1.0/det);
        if (m < 0.0 || m > 1.0){
            return -1.0;
        }
        Vector3 qvec = cross(tvec, vec1);
        Real n = (1.0/det) * dot(ray, qvec);

        if (n < 0.0 || m+n > 1.0){
            return -1.0;
        }

        Real t = (1.0/det) * dot(vec2, qvec);
        if (t > eps){
            uv.x = m;
            uv.y = n;
            return t;
        } else {
            return -1.0;
        }
}

double hit_shape(const Shape* shape, const Vector3 &ray, const Vector3 &origin, const Real eps, Vector2 &uv){
    double t_val = -1.0;
    if (auto *sph = std::get_if<Sphere>(shape)) {
        t_val = hit_sphere(*sph, ray, origin);
    } else if (auto *tri = std::get_if<Triangle>(shape)) {
        Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
        Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
        Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
        t_val = get_tri_intersect(p0, p1, p2, ray, origin, eps, uv);
    } else {
        assert(false);
    }
    if (distance(origin, origin + ray * t_val)> eps){
        return t_val;
    } else {
        return -1.0;
    }
}
double hit_shape(const Shape* shape, const Vector3 &ray, const Vector3 &origin, const Real eps){
    Vector2 uv;
    return hit_shape(shape, ray, origin, eps, uv);
}


Vector2 sphere_uv(const Sphere* sphere, const Vector3 &pt){
    Vector3 n_coords = (pt - sphere->center)/sphere->radius;
    Real theta = acos(-n_coords.y);
    Real phi = atan2(-n_coords.z, n_coords.x)+c_PI;
    return Vector2{phi/c_TWOPI, 1.0-(theta/c_PI)};
}

Vector2 triangle_uv(const Triangle* triangle, const Vector2 &uv){
    if (triangle->mesh->uvs.size()==0){
        // use b_coord
        return uv;
    }
    else {
        Vector3i indices = triangle->mesh->indices.at(triangle->face_index);
        Vector2 u0 = triangle->mesh->uvs.at(indices.x);
        Vector2 u1 = triangle->mesh->uvs.at(indices.y);
        Vector2 u2 = triangle->mesh->uvs.at(indices.z);
        return (1.0- uv.x - uv.y)*u0 + uv.x*u1 + uv.y * u2;
    }
}

Vector3 shading_norm(const Triangle* tri, const Vector2 &uv){
    Vector3i indices = tri->mesh->indices.at(tri->face_index);
    if (tri->mesh->normals.size()==0){
        Vector3 p0 = tri->mesh->positions.at(indices.x);
        Vector3 p1 = tri->mesh->positions.at(indices.y);
        Vector3 p2 = tri->mesh->positions.at(indices.z);
        return normalize(cross(p1 - p0, p2 - p1));
    } else {
        Vector3 n0 = tri->mesh->normals.at(indices.x);
        Vector3 n1 = tri->mesh->normals.at(indices.y);
        Vector3 n2 = tri->mesh->normals.at(indices.z);
        return normalize((1.0- uv.x - uv.y)*n0 + uv.x*n1 + uv.y * n2);
    }
}

void triangle_points(const Triangle* tri, Vector3 &p0, Vector3 &p1, Vector3 &p2){
    p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
    p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
    p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
}

Vector3 sample_shape_point(const Shape* shape, pcg32_state &pcg_state, Vector3 &norm){
    if (auto *sph = std::get_if<Sphere>(shape)){
        Real theta = acos(1.0-(2.0*next_pcg32_real<double>(pcg_state)));
        Real phi = c_TWOPI * next_pcg32_real<double>(pcg_state);
        Vector3 pt = sph->center + Vector3{sph->radius*sin(theta)*cos(phi), sph->radius*sin(theta)*sin(phi),sph->radius*cos(theta)};
        norm = (pt - sph->center) / sph->radius;
        return pt;
    } else if (auto *tri = std::get_if<Triangle>(shape)){
        Vector3 p0, p1, p2;
        triangle_points(tri, p0, p1, p2);
        Real u1 = next_pcg32_real<double>(pcg_state);
        Real u2 = next_pcg32_real<double>(pcg_state);
        Real b1 = 1 - sqrt(u1);
        Real b2 = u2 * sqrt(u1);
        norm = shading_norm(tri, Vector2{b1,b2});
        return (1-b1-b2)*p0 + b1*p1 + b2*p2;
    } else {
        assert(false);
    }
}

int shape_matid(const Shape* shape){
    if (auto *sph = std::get_if<Sphere>(shape)){
        return sph->material_id;
    } else if (auto *tri = std::get_if<Triangle>(shape)){
        return tri->mesh->material_id;
    } else {
        assert(false);
    }
}
Vector3 shape_shade_norm(const Shape* shape, const Vector3 &pt, const Vector2 &uv){
    if (auto *sph = std::get_if<Sphere>(shape)){
        return (pt - sph->center) / sph->radius;
    } else if (auto *tri = std::get_if<Triangle>(shape)){
        return shading_norm(tri, uv);
    } else {
        assert(false);
    }
}
Vector3 shape_geo_norm(const Shape* shape, const Vector3 &pt, const Vector2 &uv){
    if (auto *sph = std::get_if<Sphere>(shape)){
        return (pt - sph->center) / sph->radius;
    } else if (auto *tri = std::get_if<Triangle>(shape)){
        Vector3i indices = tri->mesh->indices.at(tri->face_index);
        Vector3 p0 = tri->mesh->positions.at(indices.x);
        Vector3 p1 = tri->mesh->positions.at(indices.y);
        Vector3 p2 = tri->mesh->positions.at(indices.z);
        return normalize(cross(p1 - p0, p2 - p1));
    } else {
        assert(false);
    }
}