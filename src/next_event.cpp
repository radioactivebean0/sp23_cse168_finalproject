#include "next_event.h"
#include "bvh.h"

Vector3 get_diffuse_color_v2(const Scene &scene, const Vector3 &pt, const Real eps, const Vector3 &kd, const Vector3 &n, pcg32_state &pcg_state){
    Vector3 light_pos, light_ray, temp_b;
    Real t_shadow, hit_dist;
    Shape *temp_shape;
    Vector3 l, I;
    Real d_squared, nl;
    Vector3 color{0.0,0.0,0.0};
    Vector2 uv;
    for (int lit = 0; lit < scene.lights.size(); lit++){ // iterate all lights
        if (auto *light = std::get_if<PointLight>(&scene.lights.at(lit))){
            light_pos = light->position;
            light_ray = normalize(pt-light_pos);
            t_shadow = -1.0;
            if (hit_cbvh(scene.cbvh, light_ray, light_pos, eps, eps, infinity<Real>(), &temp_shape, t_shadow, uv)){
                hit_dist = distance(light_pos, light_pos + light_ray * t_shadow);
                if (hit_dist > eps &&
                        hit_dist < (1-eps)*distance(light_pos, pt)){
                    continue;
                } else {
                    l = Real(-1.0) * light_ray;
                    I = light->intensity;
                    d_squared = distance_squared(pt, light_pos);
                    nl = dot(n, l);
                    if (nl < 0.0){
                        continue;
                    }
                    color += ((kd * nl)/c_PI)*(I/d_squared);
                }
            }
        } else if (auto *alight = std::get_if<AreaLight>(&scene.lights.at(lit))){
            // do area light sampling here
            if (auto *sph = std::get_if<Sphere>(&scene.shapes.at(alight->shape_idx))){
                // sample a point on the sphere
                Real theta = acos(1.0-(2.0*next_pcg32_real<double>(pcg_state)));
                Real phi = c_TWOPI * next_pcg32_real<double>(pcg_state);
                light_pos = sph->center + Vector3{sph->radius*sin(theta)*cos(phi), sph->radius*sin(theta)*sin(phi),sph->radius*cos(theta)};
                light_ray = normalize(pt-light_pos);
                t_shadow = -1.0;
                if (hit_cbvh(scene.cbvh, light_ray, light_pos, eps, eps, infinity<Real>(), &temp_shape, t_shadow, uv)){
                    hit_dist = distance(light_pos, light_pos + light_ray * t_shadow);
                    if (hit_dist > eps &&
                            hit_dist < (1-eps)*distance(light_pos, pt)){
                        continue;
                    } else {
                        l = Real(-1.0) * light_ray;
                        I = alight->radiance;
                        d_squared = distance_squared(pt, light_pos);
                        Real nsl = dot(n, l);
                        if (nsl < 0.0){
                            continue;
                            //nl = dot(Real(-1.0)*n,l);
                        }
                        Real nxl = dot(normalize(sph->center - light_pos),l);
                        if (nxl < 0.0){
                            continue;
                        }
                        color += ((kd * nsl)/c_PI)*((I*nxl)/d_squared)*(c_FOURPI*(sph->radius*sph->radius));
                    }
                }
            } else if (auto *tri = std::get_if<Triangle>(&scene.shapes.at(alight->shape_idx))){
                Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
                Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
                Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
                Real u1 = next_pcg32_real<double>(pcg_state);
                Real u2 = next_pcg32_real<double>(pcg_state);
                Real b1 = 1 - sqrt(u1);
                Real b2 = u2 * sqrt(u1);
                light_pos =(1-b1-b2)*p0 + b1*p1 + b2*p2;
                light_ray = normalize(pt-light_pos);
                t_shadow = -1.0;
                if (hit_cbvh(scene.cbvh, light_ray, light_pos, eps, eps, infinity<Real>(), &temp_shape, t_shadow, uv)){
                    hit_dist = distance(light_pos, light_pos + light_ray * t_shadow);
                    if (hit_dist > eps &&
                        hit_dist < (1-eps)*distance(light_pos, pt)){
                        continue;
                    } else {
                        l = Real(-1.0) * light_ray;
                        I = alight->radiance;
                        d_squared = distance_squared(pt, light_pos);
                        Real nsl = dot(n, l);
                        if (nsl < 0.0){
                            continue;
                            //nl = dot(Real(-1.0)*n,l);
                        }
                        Vector3 shadenorm = shading_norm(tri, Vector2{b1,b2});
                        Vector3 nx = normalize(cross(p1 - p0, p2 - p1));
                        if (dot(nx, shadenorm) < 0.0){
                            nx = - nx;
                        }
                        Real nxl = dot(-nx, l);
                        if (nxl < 0.0){
                            continue;
                        }
                        color += (((kd * nsl)/c_PI)*((I*nxl)/d_squared))*(0.5*length(cross(p2-p0,p1-p0)));
                    }
                }
            } else {
                assert(false);
            }
        } else {
            assert(false);
        }
    }
    return color;
}

Vector3 get_mirror_color_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                            const Real eps, const Vector3 &kd, const Vector3 &n,
                            pcg32_state &pcg_state){
    // get reflected ray intersecting object and then get color
    Real t_reflect = -1.0;
    Vector3 ray_reflect =  ray_in - 2.0*dot(ray_in,n)*n;
    Shape *temp_shape;
    Vector2 uv;
    if (hit_cbvh(scene.cbvh, ray_reflect, pt, eps, eps, infinity<Real>(), &temp_shape, t_reflect, uv)){
        Vector3 mirror_pt = pt + t_reflect * ray_reflect;
        return (kd + (1.0-kd)* pow((1.0 - dot(n,ray_reflect)), 5)) * radiance_v2(scene, ray_reflect, mirror_pt, eps, temp_shape, uv, pcg_state);
    } else {
        return (kd + (1.0-kd)* pow((1.0 - dot(n,ray_reflect)), 5)) * scene.background_color;
    }
}

Vector3 get_plastic_color_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                            const Real eps, const Vector3 &kd, const Real eta, const Vector3 &n,
                            pcg32_state &pcg_state){
    Real t_reflect = -1.0;
    Vector3 ray_reflect =  ray_in - 2.0*dot(ray_in,n)*n;
    Shape *temp_shape;
    Vector2 uv;
    const Real fnot = pow((eta - 1)/(eta+1),2);
    const Real F = fnot + ((1.0 - fnot) * pow(1 - dot(n,ray_reflect), 5));
    if (hit_cbvh(scene.cbvh, ray_reflect, pt, eps, eps, infinity<Real>(), &temp_shape, t_reflect, uv)){
        Vector3 mirror_pt = pt + t_reflect * ray_reflect;
        return (F*radiance_v2(scene, ray_reflect, mirror_pt, eps, temp_shape, uv, pcg_state) )+ ((Vector3{1.0,1.0,1.0}-F)*get_diffuse_color_v2(scene, pt, eps, kd, n, pcg_state));
    } else {
        return (F*scene.background_color) + ((Vector3{1.0,1.0,1.0}-F)*get_diffuse_color_v2(scene, pt, eps, kd, n, pcg_state));
    }
}

Vector3 get_color_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, pcg32_state &pcg_state){
    Vector3 kd, n;
    Vector2 uvt;
    if (auto *sph = std::get_if<Sphere>(shape)){ // get sphere color
        // calc uv
        uvt = sphere_uv(sph, pt);
        kd = get_texture_kd(scene.materials.at(sph->material_id).reflectance, uvt);
        n = normalize(pt - sph->center);
        if (scene.materials.at(sph->material_id).material_type == material_e::MirrorType) {
            return get_mirror_color_v2(scene, ray_in, pt, eps, kd, n, pcg_state);
        } else if (scene.materials.at(sph->material_id).material_type == material_e::PlasticType){
            Real eta = scene.materials.at(sph->material_id).ref_index;
            return get_plastic_color_v2(scene, ray_in, pt, eps, kd, eta, n, pcg_state);
        } else {
            return get_diffuse_color_v2(scene, pt, eps, kd, n, pcg_state);
        }
    } else if (auto *tri = std::get_if<Triangle>(shape)) { // triangle intersection test
        // calc uv
        uvt = triangle_uv(tri, uv);
        // shading normals
        n = shading_norm(tri, uv);
        if (dot(n, ray_in) > 0.0){
            n = -n;
        }
        kd = get_texture_kd(scene.materials.at(tri->mesh->material_id).reflectance, uvt);
        if (scene.materials.at(tri->mesh->material_id).material_type == material_e::MirrorType) {
            return get_mirror_color_v2(scene, ray_in, pt, eps, kd, n, pcg_state);
        } else if (scene.materials.at(tri->mesh->material_id).material_type == material_e::PlasticType){
            Real eta = scene.materials.at(tri->mesh->material_id).ref_index;
            return get_plastic_color_v2(scene, ray_in, pt, eps, kd, eta, n, pcg_state);
        } else {
            return get_diffuse_color_v2(scene, pt, eps, kd, n, pcg_state);
        }
    } else {
        assert(false);
    }
}

Vector3 radiance_v2(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, pcg32_state &pcg_state){
    if (auto *sph = std::get_if<Sphere>(shape)){ // get sphere color
        // check if area light
        if (sph->area_light_id != -1){
            return std::get<AreaLight>(scene.lights.at(sph->area_light_id)).radiance;
        }
    } else if (auto *tri = std::get_if<Triangle>(shape)) { // triangle intersection test
        if (tri->mesh->area_light_id != -1){
            return std::get<AreaLight>(scene.lights.at(tri->mesh->area_light_id)).radiance;
        }
    } else {
        assert(false);
    }
    return get_color_v2(scene, ray_in, pt, eps, shape, uv, pcg_state);

}