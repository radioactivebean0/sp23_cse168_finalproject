#include "scene.h"
#include "bvh.h"
#include <unordered_map>
Scene::Scene(const ParsedScene &scene) :
        camera(from_parsed_camera(scene.camera)),
        width(scene.camera.width),
        height(scene.camera.height),
        background_color(scene.background_color),
        samples_per_pixel(scene.samples_per_pixel){
    // Extract triangle meshes from the parsed scene.
    int tri_mesh_count = 0;
    for (const ParsedShape &parsed_shape : scene.shapes) {
        if (std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            tri_mesh_count++;
        }
    }
    meshes.resize(tri_mesh_count);
    // Extract the shapes
    std::unordered_map<int, std::vector<int>> arealight_map;
    tri_mesh_count = 0;
    for (int i = 0; i < (int)scene.shapes.size(); i++) {
        const ParsedShape &parsed_shape = scene.shapes[i];
        if (auto *sph = std::get_if<ParsedSphere>(&parsed_shape)) {
            shapes.push_back( Sphere{
                .center = sph->position,
                .material_id = sph->material_id,
                .radius = sph->radius,
                .area_light_id = sph->area_light_id
            });
            if (std::get<Sphere>(shapes.back()).area_light_id != -1){
                arealight_map.emplace(std::make_pair(i,std::vector<int>{(int)(shapes.size()-1)}));
            }
        } else if (auto *parsed_mesh = std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            meshes[tri_mesh_count] = TriangleMesh{.material_id= parsed_mesh->material_id,
                                                  .positions= parsed_mesh->positions,
                                                  .indices= parsed_mesh->indices, 
                                                  .normals = parsed_mesh->normals,
                                                  .uvs = parsed_mesh->uvs,
                                                  .area_light_id= parsed_mesh->area_light_id};
            // Extract all the individual triangles
            if (meshes[tri_mesh_count].area_light_id!=-1){
                arealight_map.emplace(std::make_pair(i, std::vector<int>{}));
            }
            for (int face_index = 0; face_index < (int)parsed_mesh->indices.size(); face_index++) {
                shapes.push_back(Triangle{.face_index = face_index, .mesh = &meshes[tri_mesh_count]});
                if (meshes[tri_mesh_count].area_light_id!=-1){
                    arealight_map[i].push_back(shapes.size()-1);
                }
            }
            tri_mesh_count++;
        } else {
            assert(false);
        }
    }
    // Copy the materials
    for (const ParsedMaterial &parsed_mat : scene.materials) {
        if (auto *diffuse = std::get_if<ParsedDiffuse>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&diffuse->reflectance)){
                materials.push_back(Material{.material_type = material_e::DiffuseType, .reflectance = SolidTexture{.reflectance = *tex}});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&diffuse->reflectance)){
                materials.push_back(Material{.material_type = material_e::DiffuseType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset)});
            }
        } else if (auto *mirror = std::get_if<ParsedMirror>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&mirror->reflectance)){
                materials.push_back(Material{.material_type = material_e::MirrorType, .reflectance = SolidTexture{.reflectance = *tex}});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&mirror->reflectance)){
                materials.push_back(Material{.material_type = material_e::MirrorType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset)});
            }
        } else if (auto *plastic = std::get_if<ParsedPlastic>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&plastic->reflectance)){
                materials.push_back(Material{.material_type = material_e::PlasticType, .reflectance = SolidTexture{.reflectance = *tex}, .ref_index = plastic->eta});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&plastic->reflectance)){
                materials.push_back(Material{.material_type = material_e::PlasticType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .ref_index = plastic->eta});
            }
        } else if (auto *phong = std::get_if<ParsedPhong>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&phong->reflectance)){
                materials.push_back(Material{.material_type = material_e::PhongType, .reflectance = SolidTexture{.reflectance = *tex}, .exponent=phong->exponent});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&phong->reflectance)){
                materials.push_back(Material{.material_type = material_e::PhongType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .exponent=phong->exponent});
            }
        } else if (auto *blinnphong = std::get_if<ParsedBlinnPhong>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&blinnphong->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongType, .reflectance = SolidTexture{.reflectance = *tex}, .exponent=blinnphong->exponent});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&blinnphong->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .exponent=blinnphong->exponent});
            }
        } else if (auto *blinnphongmicrofacet = std::get_if<ParsedBlinnPhongMicrofacet>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&blinnphongmicrofacet->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongMicrofacetType, .reflectance = SolidTexture{.reflectance = *tex}, .exponent=blinnphongmicrofacet->exponent});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&blinnphongmicrofacet->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongMicrofacetType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .exponent=blinnphongmicrofacet->exponent});
            }
        } else {
            assert(false);
        }
    }
    // Copy the lights
    for (const ParsedLight &parsed_light : scene.lights) {
        // We assume all lights are point lights for now.
        if (auto *point_light = std::get_if<ParsedPointLight>(&parsed_light)){
            lights.push_back(PointLight{.position= point_light->position , .intensity= point_light->intensity});
        } else if (auto *area_light = std::get_if<ParsedDiffuseAreaLight>(&parsed_light)){
            for (int s: arealight_map.at(area_light->shape_id)){
                lights.push_back(AreaLight{.radiance = area_light->radiance, .shape_idx = s});
            }
        } else {
            assert(false);
        }
    }
}

Vector3 get_diffuse_color(const Scene &scene, const Vector3 &pt, const Real eps, const Vector3 &kd, const Vector3 &n, const bool with_simd, pcg32_state &pcg_state){
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
            if (hit_bvh(scene.bvh, light_ray, light_pos, eps, eps, infinity<Real>(), &temp_shape, t_shadow, uv, with_simd)){
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
                if (hit_bvh(scene.bvh, light_ray, light_pos, eps, eps, infinity<Real>(), &temp_shape, t_shadow, uv, with_simd)){
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
                if (hit_bvh(scene.bvh, light_ray, light_pos, eps, eps, infinity<Real>(), &temp_shape, t_shadow, uv, with_simd)){
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

Vector3 get_mirror_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                            const Real eps, const Vector3 &kd, const Vector3 &n,
                            const bool with_simd, const bool shading_norms, const bool fresnel,
                            pcg32_state &pcg_state){
    // get reflected ray intersecting object and then get color
    Real t_reflect = -1.0;
    Vector3 ray_reflect =  ray_in - 2.0*dot(ray_in,n)*n;
    Shape *temp_shape;
    Vector2 uv;
    if (hit_bvh(scene.bvh, ray_reflect, pt, eps, eps, infinity<Real>(), &temp_shape, t_reflect, uv, with_simd)){
        Vector3 mirror_pt = pt + t_reflect * ray_reflect;
        if (fresnel){
            return (kd + (1.0-kd)* pow((1.0 - dot(n,ray_reflect)), 5)) * radiance(scene, ray_reflect, mirror_pt, eps, temp_shape, uv, with_simd, shading_norms, fresnel, pcg_state);
        } else {
            return kd * radiance(scene, ray_reflect, mirror_pt, eps, temp_shape, uv, with_simd, shading_norms, fresnel, pcg_state);
        }
    } else {
        if (fresnel){
            return (kd + (1.0-kd)* pow((1.0 - dot(n,ray_reflect)), 5)) * scene.background_color;
        } else {
            return kd*scene.background_color;
        }
    }
}

Vector3 get_plastic_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                            const Real eps, const Vector3 &kd, const Real eta, const Vector3 &n,
                            const bool with_simd, const bool shading_norms, const bool fresnel,
                            pcg32_state &pcg_state){
    Real t_reflect = -1.0;
    Vector3 ray_reflect =  ray_in - 2.0*dot(ray_in,n)*n;
    Shape *temp_shape;
    Vector2 uv;
    const Real fnot = pow((eta - 1)/(eta+1),2);
    const Real F = fnot + ((1.0 - fnot) * pow(1 - dot(n,ray_reflect), 5));
    if (hit_bvh(scene.bvh, ray_reflect, pt, eps, eps, infinity<Real>(), &temp_shape, t_reflect, uv, with_simd)){
        Vector3 mirror_pt = pt + t_reflect * ray_reflect;
        return (F*radiance(scene, ray_reflect, mirror_pt, eps, temp_shape, uv, with_simd, shading_norms, true, pcg_state) )+ ((Vector3{1.0,1.0,1.0}-F)*get_diffuse_color(scene, pt, eps, kd, n, true, pcg_state));
    } else {
        return (F*scene.background_color) + ((Vector3{1.0,1.0,1.0}-F)*get_diffuse_color(scene, pt, eps, kd, n, true, pcg_state));
    }
}

Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms, const bool fresnel, pcg32_state &pcg_state){
    Vector3 kd, n;
    Vector2 uvt;
    if (auto *sph = std::get_if<Sphere>(shape)){ // get sphere color
        // calc uv
        uvt = sphere_uv(sph, pt);
        kd = get_texture_kd(scene.materials.at(sph->material_id).reflectance, uvt);
        n = normalize(pt - sph->center);
        if (scene.materials.at(sph->material_id).material_type == material_e::MirrorType) {
            return get_mirror_color(scene, ray_in, pt, eps, kd, n, with_simd, shading_norms, fresnel, pcg_state);
        } else if (scene.materials.at(sph->material_id).material_type == material_e::PlasticType){
            Real eta = scene.materials.at(sph->material_id).ref_index;
            return get_plastic_color(scene, ray_in, pt, eps, kd, eta, n, with_simd, shading_norms, fresnel, pcg_state);
        } else {
            return get_diffuse_color(scene, pt, eps, kd, n, with_simd, pcg_state);
        }
    } else if (auto *tri = std::get_if<Triangle>(shape)) { // triangle intersection test
        Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
        Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
        Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
        // calc uv
        uvt = triangle_uv(tri, uv);
        // shading normals
        if (shading_norms){
            n = shading_norm(tri, uv);
        }else {
            n = normalize(cross(p1 - p0, p2 - p1));
        }
        if (dot(n, ray_in) > 0.0){
            n = -n;
        }
        kd = get_texture_kd(scene.materials.at(tri->mesh->material_id).reflectance, uvt);
        if (scene.materials.at(tri->mesh->material_id).material_type == material_e::MirrorType) {
            return get_mirror_color(scene, ray_in, pt, eps, kd, n, with_simd, shading_norms, fresnel, pcg_state);
        } else if (scene.materials.at(tri->mesh->material_id).material_type == material_e::PlasticType){
            Real eta = scene.materials.at(tri->mesh->material_id).ref_index;
            return get_plastic_color(scene, ray_in, pt, eps, kd, eta, n, with_simd, shading_norms, fresnel, pcg_state);
        } else {
            return get_diffuse_color(scene, pt, eps, kd, n, with_simd, pcg_state);
        }
    } else {
        assert(false);
    }
}
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms, const bool fresnel){
    pcg32_state t_state = pcg32_state{};
    return get_color(scene, ray_in, pt, eps, shape, uv, with_simd, shading_norms, fresnel, t_state);
}
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms){
    return get_color(scene, ray_in, pt, eps, shape, uv, with_simd, shading_norms, false);
}
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd){
    return get_color(scene, ray_in, pt, eps, shape, uv, with_simd, false, false);
}
Vector3 get_color(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt, const Real eps, Shape *shape, const Vector2 &uv){
    return get_color(scene, ray_in, pt, eps, shape, uv, true);
}

Vector3 radiance(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt,
                    const Real eps, Shape *shape, const Vector2 &uv, const bool with_simd,
                    const bool shading_norms, const bool fresnel, pcg32_state &pcg_state){
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
    return get_color(scene, ray_in, pt, eps, shape, uv, with_simd, shading_norms, fresnel, pcg_state);
}

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

bool closest_hit(const std::vector<Shape> &shapes, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, int &obj_id, Real &t_val){
    t_val = -1.0;
    obj_id = -1;
    int shape;
    Real t_temp;
    Vector3 intersection, p0, p1, p2;
    Vector2 uv;
    for (shape = 0; shape < shapes.size(); ++shape){
        if (auto *sph = std::get_if<Sphere>(&shapes.at(shape))){ // sphere intersection test
            t_temp = hit_sphere(*sph, ray, ray_origin);
            if (t_temp >= 0.0 && (t_val < 0 || t_val > t_temp) && distance(ray_origin, ray_origin + ray * t_temp)> eps){ // new closest obj
                t_val = t_temp;
                obj_id = shape;
            }
        } else if (auto *tri = std::get_if<Triangle>(&shapes.at(shape))) { // triangle intersection test
            p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
            p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
            p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
            t_temp = get_tri_intersect(p0, p1, p2, ray, ray_origin, eps, uv);

            if (t_temp >= 0.0 && (t_val < 0 || t_val > t_temp) && distance(ray_origin, ray_origin + ray * t_temp)> eps){ // new closest obj
                t_val = t_temp;
                obj_id = shape;
            }
        }
    }
    if (t_val < 0.0){
        return false;
    } else {
        return true;
    }
}

Vector3 get_diffuse_color_nobvh(const Scene &scene, const Vector3 &pt, const Real eps, const Vector3 &kd, const Vector3 &n){
    Vector3 light_pos, light_ray, temp_b;
    Real t_shadow, hit_dist;
    int temp_obj;
    Vector3 l, I;
    Real d_squared, nl;
    Vector3 color{0.0,0.0,0.0};
    for (int lit = 0; lit < scene.lights.size(); lit++){ // iterate all lights
        if (auto *light = std::get_if<PointLight>(&scene.lights.at(lit))){
            light_pos = light->position;
            light_ray = normalize(pt-light_pos);
            t_shadow = -1.0;
            if (closest_hit(scene.shapes, light_ray, light_pos, eps, temp_obj, t_shadow)){
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
                        nl = dot(Real(-1.0)*n,l);
                    }
                    color += ((kd * nl)/c_PI)*(I/d_squared);
                }
            } else {
                continue;
            }
        }
    }
    return color;
}

Vector3 get_mirror_color_nobvh(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt, const Real eps, const Vector3 &kd, const Vector3 &n){
    // get reflected ray intersecting object and then get color
    Real t_reflect = -1.0;
    Vector3 ray_reflect =  ray_in - 2.0*dot(ray_in,n)*n;
    int mirror_obj = -1;
    Vector3 b_coord;
        if (closest_hit(scene.shapes, ray_reflect, pt, eps, mirror_obj, t_reflect)){
            Vector3 mirror_pt = pt + t_reflect * ray_reflect;
            return kd*get_color_nobvh(scene, ray_reflect, mirror_pt, eps, mirror_obj);
        } else {
            return kd*scene.background_color;
        }
}

Vector3 get_color_nobvh(const Scene &scene, const Vector3 &ray_in, const Vector3 &pt, const Real eps, const int obj_id){
    Vector3 kd, n;
    if (auto *sph = std::get_if<Sphere>(&scene.shapes.at(obj_id))){ // get sphere color
        kd = std::get<SolidTexture>(scene.materials.at(sph->material_id).reflectance).reflectance;
        n = normalize(pt - sph->center);
        if (scene.materials.at(sph->material_id).material_type == material_e::MirrorType) {
            return get_mirror_color_nobvh(scene, ray_in, pt, eps, kd, n);
        } else {
            return get_diffuse_color_nobvh(scene, pt, eps, kd, n);
        }
    } else if (auto *tri = std::get_if<Triangle>(&scene.shapes.at(obj_id))) { // triangle intersection test
        Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
        Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
        Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
        n = normalize(cross(p1 - p0, p2 - p1));
        kd = std::get<SolidTexture>(scene.materials.at(tri->mesh->material_id).reflectance).reflectance;
        if (scene.materials.at(tri->mesh->material_id).material_type == material_e::MirrorType) {
            return get_mirror_color_nobvh(scene, ray_in, pt, eps, kd, n);
        } else {
            return get_diffuse_color_nobvh(scene, pt, eps, kd, n);
        }
    } else {
        assert(false);
    }
}

inline Vector3 rand_cos(pcg32_state &pcg_state){
    Real r1 = next_pcg32_real<Real>(pcg_state);
    Real r2 = next_pcg32_real<Real>(pcg_state);
    Real z = sqrt(1-r2);

    Real phi = c_TWOPI*r1;
    Real x = cos(phi)*sqrt(r2);
    Real y = sin(phi)*sqrt(r2);


    return Vector3{x, y, z};
}

inline Vector3 rand_phong_cos(pcg32_state &pcg_state, Real alpha){
    Real u1 = next_pcg32_real<Real>(pcg_state);
    Real u2 = next_pcg32_real<Real>(pcg_state);
    Real theta = acos(pow(1.0 - u1, 1/(alpha+1)));
    Real phi = c_TWOPI*u2;

    Real x = sin(theta)*cos(phi);
    Real y = sin(theta)*sin(phi);


    return normalize(Vector3{x, y, cos(theta)});
}

// returning pdf of -1.0 means pure specular
Vector3 brdf_sample(const Scene &scene, const Vector3 &ray, const material_e mat, const int mat_id, const Vector3 &sn, const Vector3 &gn, pcg32_state &pcg_state, int &specular){
    if (mat == material_e::MirrorType){ // mirrors
            specular = 1;
            return ray - 2.0*dot(ray, gn)*gn;
    } else if (mat == material_e::PlasticType) { // plastic stocastically choose diffuse and reflect
        const Vector3 ray_reflect = ray - 2.0*dot(ray, gn)*gn;
        const Real fnot = pow((scene.materials.at(mat_id).ref_index - 1)/( scene.materials.at(mat_id).ref_index +1),2);
        const Real F = fnot + ((1.0 - fnot) * pow(1 - dot(gn,ray_reflect), 5));
        if (next_pcg32_real<Real>(pcg_state) < F){
            specular = 1;
            return ray_reflect;
        } else {
            Vector3 scatter = rand_cos(pcg_state);
            Vector3 a = Vector3{1.0,0.0,0.0};
            if (abs(gn.x) > 0.9){
                a = Vector3{0.0,1.0,0.0};
            }
            Vector3 t = normalize(cross(a,gn));
            Vector3 s = cross(t,gn);
            Vector3 bounce = scatter.x*s + scatter.y*t + scatter.z*gn;
            specular = 0;
            return bounce;
        }
    } else if (mat == material_e::PhongType) { // phong material
        const Vector3 ray_reflect = ray - 2.0*dot(ray, sn)*sn;
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 scatter = rand_phong_cos(pcg_state, exponent);
        Vector3 a = Vector3{1.0,0.0,0.0};
        if (abs(ray_reflect.x) > 0.9){
            a = Vector3{0.0,1.0,0.0};
        }
        Vector3 t = normalize(cross(a,ray_reflect));
        Vector3 s = cross(t,ray_reflect);
        Vector3 omeganot = scatter.x*s + scatter.y*t + scatter.z*ray_reflect;
        specular = 0;
        return omeganot;
    } else if (mat == material_e::BlinnPhongType){ // blinn phong
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 scatter = rand_phong_cos(pcg_state, exponent);
        Vector3 a = Vector3{1.0,0.0,0.0};
        if (abs(sn.x) > 0.9){
            a = Vector3{0.0,1.0,0.0};
        }
        Vector3 t = normalize(cross(a,sn));
        Vector3 s = cross(t,sn);
        Vector3 halfvec = scatter.x*s + scatter.y*t + scatter.z*sn;
        specular = 0;
        Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
        return omeganot;
    } else if (mat == material_e::BlinnPhongMicrofacetType){
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 scatter = rand_phong_cos(pcg_state, exponent);
        Vector3 a = Vector3{1.0,0.0,0.0};
        if (abs(sn.x) > 0.9){
            a = Vector3{0.0,1.0,0.0};
        }
        Vector3 t = normalize(cross(a,sn));
        Vector3 s = cross(t,sn);
        Vector3 halfvec = scatter.x*s + scatter.y*t + scatter.z*sn;
        specular = 0;
        Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
        return omeganot;
    } else {        // scattering, cosine hemisphere sampling and diffuse
        Vector3 scatter = rand_cos(pcg_state);
        Vector3 a = Vector3{1.0,0.0,0.0};
        if (abs(gn.x) > 0.9){
            a = Vector3{0.0,1.0,0.0};
        }
        Vector3 t = normalize(cross(a,gn));
        Vector3 s = cross(t,gn);
        Vector3 bounce = scatter.x*s + scatter.y*t + scatter.z*gn;
        specular = 0;
        return bounce;
    }
}

bool brdf_eval(const Scene &scene, const Vector3 &ray, const Vector3 &omeganot, const material_e mat, const int mat_id, const Vector3 &sn, const Vector3 &gn, const Vector3 &kd, pcg32_state &pcg_state, Vector3 &value, Real &pdf){
    if (mat == material_e::PhongType) { // phong material
        const Vector3 ray_reflect = ray - 2.0*dot(ray, sn)*sn;
        Real exponent = scene.materials.at(mat_id).exponent;
        // Real nwo = dot(sn, omeganot);
        // if (nwo <=0.0){
        //     return false;
        // }
        Real nwo = dot(ray_reflect, omeganot);
        if (nwo <= 0.0){
            return false;
        }
        nwo = pow(nwo, exponent);
        value = (kd*nwo*((exponent+1)*c_INVTWOPI));
        pdf = ((exponent+1)*c_INVTWOPI*(pow(dot(ray_reflect, omeganot), exponent)));
        return true;
    } else if (mat == material_e::BlinnPhongType){ // blinn phong
        // sample omega is half vec
        const Vector3 halfvec = normalize(-ray + omeganot);
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
        // if (dot(sn,omeganot)<0.0){
        //     return false;
        // }
        Vector3 fh = kd + (1.0-kd)*pow(1-dot(halfvec,omeganot),5);
        pdf = ((exponent+1.0)*pow(dot(sn,halfvec),exponent))/(c_TWOPI*4.0*dot(omeganot,halfvec));
        Real normc = (exponent + 2)/(c_FOURPI*(2.0 - pow(2.0, -exponent/2.0)));
        value = normc*fh*pow(dot(sn,halfvec),exponent);
        return true;
    } else if (mat == material_e::BlinnPhongMicrofacetType){
        const Vector3 halfvec = normalize(-ray + omeganot);
        Real exponent = scene.materials.at(mat_id).exponent;
        // if (dot(sn,omeganot)<0.0){
        //     return false;
        // }
        pdf = ((exponent+1.0)*pow(dot(sn,halfvec),exponent))/(c_TWOPI*4.0*dot(omeganot,halfvec));
        Vector3 fh = kd + (1.0-kd)*pow(1-dot(halfvec,omeganot),5);
        Real d = (exponent + 2)/(c_TWOPI)*pow(dot(sn,halfvec),exponent);
        Real alphanot = sqrt(0.5*exponent + 1.0)/tan(acos(dot(omeganot,sn)));
        Real alphai = sqrt(0.5*exponent + 1.0)/tan(acos(dot(-ray,sn)));
        Real Gnot; Real Gi;
        Real G;
        if (alphanot < 1.6){
            Gnot = (3.535*alphanot + 2.181*pow(alphanot,2.0))/(1+2.276*alphanot+2.577*pow(alphanot,2.0));
        } else {
            Gnot = 1.0;
        }
        if (alphai < 1.6){
            Gi = (3.535*alphai + 2.181*pow(alphai,2.0))/(1+2.276*alphai+2.577*pow(alphai,2.0));
        } else {
            Gi = 1.0;
        }
        G = Gnot*Gi;
        value = (fh*d*G)/(4*dot(sn,-ray));
        return true;
    } else {        // scattering, cosine hemisphere sampling and diffuse
        const Vector3 &scatter = omeganot;
        Real nwo = dot(sn,scatter);
        if (nwo <= 0.0){
            return false;
        }
        pdf = (dot(gn,scatter)*c_INVPI);
        value = (kd*nwo*c_INVPI);
        return true;
    }
}
//pdf -1.0 means no hit;
Vector3 light_sample(const Scene &scene, const Vector3 &ray, const Vector3 &pt, const Real &eps, const material_e mat, const int mat_id, const Vector3 &sn, const Vector3 &gn, pcg32_state &pcg_state){
    Vector3 light_pos, light_ray;
    Vector3 l,I;
    Vector3 color{0.0,0.0,0.0};
    Vector2 uv;
    int lit = next_pcg32_real<Real>(pcg_state)*scene.lights.size();
    if (auto *alight = std::get_if<AreaLight>(&scene.lights.at(lit))){
        if (auto *sph = std::get_if<Sphere>(&scene.shapes.at(alight->shape_idx))){
            // sample a point on the sphere
            Real theta = acos(1.0-(2.0*next_pcg32_real<double>(pcg_state)));
            Real phi = c_TWOPI * next_pcg32_real<double>(pcg_state);
            light_pos = sph->center + Vector3{sph->radius*sin(theta)*cos(phi), sph->radius*sin(theta)*sin(phi),sph->radius*cos(theta)};
            light_ray = normalize(pt-light_pos);
            return -light_ray;
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
            return -light_ray;
        } else {
            assert(false);
        }
    } else {
        assert (false);
    }
}

Real light_pdf(const Scene &scene, const Vector3 &omeganot, const Vector3 &pt){
    Real d_squared;
    Vector2 uv;
    const Real eps = 0.000001;
    for (int lit = 0; lit < scene.lights.size(); lit++){ // iterate all lights
        if (auto *alight = std::get_if<AreaLight>(&scene.lights.at(lit))){
            // do area light sampling here
            if (auto *sph = std::get_if<Sphere>(&scene.shapes.at(alight->shape_idx))){
                Real t_temp = hit_sphere(*sph, omeganot, pt);
                if (t_temp > 0.0){
                    Vector3 hit_pt = pt + t_temp*omeganot;
                    d_squared = distance_squared(pt, hit_pt);
                    // Real cos_theta_max = sqrt(1 - (sph->radius*sph->radius)/length_squared(sph->center-pt));
                    // Real solid_angle = c_TWOPI*(1-cos_theta_max);
                    // return  d_squared*1.0 / solid_angle;
                    Real cosine = dot(omeganot,normalize(sph->center-hit_pt));
                    Real pdf = 2.0*d_squared/(cosine*(c_FOURPI*(sph->radius*sph->radius)));
                    return pdf;
                } else {
                    continue;
                }
            } else if (auto *tri = std::get_if<Triangle>(&scene.shapes.at(alight->shape_idx))){
                Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
                Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
                Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
                Real t_shadow = get_tri_intersect(p0, p1, p2, omeganot, pt, eps, uv);
                if (t_shadow > 0.0){
                    Vector3 hit_pt = pt + t_shadow*omeganot;
                    //std::cout << "asdfasdfasd" << std::endl;
                    d_squared = distance_squared(pt, hit_pt);
                    Vector3 norm =  shading_norm(tri, uv);
                    Real cosine = dot(omeganot,-norm);
                    return d_squared/(cosine*(0.5*length(cross(p2-p0,p1-p0))));
                } else {
                    continue;
                }
            } else {
                assert(false);
            }
        } else {
            assert(false);
        }
    }
    return 0.0;
}
Vector3 mis_path_trace(const Scene &scene, const Vector3 &ray, const Vector3 &ray_origin, pcg32_state &pcg_state, int max_depth){
    if (max_depth <= 0){return Vector3{0.0,0.0,0.0};}
    const Real eps = 0.00001;
    Shape *hs;
    Real t_val = -1.0;
    Vector2 uv;

    if (hit_cbvh(scene.cbvh, ray, ray_origin, eps, eps, infinity<Real>(), &hs, t_val, uv)){ // was a hit
        Vector3 pt = ray_origin + t_val * ray;
        Vector3 emission = Vector3{0.0,0.0,0.0};
        Vector3 kd;
        Vector3 gn;
        Vector3 sn;
        material_e mat;
        int mat_id;
        if (auto *sph = std::get_if<Sphere>(hs)){
            // check if area light
            if (sph->area_light_id != -1){
                emission = std::get<AreaLight>(scene.lights.at(sph->area_light_id)).radiance;
            }
            // get the color
            uv = sphere_uv(sph, pt);
            gn = normalize(pt - sph->center);
            sn = gn;
            kd = get_texture_kd(scene.materials.at(sph->material_id).reflectance, uv);
            mat_id = sph->material_id;
        } else if (auto *tri = std::get_if<Triangle>(hs)) {
            // check if area light
            if (tri->mesh->area_light_id != -1){
                emission =  std::get<AreaLight>(scene.lights.at(tri->mesh->area_light_id)).radiance;
            }
            // calc uv
            Vector2 uvt = triangle_uv(tri, uv);
            Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
            Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
            Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
            gn = normalize(cross(p1 - p0, p2 - p1));
            // shading normals
            sn = shading_norm(tri, uv);
            if (dot(sn, ray) > 0.0){
                sn = -sn;
            }
            if (dot(gn, sn) < 0.0){
                gn = -gn;
            }
            kd = get_texture_kd(scene.materials.at(tri->mesh->material_id).reflectance, uvt);
            mat_id = tri->mesh->material_id;
        } else {
            assert(false);
        }
        mat = scene.materials.at(mat_id).material_type;
        // skip direct lighting
        // brdf calculations
        int spec;
        Vector3 omeganot = brdf_sample(scene, ray, mat, mat_id, sn, gn, pcg_state, spec);
        if (spec==1){ // purely specular
            return emission + (kd + (1.0-kd)* pow((1.0 - dot(gn,omeganot)), 5))*mis_path_trace(scene, omeganot, pt, pcg_state, max_depth-1);
        } else {
            Vector3 val;
            Real brdfpdf;
            Vector3 light_ray = light_sample(scene, ray, pt, eps, mat, mat_id, sn, gn, pcg_state);
            Real coinflip = next_pcg32_real<Real>(pcg_state);
            if (coinflip < 0.5){ // use light ray
                if (brdf_eval(scene, ray, light_ray, mat, mat_id, sn, gn, kd, pcg_state, val, brdfpdf)){
                    Real combinedpdf = 0.5*brdfpdf + 0.5*light_pdf(scene, light_ray, pt);
                    Vector3 temp_vec = val*mis_path_trace(scene, light_ray, pt, pcg_state, max_depth-1)/combinedpdf;
                    return emission + temp_vec;
                } else {
                    return emission;
                }
            } else {
                if (brdf_eval(scene, ray, omeganot, mat, mat_id, sn, gn, kd, pcg_state, val, brdfpdf)){
                    Real combinedpdf = 0.5*brdfpdf + 0.5*light_pdf(scene, omeganot, pt);
                    return emission + val*mis_path_trace(scene, omeganot, pt, pcg_state, max_depth-1)/combinedpdf;
                } else {
                    return emission;
                }
            }
        }
    } else {
        return scene.background_color;
    }
}

Vector3 path_trace(const Scene &scene, const Vector3 &ray, const Vector3 &ray_origin, pcg32_state &pcg_state, int max_depth){
    if (max_depth == 0){return Vector3{0.0,0.0,0.0};}
    const Real eps = 0.00001;
    Shape *hs;
    Real t_val = -1.0;
    Vector2 uv;

    if (hit_cbvh(scene.cbvh, ray, ray_origin, eps, eps, infinity<Real>(), &hs, t_val, uv)){ // was a hit
        Vector3 pt = ray_origin + t_val * ray;
        Vector3 emission = Vector3{0.0,0.0,0.0};
        Vector3 kd;
        Vector3 gn;
        Vector3 sn;
        material_e mat;
        int mat_id;
        if (auto *sph = std::get_if<Sphere>(hs)){
            // check if area light
            if (sph->area_light_id != -1){
                emission = std::get<AreaLight>(scene.lights.at(sph->area_light_id)).radiance;
            }
            // get the color
            uv = sphere_uv(sph, pt);
            gn = normalize(pt - sph->center);
            sn = gn;
            kd = get_texture_kd(scene.materials.at(sph->material_id).reflectance, uv);
            mat_id = sph->material_id;
        } else if (auto *tri = std::get_if<Triangle>(hs)) {
            // check if area light
            if (tri->mesh->area_light_id != -1){
                emission =  std::get<AreaLight>(scene.lights.at(tri->mesh->area_light_id)).radiance;
            }
            // calc uv
            Vector2 uvt = triangle_uv(tri, uv);
            Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
            Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
            Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
            gn = normalize(cross(p1 - p0, p2 - p1));
            // shading normals
            sn = shading_norm(tri, uv);
            if (dot(sn, ray) > 0.0){
                sn = -sn;
            }
            if (dot(gn, sn) < 0.0){
                gn = -gn;
            }
            kd = get_texture_kd(scene.materials.at(tri->mesh->material_id).reflectance, uvt);
            mat_id = tri->mesh->material_id;
        } else {
            assert(false);
        }
        mat = scene.materials.at(mat_id).material_type;
        // skip direct lighting
        // brdf calculations
        if (mat == material_e::MirrorType){ // mirrors
            const Vector3 ray_reflect = ray - 2.0*dot(ray, gn)*gn;
            return emission + (kd + (1.0-kd)* pow((1.0 - dot(gn,ray_reflect)), 5)) * path_trace(scene, ray_reflect, pt, pcg_state, max_depth-1);
        } else if (mat == material_e::PlasticType) { // plastic stocastically choose diffuse and reflect
            const Vector3 ray_reflect = ray - 2.0*dot(ray, gn)*gn;
            const Real fnot = pow((scene.materials.at(mat_id).ref_index - 1)/( scene.materials.at(mat_id).ref_index +1),2);
            const Real F = fnot + ((1.0 - fnot) * pow(1 - dot(gn,ray_reflect), 5));
            if (next_pcg32_real<Real>(pcg_state) < F){
                return emission + path_trace(scene, ray_reflect, pt, pcg_state, max_depth-1);
            } else {
                Vector3 scatter = rand_cos(pcg_state);
                Vector3 a = Vector3{1.0,0.0,0.0};
                if (abs(gn.x) > 0.9){
                    a = Vector3{0.0,1.0,0.0};
                }
                Vector3 t = normalize(cross(a,gn));
                Vector3 s = cross(t,gn);
                Vector3 bounce = scatter.x*s + scatter.y*t + scatter.z*gn;
                Real nwo = dot(sn,bounce);
                if (nwo < 0.0){
                    return emission;
                }
                return emission + (kd*nwo*c_INVPI)*(path_trace(scene, bounce, pt, pcg_state, max_depth-1))/(dot(gn,bounce)*c_INVPI);
            }
        } else if (mat == material_e::PhongType) { // phong material
            const Vector3 ray_reflect = ray - 2.0*dot(ray, sn)*sn;
            Real exponent = scene.materials.at(mat_id).exponent;
            Vector3 scatter = rand_phong_cos(pcg_state, exponent);
            Vector3 a = Vector3{1.0,0.0,0.0};
            if (abs(ray_reflect.x) > 0.9){
                a = Vector3{0.0,1.0,0.0};
            }
            Vector3 t = normalize(cross(a,ray_reflect));
            Vector3 s = cross(t,ray_reflect);
            Vector3 omeganot = scatter.x*s + scatter.y*t + scatter.z*ray_reflect;
            Real nwo = dot(sn, omeganot);
            if (nwo < 0.0){
                return emission;
            }
            nwo = dot(ray_reflect, omeganot);
            if (nwo < 0.0){
                return emission;
            }
            nwo = pow(nwo, exponent);
            return emission + (kd*nwo*((exponent+1)*c_INVTWOPI))*path_trace(scene, omeganot, pt, pcg_state, max_depth-1)/((exponent+1)*c_INVTWOPI*(pow(dot(ray_reflect, omeganot), exponent)));
        } else if (mat == material_e::BlinnPhongType){ // blinn phong
            Real exponent = scene.materials.at(mat_id).exponent;
            Vector3 scatter = rand_phong_cos(pcg_state, exponent);
            Vector3 a = Vector3{1.0,0.0,0.0};
            if (abs(sn.x) > 0.9){
                a = Vector3{0.0,1.0,0.0};
            }
            Vector3 t = normalize(cross(a,sn));
            Vector3 s = cross(t,sn);
            Vector3 halfvec = scatter.x*s + scatter.y*t + scatter.z*sn;
            Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
            if (dot(sn,omeganot)<0.0){
                return emission;
            }
            Vector3 fh = kd + (1.0-kd)*pow(1-dot(halfvec,omeganot),5);
            Real pdf = ((exponent+1.0)*pow(dot(sn,halfvec),exponent))/(c_TWOPI*4.0*dot(omeganot,halfvec));
            Real normc = (exponent + 2)/(c_FOURPI*(2.0 - pow(2.0, -exponent/2.0)));
            return emission + normc*fh*pow(dot(sn,halfvec),exponent)*path_trace(scene, omeganot, pt, pcg_state, max_depth-1)/pdf;
        } else if (mat == material_e::BlinnPhongMicrofacetType){
            Real exponent = scene.materials.at(mat_id).exponent;
            Vector3 scatter = rand_phong_cos(pcg_state, exponent);
            Vector3 a = Vector3{1.0,0.0,0.0};
            if (abs(sn.x) > 0.9){
                a = Vector3{0.0,1.0,0.0};
            }
            Vector3 t = normalize(cross(a,sn));
            Vector3 s = cross(t,sn);
            Vector3 halfvec = scatter.x*s + scatter.y*t + scatter.z*sn;
            Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
            if (dot(sn,omeganot)<0.0){
                return emission;
            }
            Real pdf = ((exponent+1.0)*pow(dot(sn,halfvec),exponent))/(c_TWOPI*4.0*dot(omeganot,halfvec));
            Vector3 fh = kd + (1.0-kd)*pow(1-dot(halfvec,omeganot),5);
            Real d =(exponent + 2)/(c_TWOPI)*pow(dot(sn,halfvec),exponent);
            Real alphanot = sqrt(0.5*exponent + 1.0)/tan(acos(dot(omeganot,sn)));
            Real alphai = sqrt(0.5*exponent + 1.0)/tan(acos(dot(-ray,sn)));
            Real Gnot; Real Gi;
            Real G;
            if (alphanot < 1.6){
                Gnot = (3.535*alphanot + 2.181*pow(alphanot,2.0))/(1+2.276*alphanot+2.577*pow(alphanot,2.0));
            } else {
                Gnot = 1.0;
            }
            if (alphai < 1.6){
                Gi = (3.535*alphai + 2.181*pow(alphai,2.0))/(1+2.276*alphai+2.577*pow(alphai,2.0));
            } else {
                Gi = 1.0;
            }
            G = Gnot*Gi;

            return emission + (fh*d*G)/(4*dot(sn,-ray)) *path_trace(scene, omeganot, pt, pcg_state, max_depth-1)/pdf;
        } else {        // scattering, cosine hemisphere sampling and diffuse
            Vector3 scatter = rand_cos(pcg_state);
            Vector3 a = Vector3{1.0,0.0,0.0};
            if (abs(gn.x) > 0.9){
                a = Vector3{0.0,1.0,0.0};
            }
            Vector3 t = normalize(cross(a,gn));
            Vector3 s = cross(t,gn);
            Vector3 bounce = scatter.x*s + scatter.y*t + scatter.z*gn;
            Real nwo = dot(sn,bounce);
            if (nwo < 0.0){
                return emission;
            }
            return emission + (kd*nwo*c_INVPI)*(path_trace(scene, bounce, pt, pcg_state, max_depth-1))/(dot(gn,bounce)*c_INVPI);
        }
    } else {
        return scene.background_color;
    }
}