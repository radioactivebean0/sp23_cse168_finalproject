#include "ppm.h"
#include "bvh.h"

PPMHitPoint generate_visible_point(const Scene &scene, const Vector3 &ray, pcg32_state &pcg_state, int max_depth, const Real default_radius){
    int depth = 0;
    const Real eps = 0.00001;
    Vector3 current_ray = ray;
    Vector3 current_origin = scene.camera.lookfrom;
    Vector3 transmission = Vector3{1.0,1.0,1.0};
    Shape *hs;
    Real t_val = -1.0;
    Vector2 uv;
    PPMHitPoint vis_pt;
    vis_pt.r=default_radius;
    vis_pt.n = 0;
    vis_pt.tau = Vector3{0.0,0.0,0.0};
    Vector3 pt, sn, kd;
    int mat_id;
    while (depth < max_depth){
        if (hit_cbvh(scene.cbvh, current_ray, current_origin, eps, eps, infinity<Real>(), &hs, t_val, uv)){
            // was a hit
            pt = current_origin + t_val * current_ray;
            Vector3 emission = Vector3{0.0,0.0,0.0};
            if (auto *sph = std::get_if<Sphere>(hs)){
                // check if area light
                if (sph->area_light_id != -1){
                    emission = std::get<AreaLight>(scene.lights.at(sph->area_light_id)).radiance;
                }
                // get the color
                uv = sphere_uv(sph, pt);
                sn = normalize(pt - sph->center);
                kd = get_texture_kd(scene.materials.at(sph->material_id).reflectance, uv);
                mat_id = sph->material_id;
            } else if (auto *tri = std::get_if<Triangle>(hs)) {
                // check if area light
                if (tri->mesh->area_light_id != -1){
                    emission = std::get<AreaLight>(scene.lights.at(tri->mesh->area_light_id)).radiance;
                }
                // calc uv
                Vector2 uvt = triangle_uv(tri, uv);
                Vector3 p0, p1, p2;
                triangle_points(tri, p0, p1, p2);
                // shading normals
                sn = shading_norm(tri, uv);
                if (dot(sn, current_ray) > 0.0){
                    sn = -sn;
                }
                kd = get_texture_kd(scene.materials.at(tri->mesh->material_id).reflectance, uvt);
                mat_id = tri->mesh->material_id;
            } else {
                std::cerr << "Error: not sphere or triangle" << std::endl;
                exit(1);
            }
            material_e mat = scene.materials.at(mat_id).material_type;
            if (mat == material_e::DiffuseType){ // hit diffuse, create a visible point
                vis_pt.beta = transmission*kd;
                vis_pt.emission = emission;
                vis_pt.normal = sn;
                vis_pt.position = pt;
                vis_pt.omeganot = current_ray;
                vis_pt.mat = material_e::DiffuseType;
                return vis_pt;
            } else if (mat == material_e::MirrorType) { // reflect
                Vector3 ray_reflect =  current_ray - 2.0*dot(current_ray, sn)*sn;
                transmission *= (kd + (1.0-kd)* pow((1.0 - dot(sn,ray_reflect)), 5));
                vis_pt.beta = transmission;
                vis_pt.emission = emission;
                vis_pt.normal = sn;
                vis_pt.position = pt;
                vis_pt.omeganot = current_ray;
                vis_pt.mat = material_e::MirrorType;
                current_origin = pt;
                current_ray = current_ray - 2.0*dot(current_ray, sn)*sn;
                depth++;
            } else if (mat == material_e::DielectricType) { // reflect or refract
                const Real inIR = scene.materials.at(mat_id).ref_index;
                const Real outIR = scene.materials.at(mat_id).exponent;
                Vector3 n = sn;
                Real cos_theta = fmin(dot(-current_ray, sn), 1.0);
                Real ref_ratio = (cos_theta > 0.0) ? outIR/inIR : inIR/outIR; // depends if ray points in or out
                if (cos_theta < 0.0){ // ray inside
                    cos_theta = -cos_theta;
                    n = -n;
                }
                Real sin_theta = sqrt(Real(1.0) - cos_theta*cos_theta);
                Vector3 direction;
                Real schlick_fresnel = pow((ref_ratio-1.0)/(ref_ratio+1.0), 2.0);
                schlick_fresnel = schlick_fresnel + (1.0-schlick_fresnel)*pow(1.0-cos_theta, 5);
                if (ref_ratio* sin_theta > 1.0 || schlick_fresnel > next_pcg32_real<Real>(pcg_state)){ // reflect
                    current_origin = pt;
                    current_ray = current_ray - 2.0*dot(current_ray,n)*n;
                } else { // refract
                    Vector3 r_out_perp = ref_ratio * (current_ray + cos_theta*n);
                    Vector3 r_out_parallel = -sqrt(fabs(1.0 - length_squared(r_out_perp))) * n;
                    current_origin = current_origin + (t_val+eps) * current_ray; // move to avoid moire patterns
                    current_ray = r_out_perp + r_out_parallel;
                }
                vis_pt.beta = transmission;
                vis_pt.emission = emission;
                vis_pt.normal = n;
                vis_pt.position = pt;
                vis_pt.omeganot = current_ray;
                vis_pt.mat = material_e::DielectricType;
                depth++;
            } else {
                std::cerr << "Error: unsupported material" << std::endl;
                exit(1);
            }
        } else {
            // was not a hit
            if (depth == 0){
                return PPMHitPoint{
                    .position=Vector3{0.0,0.0,0.0},
                    .normal=Vector3{0.0,0.0,0.0},
                    .omeganot=Vector3{0.0,0.0,0.0},
                    .emission=Vector3{0.0,0.0,0.0},
                    .beta=scene.background_color,
                    .mat=material_e::DiffuseType, // just use diffuse as default, shouldnt check this anyways
                    .r=-1.0,
                    .n=0,
                    .tau=Vector3{0.0,0.0,0.0}
                };
            } else {
                return vis_pt;
            }
        }
    }
    return PPMHitPoint{ // the case where no diffuse surface is found before reaching max depth
        .position=Vector3{0.0,0.0,0.0},
        .normal=Vector3{0.0,0.0,0.0},
        .omeganot=Vector3{0.0,0.0,0.0},
        .emission=Vector3{0.0,0.0,0.0},
        .beta=scene.background_color,
        .mat=material_e::DiffuseType, // just use diffuse as default, shouldnt check this anyways
        .r=-1.0,
        .n=0,
        .tau=Vector3{0.0,0.0,0.0}
    };
}