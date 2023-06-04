#include "mis.h"
#include "bvh.h"
#include "sampling.h"
#include "ortho_basis.h"

// returning pdf of -1.0 means pure specular
Vector3 brdf_sample(const Scene &scene, const Vector3 &ray, const material_e mat, const int mat_id, const Vector3 &sn, const Vector3 &gn, pcg32_state &pcg_state, int &specular){
    if (mat == material_e::MirrorType){ // mirrors
        specular = 1;
        return ray - 2.0*dot(ray, gn)*gn;
    } else if (mat == material_e::DielectricType){
        const Real inIR = scene.materials.at(mat_id).ref_index;
        const Real outIR = scene.materials.at(mat_id).exponent;
        Vector3 n = sn;
        Real cos_theta = fmin(dot(-ray, sn), 1.0);
        Real ref_ratio = (cos_theta > 0.0) ? outIR/inIR : inIR/outIR; // depends if ray points in or out
        if (cos_theta < 0.0){ // ray inside
            cos_theta = -cos_theta;
            //ref_ratio = outIR/inIR;
            n = -n;
        }
        Real sin_theta = sqrt(Real(1.0) - cos_theta*cos_theta);
        Vector3 direction;
        Real schlick_fresnel = pow((ref_ratio-1.0)/(ref_ratio+1.0), 2.0);
        schlick_fresnel = schlick_fresnel + (1.0-schlick_fresnel)*pow(1.0-cos_theta, 5);
        if (ref_ratio* sin_theta > 1.0 || schlick_fresnel > next_pcg32_real<Real>(pcg_state)){ // reflect
            specular = 2;
            direction = ray - 2.0*dot(ray,n)*n;
        } else { // refract
            specular = 3;
            Vector3 r_out_perp = ref_ratio * (ray + cos_theta*n);
            Vector3 r_out_parallel = -sqrt(fabs(1.0 - length_squared(r_out_perp))) * n;
            direction = r_out_perp + r_out_parallel;
        }
        return direction;
    } else if (mat == material_e::PlasticType) { // plastic stocastically choose diffuse and reflect
        const Vector3 ray_reflect = ray - 2.0*dot(ray, gn)*gn;
        const Real fnot = pow((scene.materials.at(mat_id).ref_index - 1)/( scene.materials.at(mat_id).ref_index +1),2);
        const Real F = fnot + ((1.0 - fnot) * pow(1 - dot(gn,ray_reflect), 5));
        if (next_pcg32_real<Real>(pcg_state) < F){
            specular = 1;
            return ray_reflect;
        } else {
            Vector3 scatter = rand_cos(pcg_state);
            Vector3 bounce = ortho_basis(scatter, gn);
            specular = 0;
            return bounce;
        }
    } else if (mat == material_e::PhongType) { // phong material
        const Vector3 ray_reflect = ray - 2.0*dot(ray, sn)*sn;
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 scatter = rand_phong_cos(pcg_state, exponent);
        Vector3 omeganot = ortho_basis(scatter, ray_reflect);
        specular = 0;
        return omeganot;
    } else if (mat == material_e::BlinnPhongType){ // blinn phong
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 scatter = rand_phong_cos(pcg_state, exponent);
        Vector3 halfvec = ortho_basis(scatter, sn);
        specular = 0;
        Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
        return omeganot;
    } else if (mat == material_e::BlinnPhongMicrofacetType){
        Real exponent = scene.materials.at(mat_id).exponent;
        Vector3 scatter = rand_phong_cos(pcg_state, exponent);
        Vector3 halfvec = ortho_basis(scatter, sn);
        specular = 0;
        Vector3 omeganot = ray - 2.0*dot(ray, halfvec)*halfvec;
        return omeganot;
    } else {        // scattering, cosine hemisphere sampling and diffuse
        Vector3 scatter = rand_cos(pcg_state);
        Vector3 bounce = ortho_basis(scatter, gn);
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
            Vector3 p0, p1, p2;
            triangle_points(tri, p0, p1, p2);
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
    Real pdf = 0.0;
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
                    pdf += (2.0*d_squared/(cosine*(c_FOURPI*(sph->radius*sph->radius))))/scene.lights.size();
                } else {
                    continue;
                }
            } else if (auto *tri = std::get_if<Triangle>(&scene.shapes.at(alight->shape_idx))){
                Vector3 p0, p1, p2;
                triangle_points(tri, p0, p1, p2);
                Real t_shadow = get_tri_intersect(p0, p1, p2, omeganot, pt, eps, uv);
                if (t_shadow > 0.0){
                    Vector3 hit_pt = pt + t_shadow*omeganot;
                    //std::cout << "asdfasdfasd" << std::endl;
                    d_squared = distance_squared(pt, hit_pt);
                    Vector3 norm =  shading_norm(tri, uv);
                    Real cosine = dot(omeganot,-norm);
                    pdf += (d_squared/(cosine*(0.5*length(cross(p2-p0,p1-p0)))))/scene.lights.size();
                } else {
                    continue;
                }
            } else {
                assert(false);
            }
        }
    }
    return pdf;
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
            Vector3 p0, p1, p2;
            triangle_points(tri, p0, p1, p2);
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
            return emission + 2.0*(kd + (1.0-kd)* pow((1.0 - dot(gn,omeganot)), 5))*mis_path_trace(scene, omeganot, pt, pcg_state, max_depth-1);
        } else if (spec==2) { // dielectric
            return emission + mis_path_trace(scene, omeganot, pt, pcg_state, max_depth-1);
        } else if (spec==3) { // refracting dielectric
            return emission + mis_path_trace(scene, omeganot, pt+(ray*eps), pcg_state, max_depth-1); // move by epsilon avoid moire patterns
        } else {
            Vector3 val;
            Real brdfpdf;
            Vector3 light_ray = light_sample(scene, ray, pt, eps, mat, mat_id, sn, gn, pcg_state);
            Real coinflip = next_pcg32_real<Real>(pcg_state);
            if (coinflip < 0.5){ // use light ray
                if (brdf_eval(scene, ray, light_ray, mat, mat_id, sn, gn, kd, pcg_state, val, brdfpdf)){
                    Real combinedpdf = 0.5*brdfpdf + 0.5*light_pdf(scene, light_ray, pt);
                    if (max(val) > Real(0) && combinedpdf > Real(0)){
                        return emission + val*mis_path_trace(scene, light_ray, pt, pcg_state, max_depth-1)/combinedpdf;
                    }
                }
            } else {
                if (brdf_eval(scene, ray, omeganot, mat, mat_id, sn, gn, kd, pcg_state, val, brdfpdf)){
                    Real combinedpdf = 0.5*brdfpdf + 0.5*light_pdf(scene, omeganot, pt);
                    if (max(val) > Real(0) && combinedpdf > Real(0)){
                        return emission + val*mis_path_trace(scene, omeganot, pt, pcg_state, max_depth-1)/combinedpdf;
                    }
                }
            }
            return emission;
        }
    } else {
        return scene.background_color;
    }
}