#include "path_trace.h"
#include "bvh.h"
#include "sampling.h"

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