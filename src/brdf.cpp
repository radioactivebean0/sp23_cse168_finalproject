#include "brdf.h"
#include "sampling.h"
#include "ortho_basis.h"

// calulates brdf bounce direction, if dielectric will move the hitpoint to avoid moire patterns
Vector3 brdf_bounce(const Scene &scene, Shape *hs, const Vector3 &omeganot, Vector3 &hit_pt, const Vector2 &uv, pcg32_state &pcg_state){
    int mat_id = shape_matid(hs);
    material_e mat = scene.materials.at(mat_id).material_type;
    Vector3 n = shape_shade_norm(hs, hit_pt, uv);
    if (mat == material_e::MirrorType){ // mirrors
        return omeganot - 2.0*dot(omeganot, n)*n;
    } else if (mat == material_e::DielectricType){
        const Real inIR = scene.materials.at(mat_id).ref_index;
        const Real outIR = scene.materials.at(mat_id).exponent;
        Real cos_theta = fmin(dot(-omeganot, n), 1.0);
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
            direction = omeganot - 2.0*dot(omeganot,n)*n;
        } else { // refract
            Vector3 r_out_perp = ref_ratio * (omeganot + cos_theta*n);
            Vector3 r_out_parallel = -sqrt(fabs(1.0 - length_squared(r_out_perp))) * n;
            direction = r_out_perp + r_out_parallel;
            hit_pt = hit_pt + 0.00001*omeganot; // adjust to avoid moire
        }
        return direction;
    } else {        // scattering, cosine hemisphere sampling and diffuse
        Vector3 scatter = rand_cos(pcg_state);
        Vector3 bounce = ortho_basis(scatter, n);
        return bounce;
    }
}

// evaluates brdf for omega not, normalizes with pdf
Vector3 brdf_color(const Scene &scene, Shape *hs, const Vector3 &omeganot, const Vector3 &hit_pt, const Vector2 &uv) {
    int mat_id = shape_matid(hs);
    material_e mat = scene.materials.at(mat_id).material_type;
    Vector3 kd = get_texture_kd(scene.materials.at(mat_id).reflectance, uv);
    if (mat==material_e::MirrorType){ // purely specular
        Vector3 gn = shape_geo_norm(hs, hit_pt, uv);
        return (kd + (1.0-kd)* pow((1.0 - dot(gn,omeganot)), 5));
    } else if (mat == material_e::DielectricType) { // dielectric
        return Vector3{1.0,1.0,1.0};
    } else {        // scattering, cosine hemisphere sampling and diffuse
        Vector3 sn = shape_shade_norm(hs, hit_pt, uv);
        Real nwo = dot(sn,omeganot);
        if (nwo <= 0.0){
            return Vector3{0.0,0.0,0.0};
        }
        Real pdf = dot(sn,omeganot)*c_INVPI;
        if (pdf == 0.0){
            return Vector3{0.0,0.0,0.0};
        }
        return (kd*nwo*c_INVPI)/pdf;
    }
}