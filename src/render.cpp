#include "render.h"
#include "parse_scene.h"
#include "scene.h"
#include "timer.h"
#include "progressreporter.h"
#include "parallel.h"
#include "bvh.h"
#include "mis.h"
#include "next_event.h"
#include "path_trace.h"
#include "ppm.h"
#include "sampling.h"
#include "ortho_basis.h"
#include "3rdparty/nanoflann.hpp"

#include <stdlib.h>

void multiple_importance(Scene &scene, Image3 &img, int max_depth, const int num_tiles_x, const int num_tiles_y, const int tile_size, const int w, const int h, const int spp){
    const Real real_spp = Real(spp);
    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        pcg32_state pcg_state = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        const int x0 = tile[0] * tile_size;
        const int x1 = min(x0 + tile_size, w);
        const int y0 = tile[1] * tile_size;
        const int y1 = min(y0 + tile_size, h);
        int x, y, sample;
        Vector3 color_sum, ray;
        for (y = y0; y < y1; ++y){
            for (x = x0; x < x1; ++x){
                color_sum = Vector3{0.0,0.0,0.0};
                for (sample = 0; sample < spp; ++sample){
                    ray = get_ray(scene.camera, x, y, pcg_state);
                    color_sum += mis_path_trace(scene, ray, scene.camera.lookfrom, pcg_state, max_depth);
                }
                img(x,y) = color_sum/real_spp;
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
}


void next_event(Scene &scene, Image3 &img, const int num_tiles_x, const int num_tiles_y, const int tile_size, const int w, const int h, const int spp){
    const Real real_spp = Real(spp);
    const Real eps = 0.00000001;
    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        pcg32_state pcg_state = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        const int x0 = tile[0] * tile_size;
        const int x1 = min(x0 + tile_size, w);
        const int y0 = tile[1] * tile_size;
        const int y1 = min(y0 + tile_size, h);
        int x, y, sample;
        Real t;
        Vector3 color_sum, ray, color, pt;
        Vector2 uv;
        Shape *hs;
        for (y = y0; y < y1; ++y){
            for (x = x0; x < x1; ++x){
                color_sum = Vector3{0.0,0.0,0.0};
                for (sample = 0; sample < spp; ++sample){
                    ray = get_ray(scene.camera, x, y, pcg_state);
                    t = -1.0;
                    if (hit_cbvh(scene.cbvh, ray, scene.camera.lookfrom, eps, eps, infinity<Real>(), &hs, t, uv)){ // was a hit
                        pt = scene.camera.lookfrom + t * ray;
                        color = radiance_v2(scene, ray, pt, eps, hs, uv, pcg_state);
                    } else {
                        color = scene.background_color;
                    }
                    color_sum += color;
                }
                img(x,y) = color_sum/real_spp;
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
}

void path(Scene &scene, Image3 &img, int max_depth, const int num_tiles_x, const int num_tiles_y, const int tile_size, const int w, const int h, const int spp){
    const Real real_spp = Real(spp);
    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        pcg32_state pcg_state = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        const int x0 = tile[0] * tile_size;
        const int x1 = min(x0 + tile_size, w);
        const int y0 = tile[1] * tile_size;
        const int y1 = min(y0 + tile_size, h);
        int x, y, sample;
        Vector3 color_sum, ray;
        for (y = y0; y < y1; ++y){
            for (x = x0; x < x1; ++x){
                color_sum = Vector3{0.0,0.0,0.0};
                for (sample = 0; sample < spp; ++sample){
                    ray = get_ray(scene.camera, x, y, pcg_state);
                    color_sum += path_trace(scene, ray, scene.camera.lookfrom, pcg_state, max_depth);
                }
                img(x,y) = color_sum/real_spp;
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
}

// TODO clean this up and move elsewhere
struct PointCloud {
    struct Point {
        Real x, y, z;
    };

    using coord_t = Real;  //!< The type of each coordinate

    std::vector<Point> pts;

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline Real kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0)
            return pts[idx].x;
        else if (dim == 1)
            return pts[idx].y;
        else
            return pts[idx].z;
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template<class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
        return false;
    }
};

// todo, move this out
Vector3 photon_bounce(Scene &scene, Shape *hs, const Vector3 &omeganot, Vector3 &hit_pt, const Vector2 &uv, pcg32_state &pcg_state){
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

Vector3 photon_color(Scene &scene, Shape *hs, const Vector3 &omegai, const Vector3 &hit_pt, const Vector2 &uv) {
    int mat_id = shape_matid(hs);
    material_e mat = scene.materials.at(mat_id).material_type;
    Vector3 kd = get_texture_kd(scene.materials.at(mat_id).reflectance, uv);
    if (mat==material_e::MirrorType){ // purely specular
        Vector3 gn = shape_geo_norm(hs, hit_pt, uv);
        return 2.0*(kd + (1.0-kd)* pow((1.0 - dot(gn,omegai)), 5));
    } else if (mat == material_e::DielectricType) { // dielectric
        return Vector3{1.0,1.0,1.0};
    } else {        // scattering, cosine hemisphere sampling and diffuse
        Vector3 sn = shape_shade_norm(hs, hit_pt, uv);
        Real nwo = dot(sn,omegai);
        if (nwo <= 0.0){
            return Vector3{0.0,0.0,0.0};
        }
        Real pdf = dot(sn,omegai)*c_INVPI;
        return (kd*nwo*c_INVPI)/pdf;
    }
}

void ppm(Scene &scene, Image3 &img, int max_depth, const int num_tiles_x, const int num_tiles_y, const int tile_size, const int w, const int h, const int spp, const long photon_count, const Real alpha, const int passes, const Real default_radius){
    // STEP 1: Trace rays into the scene and get hit points
    PPMGrid ppm_pixels(w,h,spp);
    ProgressReporter reporter(num_tiles_x * num_tiles_y);

    parallel_for([&](const Vector2i &tile) {
        pcg32_state pcg_state = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        const int x0 = tile[0] * tile_size;
        const int x1 = min(x0 + tile_size, w);
        const int y0 = tile[1] * tile_size;
        const int y1 = min(y0 + tile_size, h);
        int x, y, sample;
        Vector3 ray;
        for (y = y0; y < y1; ++y){
            for (x = x0; x < x1; ++x){
                for (sample = 0; sample < spp; ++sample){
                    ray = get_ray(scene.camera, x, y, pcg_state);
                    ppm_pixels(x,y).push_back(generate_visible_point(scene, ray, pcg_state, max_depth, default_radius));
                }
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    // TODO remove this kd tree test code

    {
        auto n = 100000;
        using std::cout;
        using std::endl;

        PointCloud cloud;
        cloud.pts.resize(n);
        for (size_t i = 0; i < n; i++) {
            cloud.pts[i].x = 10 * (rand() % 1000) / Real(1000);
            cloud.pts[i].y = 10 * (rand() % 1000) / Real(1000);
            cloud.pts[i].z = 10 * (rand() % 1000) / Real(1000);
        }

        // construct a kd-tree index:
        using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<Real, PointCloud>,
            PointCloud, 3 /* dim */
        >;

        my_kd_tree_t index(3 /*dim*/, cloud, {10 /* max leaf */});

        const Real query_pt[3] = {0.5, 0.5, 0.5};
        const Real search_radius = static_cast<Real>(0.1);
        std::vector<nanoflann::ResultItem<uint32_t, Real>> ret_matches;

        // nanoflanSearchParamsameters params;
        // params.sorted = false;

        const size_t nMatches =
            index.radiusSearch(&query_pt[0], search_radius, ret_matches);

        cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches
             << " matches\n";
        for (size_t i = 0; i < nMatches; i++)
            cout << "idx[" << i << "]=" << ret_matches[i].first << " dist[" << i
                 << "]=" << ret_matches[i].second << endl;
        cout << "\n";
    }
    pcg32_state pcg_state = init_pcg32();
    const Real eps = 0.0000001;

    for (int pass = 0; pass < passes; pass++){
        // TODO
        // STEP 2: Trace photons into the scene, keep going until photon count reached
        long photon = 0;
        int lit;
        Vector3 photon_ori, photon_dir;
        PointCloud cloud;
        cloud.pts.resize(photon_count);
        ProgressReporter reporter(photon_count);

        while(photon < photon_count){ // TODO: also parallelize, also probably move out of here
            // uniform sample and choose a light source
            lit = next_pcg32_real<Real>(pcg_state) * scene.lights.size();
            // sample ray from light source
            Vector3 luminance;
            if (auto *alight = std::get_if<AreaLight>(&scene.lights.at(lit))){
                Vector3 light_norm;
                photon_ori = sample_shape_point(&(scene.shapes.at(alight->shape_idx)), pcg_state, light_norm);
                photon_dir = ortho_basis(rand_cos(pcg_state), light_norm);
                luminance = dot(photon_dir, light_norm) * alight->radiance;
            } else if (auto *light = std::get_if<PointLight>(&scene.lights.at(lit))){
                photon_ori = light->position;
                photon_dir = rand_uniform_sphere(pcg_state);
                luminance = light->intensity;
            } else {
                assert(false);
            }
            Vector2 uv;
            Shape *hs;
            // intersect the ray with scene
            for (int depth = 0; depth < max_depth; depth++){
                Real t = -1.0;

                if (hit_cbvh(scene.cbvh, photon_dir, photon_ori, eps, eps, infinity<Real>(), &hs, t, uv)){ // was a hit, deposit photon
                    // put photon into the kd tree
                    Vector3 hit_pt = photon_ori + t*photon_dir;
                    cloud.pts[photon].x = hit_pt.x;
                    cloud.pts[photon].y = hit_pt.y;
                    cloud.pts[photon].z = hit_pt.z;
                    photon++;
                    reporter.update(1);
                    // russian roulette, reference: https://www.pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/Stochastic_Progressive_Photon_Mapping#AccumulatingVisiblePoints
                    if (luminance.y < 0.25) {
                        if (next_pcg32_real<Real>(pcg_state) > luminance.y) {
                            break;
                        }
                        luminance /= luminance.y;
                    }
                    // calculate bounce
                    Vector3 hit_pt_cpy = hit_pt;                    // copy hit point since possible it got moved
                    photon_dir = photon_bounce(scene, hs, photon_dir, hit_pt, uv, pcg_state); 
                    photon_ori = hit_pt;
                    luminance *= photon_color(scene, hs, photon_dir, hit_pt_cpy, uv);
                } else { // ray goes away
                    break;
                }
            }
        }
        using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<Real, PointCloud>,
            PointCloud, 3 /* dim */
        >;

        my_kd_tree_t photon_tree(3 /*dim*/, cloud, {10 /* max leaf */});
        reporter.done();
        // TODO
        // STEP 3: Gather photons for the hit points

        // TODO
        // STEP 4: Adjust the radius of the visible points, discard all photons and repeat
    }
    return;
}

Image3 render_img(const std::vector<std::string> &params) {
    
    if (params.size() < 1) {
        return Image3(0, 0);
    }

    int max_depth = 50;
    std::string filename;
    std::string renderer = "mis";
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-max_depth") {
            max_depth = std::stoi(params[++i]);
        } else if (filename.empty()) {
            filename = params[i];
        } else if (params[i] == "-renderer") {
            renderer = params[++i];
        }
    }

    Timer timer;
    ParsedScene pscene = parse_scene(params[0]);
    std::cout << "Scene parsing done. Took " << tick(timer) << " seconds." << std::endl;
    tick(timer);
    Scene scene = Scene(pscene);
    std::cout << "Scene data construction done. Took " << tick(timer) << " seconds." << std::endl;
    tick(timer);
    int bvhsize = 0;
    scene.bvh = build_bvh(scene.shapes, 2, &bvhsize);
    scene.cbvh = (compact_AABB*) aligned_alloc(64, sizeof(compact_AABB)*bvhsize);
    int offset = 0;
    flatten_bvh(scene.bvh, scene.cbvh, &offset);
    std::cout << "BVH construction done. Took " << tick(timer) << " seconds." << std::endl;
    tick(timer);
    const int spp = scene.samples_per_pixel;
    Image3 img(scene.width, scene.height);
    const int w = img.width;
    const int h = img.height;
    constexpr int tile_size = 16;
    const int num_tiles_x = (w + tile_size - 1) / tile_size;
    const int num_tiles_y = (h + tile_size - 1) / tile_size;
    if (renderer=="mis") {
        multiple_importance(scene, img, max_depth, num_tiles_x, num_tiles_y, tile_size, w, h, spp);
    } else if (renderer == "nee"){
        next_event(scene, img, num_tiles_x, num_tiles_y, tile_size, w, h, spp);
    } else if (renderer == "path"){
        path(scene, img, max_depth, num_tiles_x, num_tiles_y, tile_size, w, h, spp);
    } else if (renderer == "ppm") {
        ppm(
            scene,
            img,
            10,
            num_tiles_x,
            num_tiles_y,
            tile_size,
            w,
            h,
            spp,
            // FIXME use values that actually make sense
            // const long photon_count, const Real alpha, const int passes, const Real default_radius);
            50000,
            0.7, // alpha value from the paper
            5,
            0.1
        );
    } else {
        assert("unsupported render method");
    }
    std::cout << "Scene image render done. Took " << tick(timer) << " seconds." << std::endl;
    free(scene.cbvh);
    return img;
}