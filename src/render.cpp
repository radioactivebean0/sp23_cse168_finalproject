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
#include "brdf.h"
#include "3rdparty/nanoflann.hpp"
#include <unordered_map>

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

// evaluate tau which is the flux contribution of a photon that is accumulated on a hit point in the paper
Vector3 eval_tau(const PPMHitPoint &vp, const Vector3 &omegai, const Vector3 &phi_phot){
    if (vp.mat == material_e::DiffuseType){
        Real nwo = dot(vp.normal, -omegai);
        if (nwo <= 0.0){
            return Vector3{0.0,0.0,0.0}; // hitting wrong side of surface?
        }
        return phi_phot*vp.beta*nwo*c_INVPI;
    } else { // mirror, dielectric dont store flux, other mats unsupported will do the same for now
        return Vector3{0.0,0.0,0.0};
    }
}

void ppm(
    Scene &scene,
    Image3 &img,
    int max_depth,
    const int num_tiles_x,
    const int num_tiles_y,
    const int tile_size,
    const int w,
    const int h,
    const int spp,
    const long photon_count,
    const Real alpha,
    const int passes,
    const Real default_radius
) {
    // STEP 1: Trace rays into the scene and get hit points
    PPMGrid ppm_pixels(w,h,spp);
    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    std::cout << "tracing visible points" << std::endl;
    // for each pixel, generate sample visible points
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

    pcg32_state pcg_state = init_pcg32();
    const Real eps = 0.0000001;

    for (int pass = 0; pass < passes; pass++){
        // TODO
        // STEP 2: Trace photons into the scene, keep going until photon count reached
        std::cout << "accumulating photons pass " << pass << "/" << passes << std::endl;
        long photon = 0;
        int lit;
        Vector3 photon_ori, photon_dir;
        PointCloud cloud;
        cloud.pts.resize(photon_count);
        std::cout << "dispersing photons" << std::endl;
        ProgressReporter reporter(photon_count);
        std::unordered_map<Vector3, std::pair<Vector3,Vector3> > flux_map; // store location to (omegai, flux) map
        while(photon < photon_count){ // TODO: also parallelize, also probably move out of here
            // uniform sample and choose a light source
            lit = next_pcg32_real<Real>(pcg_state) * scene.lights.size();
            // sample ray from light source
            Vector3 luminance;
            Real light_pdf;
            if (auto *alight = std::get_if<AreaLight>(&scene.lights.at(lit))){
                Vector3 light_norm;
                photon_ori = sample_shape_point(&(scene.shapes.at(alight->shape_idx)), pcg_state, light_norm, light_pdf);
                photon_dir = ortho_basis(rand_cos(pcg_state), light_norm);
                luminance = (dot(photon_dir, light_norm) * alight->radiance)/(light_pdf*(1.0/scene.lights.size()));
            } else if (auto *light = std::get_if<PointLight>(&scene.lights.at(lit))){
                photon_ori = light->position;
                photon_dir = rand_uniform_sphere(pcg_state);
                luminance = light->intensity/((1.0/scene.lights.size())*c_INVFOURPI); // normalize by chance of picking light and sphere pdf of picking
            } else {
                assert(false);
            }
            Vector2 uv;
            Shape *hs;
            // intersect the ray with scene
            for (int depth = 0; depth < max_depth; depth++){
                Real t = -1.0;
                if (hit_cbvh(scene.cbvh, photon_dir, photon_ori, eps, eps, infinity<Real>(), &hs, t, uv)) {
                    // was a hit, deposit photon into the kd tree
                    Vector3 hit_pt = photon_ori + t*photon_dir;
                    cloud.pts[photon].x = hit_pt.x;
                    cloud.pts[photon].y = hit_pt.y;
                    cloud.pts[photon].z = hit_pt.z;
                    if (flux_map.find(hit_pt) != flux_map.end() ){
                        std::cerr << "Error: hash collision" << std::endl;
                        exit(1);
                    }
                    flux_map.insert(std::make_pair(hit_pt, std::make_pair(photon_dir, luminance)));
                    photon++;
                    reporter.update(1);
                    if (photon >= photon_count) {
                        break;
                    }

                    Vector3 hit_pt_cpy = hit_pt;
                    photon_dir = brdf_bounce(scene, hs, photon_dir, hit_pt, uv, pcg_state); 

                    Vector3 fr = brdf_color(scene, hs, photon_dir, hit_pt_cpy, uv);
                    if (fr==Vector3{0.0,0.0,0.0}){
                        break;
                    }
                    Vector3 newluminance = luminance * fr;

                    // russian roulette, reference: https://www.pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/Stochastic_Progressive_Photon_Mapping#AccumulatingPhotonContributions
                    Real q = max(Real(0.0), 1.0 - newluminance.y / luminance.y);
                    if (next_pcg32_real<Real>(pcg_state) < q) {
                        break;
                    }
                    luminance = newluminance / (1.0 - q);
                    // calculate bounce
                    photon_ori = hit_pt;
                } else { // ray goes away
                    break;
                }
            }
        }
        if (photon_count != photon) {
            std::cerr << "Error: unexpected photon count for n_emitted" << std::endl;
            exit(1);
        }
        using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<Real, PointCloud>,
            PointCloud, 3 /* dim */
        >;

        kd_tree_t photon_tree(3 /*dim*/, cloud, {10 /* max leaf */});
        reporter.done();
        std::cout << "calculating photon contributions" << std::endl;
        ProgressReporter reporter2(num_tiles_x*num_tiles_y);
        parallel_for([&](const Vector2i &tile) {
            const int x0 = tile[0] * tile_size;
            const int x1 = min(x0 + tile_size, w);
            const int y0 = tile[1] * tile_size;
            const int y1 = min(y0 + tile_size, h);
            int x, y;
            Vector3 acc_color;
            std::vector<PPMHitPoint> sample;
            std::vector<nanoflann::ResultItem<uint32_t, Real>> points;

            for (y = y0; y < y1; ++y){
                for (x = x0; x < x1; ++x){
                    //Real samples = Real(sample.size())/spp;
                    acc_color = Vector3{0.0,0.0,0.0};
                    for (auto& point: ppm_pixels(x,y)){
                        // STEP 3: Gather photons for the hit points
                        // TODO check if this makes sense
                        if (point.r < 0.0){
                            continue;
                        }
                        Real query[3];
                        for (size_t i = 0; i < 3; ++i) {
                            query[i] = point.position[i];
                        }
                        // TODO find way to not store into points
                        const size_t m = photon_tree.radiusSearch(&query[0], point.r, points);
                        if (m == 0){ // nothing in the radius, estimate a radius by getting 20 nearest neighbors
                            size_t find_nearest = 20;
                            std::vector<uint32_t> ret_index(find_nearest);
                            std::vector<Real> out_dist_sqr(find_nearest);

                            find_nearest = photon_tree.knnSearch(
                                &query[0], find_nearest, &ret_index[0], &out_dist_sqr[0]);

                            // In case of less points in the tree than requested:
                            ret_index.resize(find_nearest);
                            out_dist_sqr.resize(find_nearest);
                            point.r = sqrt(*max_element(out_dist_sqr.begin(), out_dist_sqr.end()));
                            continue;
                        }
                        // STEP 4: Adjust the radius of the visible points
                        // TODO check if this makes sense

                        point.r *= sqrt((point.n + alpha * m) / (point.n + m));
                        // TODO check type conversion

                        // STEP 5: accumulate flux
                        Vector3 tau_m = Vector3{0.0,0.0,0.0};
                        for (size_t i = 0; i < m; i++){
                            PointCloud::Point p_loc = cloud.pts[points[i].first];
                            std::pair<Vector3, Vector3> f_contrib = flux_map.at(Vector3{p_loc.x, p_loc.y, p_loc.z});
                            tau_m += eval_tau(point, f_contrib.first, f_contrib.second);
                            // std::cout << tau_m << std::endl;
                        }
                        // bool debug = rand() % 100000 == 0;
                        bool debug = false;
                        if (debug) {
                            std::cout << "[DEBUG] --- new iteration --- " << std::endl;
                            std::cout << "[DEBUG] original point.tau = " << point.tau << ", tau_m = " << tau_m << ", point.n = " << point.n << ", alpha = " << alpha << ", m = " << m << std::endl;
                            std::cout << "[DEBUG] so it's " << (point.tau + tau_m) << " / " << ((point.n + alpha * m) / (point.n + m)) << std::endl;
                        }
                        point.tau = (point.tau + tau_m) * ((point.n + alpha * m) / (point.n + m));
                        if (debug) {
                            std::cout << "[DEBUG] point.tau = " << point.tau << ", m = " << m << ", alpha = " << alpha << std::endl;
                            std::cout << "[DEBUG] and original point.n = " << point.n << std::endl;
                        }
                        point.n += m*alpha; // unsure about this, probably should be adding alpha*m photons since that was the point of radius reduction may run into weird rounding here
                        if (debug) {
                            std::cout << "[DEBUG] and new point.n = " << point.n << std::endl;
                        }
                        // floats? on the scale of 10
                        // if (point.n > 0 && rand() % 100000 == 0) {
                        //     std::cout << "[DEBUG] point.n = " << point.n << ", m = " << m << ", alpha = " << alpha << std::endl;
                        // }
                    }
                }
            }
            reporter2.update(1);
        }, Vector2i(num_tiles_x, num_tiles_y));
        reporter2.done();
    }
    // render img using photons
    std::cout << "rendering final img" << std::endl;
    const long n_emitted = photon_count;
    parallel_for([&](const Vector2i &tile) {
        // DEBUG
        auto tile_size = 10000000;
        const int x0 = tile[0] * tile_size;
        const int x1 = min(x0 + tile_size, w);
        const int y0 = tile[1] * tile_size;
        const int y1 = min(y0 + tile_size, h);
        int x, y;
        Vector3 acc_color;
        std::vector<PPMHitPoint> sample;
        for (y = y0; y < y1; ++y){
            for (x = x0; x < x1; ++x){
                sample = ppm_pixels(x,y);
                //Real samples = Real(sample.size())/spp;
                acc_color = Vector3{0.0,0.0,0.0};
                for (auto point: sample){
                    if (point.r < 0.0){
                        //acc_color += Vector3{1.0,1.0,1.0};
                        // acc_color += Vector3{Real(point.n),Real(point.n),Real(point.n)};
                        acc_color += scene.background_color;
                    }else {
                        acc_color += point.emission + (1.0/(c_PI*point.r*point.r*n_emitted))*(point.tau); //1000000 * (1.0/(c_PI*point.r*point.r*n_emitted))*(point.tau);
                        // acc_color += 100000 * (1.0/(c_PI*point.r*point.r* point.n))*(point.tau);
                        //std::cout << "[DEBUG] acc_color = " << acc_color << ", point.r = " << point.r << ", n_mitted = " << n_emitted << ", point.tau = " << point.tau << std::endl;
                        //std::cout << "[DEBUG] more info: point.n = " << point.n << std::endl;
                    }

                }
                img(x, y) = acc_color/Real(passes);/// Real(spp); ///Real(spp); ///((Real)sample.size());
                // auto point = img(x, y);
                // auto max = std::max({ point.x, point.y, point.z });
                // if (max > 1.0) {
                //     img(x, y).x = 1.0;
                //     img(x, y).y = 1.0;
                //     img(x, y).z = 1.0;
                //     // img(x, y) /= max;
                // }
                // std::cout << "[DEBUG] img(x, y) = " << img(x, y) << ", max = " << max << std::endl;
            }
        }
        reporter.update(1);
    // }, Vector2i(num_tiles_x, num_tiles_y));
    }, Vector2i(1, 1));
    reporter.done();
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
            15,
            num_tiles_x,
            num_tiles_y,
            tile_size,
            w,
            h,
            64, // hard code for now
            // FIXME use values that actually make sense
            // const long photon_count, const Real alpha, const int passes, const Real default_radius);
            200000,
            0.7, // alpha value from the paper
            20,
            20.0
        );
    } else {
        assert("unsupported render method");
    }
    std::cout << "Scene image render done. Took " << tick(timer) << " seconds." << std::endl;
    free(scene.cbvh);
    return img;
}