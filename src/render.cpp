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
                //std::cout << x << "," << y << std::endl;
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
                //std::cout << x << "," << y << std::endl;
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

void ppm(Scene &scene, Image3 &img, int max_depth, const int num_tiles_x, const int num_tiles_y, const int tile_size, const int w, const int h, const int spp, const long photon_count, const Real alpha, const int passes, const Real default_radius){
    const Real real_spp = Real(spp);
    // STEP 1: Trace rays into the scene and get hit points
    PPMGrid ppm_pixels(w,h,spp);
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
    }, Vector2i(num_tiles_x, num_tiles_y));

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

    for (int pass = 0; pass < passes; pass++){
        // TODO
        // STEP 2: Trace photons into the scene
        for (int photon = 0; photon < photon_count; photon++){ // TODO: also parallelize
            // uniform sample and choose a light source

            // sample ray from the light source
            for (int depth = 0; depth < max_depth; depth++){
                // intersect the ray with scene

                // put photon into the kd tree

                // calculate bounce if any
            }
        }

        // TODO
        // STEP 3: Gather photons for the hit points

        // TODO
        // STEP 4: Adjust the radius of the visible points, discard all photons and repeat
    }
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
            max_depth,
            num_tiles_x,
            num_tiles_y,
            tile_size,
            w,
            h,
            spp,
            // FIXME use values that actually make sense
            // const long photon_count, const Real alpha, const int passes, const Real default_radius);
            10,
            10,
            2,
            0.1
        );
    } else {
        assert("unsupported render method");
    }
    std::cout << "Scene image render done. Took " << tick(timer) << " seconds." << std::endl;
    free(scene.cbvh);
    return img;
}