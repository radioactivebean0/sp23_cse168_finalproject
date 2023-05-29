#include "render.h"
#include "parse_scene.h"
#include "scene.h"
#include "timer.h"
#include "progressreporter.h"
#include "parallel.h"
#include "bvh.h"
#include <stdlib.h>

Image3 render_img(const std::vector<std::string> &params) {
    
    if (params.size() < 1) {
        return Image3(0, 0);
    }

    int max_depth = 50;
    std::string filename;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-max_depth") {
            max_depth = std::stoi(params[++i]);
        } else if (filename.empty()) {
            filename = params[i];
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
    const Real real_spp = Real(spp);
    Image3 img(scene.width, scene.height);
    const int w = img.width;
    const int h = img.height;
    constexpr int tile_size = 16;
    const int num_tiles_x = (w + tile_size - 1) / tile_size;
    const int num_tiles_y = (h + tile_size - 1) / tile_size;
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
                    color_sum += mis_path_trace(scene, ray, scene.camera.lookfrom, pcg_state, max_depth);
                }
                img(x,y) = color_sum/real_spp;
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    std::cout << "Scene image render done. Took " << tick(timer) << " seconds." << std::endl;
    free(scene.cbvh);
    return img;
}