#include "camera.h"

Camera from_parsed_camera(const ParsedCamera &camera){
    Real theta = radians(camera.vfov);
    Real h = tan(theta/2);
    Real ratio = Real(camera.width)/camera.height;

    Vector3 w = normalize(camera.lookfrom - camera.lookat);
    Vector3 u = normalize(cross(camera.up, w));
    Vector3 v = cross(w, u);

    Vector3 vertical = 2.0 * h * v;
    Vector3 horizontal = ratio * 2.0 * h * u;
    Vector3 llc = camera.lookfrom - horizontal/Real(2.0) + vertical/Real(2.0) - w;

    return Camera{
        .lookfrom= camera.lookfrom,
        .lookat= camera.lookat,
        .up= camera.up,
        .vfov= camera.vfov,
        .llc=llc,
        .horizontal=horizontal,
        .vertical=vertical,
        .width= camera.width,
        .height= camera.height,
    };
}

Vector3 get_ray(Camera &camera, int x, int y, pcg32_state &pcg_state){

    Real ux;
    Real vy;
    ux = ((x + next_pcg32_real<double>(pcg_state)) / camera.width);
    vy = ((y + next_pcg32_real<double>(pcg_state)) / camera.height);

    return normalize(camera.llc + ux*camera.horizontal - vy*camera.vertical - camera.lookfrom);
}