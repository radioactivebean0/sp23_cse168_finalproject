#include "materials.h"
ImgTexture::ImgTexture(const fs::path &path, const Real us, const Real vs, const Real uo, const Real vo): 
    uoffset(uo),
    voffset(vo),
    uscale(us),
    vscale(vs)
{
    data = imread3(path);
    width = data.width;
    height = data.height;
}

Vector3 get_texture_kd(const Texture &texture, const Vector2 uv){
    if (auto *solid = std::get_if<SolidTexture>(&texture)){
        return solid->reflectance;
    } else if (auto *img = std::get_if<ImgTexture>(&texture)){
        // bilinear interpolation
        Real x = Real(img->width) * modulo(img->uscale * uv.x + img->uoffset, 1.0);
        Real y = Real(img->height) * modulo(img->vscale * uv.y + img->voffset, 1.0);
        int x1 = std::floor(x); int x2 = std::ceil(x);
        int y1 = std::floor(y); int y2 = std::ceil(y);

        if (x2>=img->width){
            x2 = 0;
        }
        if (y2 >= img->height){
            y2 = 0;
        }
        Vector3 q11 = img->data(x1,y1);
        Vector3 q21 = img->data(x2,y1);
        Vector3 q12 = img->data(x1,y2);
        Vector3 q22 = img->data(x2,y2);
        Vector3 fxy1 = ((1.0 - (x - x1))/1 * q11) +( (x - Real(x1))/1 * q21);
        Vector3 fxy2 = ((1.0 - (x - x1))/1 * q12) +( (x - Real(x1))/1 * q22);
        return ((1.0 - (y - y1))/(1) * fxy1 )+ ((y - Real(y1))/(1) * fxy2);
    } else {
        assert(false);
    }
}