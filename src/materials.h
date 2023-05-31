#pragma once
#include "vector.h"
#include "image.h"
#include <variant>
#include <filesystem>

enum material_e {
    DiffuseType,
    MirrorType,
    DielectricType,
    PlasticType,
    PhongType,
    BlinnPhongType,
    BlinnPhongMicrofacetType
};

struct SolidTexture {
    Vector3 reflectance;
};
struct ImgTexture {
    ImgTexture(const fs::path &path, const Real us, const Real vs, const Real uo, const Real vo);
    int width, height;
    Real uoffset, voffset;
    Real uscale, vscale;
    Image3 data;
};

using Texture = std::variant<SolidTexture, ImgTexture>;

struct Material {
    Texture reflectance;
    Real ref_index; // ineta for dielectrics
    Real exponent; // outeta for dielectrics
    material_e material_type;
};

Vector3 get_texture_kd(const Texture &texture, const Vector2 uv);