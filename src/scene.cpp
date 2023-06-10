#include "scene.h"
#include "bvh.h"
#include <unordered_map>
Scene::Scene(const ParsedScene &scene) :
        camera(from_parsed_camera(scene.camera)),
        width(scene.camera.width),
        height(scene.camera.height),
        background_color(scene.background_color),
        samples_per_pixel(scene.samples_per_pixel){
    // Extract triangle meshes from the parsed scene.
    int tri_mesh_count = 0;
    for (const ParsedShape &parsed_shape : scene.shapes) {
        if (std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            tri_mesh_count++;
        }
    }
    meshes.resize(tri_mesh_count);
    // Extract the shapes
    std::unordered_map<int, std::vector<int>> arealight_map;
    tri_mesh_count = 0;
    for (int i = 0; i < (int)scene.shapes.size(); i++) {
        const ParsedShape &parsed_shape = scene.shapes[i];
        if (auto *sph = std::get_if<ParsedSphere>(&parsed_shape)) {
            shapes.push_back(Sphere {
                .center = sph->position,
                .radius = sph->radius,
                .material_id = sph->material_id,
                .area_light_id = sph->area_light_id
            });
            if (std::get<Sphere>(shapes.back()).area_light_id != -1){
                arealight_map.emplace(std::make_pair(i,std::vector<int>{(int)(shapes.size()-1)}));
            }
        } else if (auto *parsed_mesh = std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            meshes[tri_mesh_count] = TriangleMesh {
                .material_id = parsed_mesh->material_id,
                .area_light_id = parsed_mesh->area_light_id,
                .positions = parsed_mesh->positions,
                .indices = parsed_mesh->indices,
                .normals = parsed_mesh->normals,
                .uvs = parsed_mesh->uvs,
            };
            // Extract all the individual triangles
            if (meshes[tri_mesh_count].area_light_id!=-1){
                arealight_map.emplace(std::make_pair(i, std::vector<int>{}));
            }
            for (int face_index = 0; face_index < (int)parsed_mesh->indices.size(); face_index++) {
                shapes.push_back(Triangle{.face_index = face_index, .mesh = &meshes[tri_mesh_count]});
                if (meshes[tri_mesh_count].area_light_id!=-1){
                    arealight_map[i].push_back(shapes.size()-1);
                }
            }
            tri_mesh_count++;
        } else {
            assert(false);
        }
    }
    // Copy the materials
    for (const ParsedMaterial &parsed_mat : scene.materials) {
        if (auto *diffuse = std::get_if<ParsedDiffuse>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&diffuse->reflectance)){
                materials.push_back(Material {
                    .reflectance = SolidTexture { .reflectance = *tex },
                    .material_type = material_e::DiffuseType,
                });
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&diffuse->reflectance)){
                materials.push_back(Material {
                    .reflectance = ImgTexture(
                        tex->filename,
                        tex->uscale,
                        tex->vscale,
                        tex->uoffset,
                        tex->voffset
                    ),
                    .material_type = material_e::DiffuseType,
                });
            }
        } else if (auto *mirror = std::get_if<ParsedMirror>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&mirror->reflectance)){
                materials.push_back(Material {
                    .reflectance = SolidTexture { .reflectance = *tex },
                    .material_type = material_e::MirrorType,
                });
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&mirror->reflectance)){
                materials.push_back(Material {
                    .reflectance = ImgTexture(
                        tex->filename,
                        tex->uscale,
                        tex->vscale,
                        tex->uoffset,
                        tex->voffset
                    ),
                    .material_type = material_e::MirrorType,
                });
            }
        } else if (auto *dielec = std::get_if<ParsedDielectric>(&parsed_mat)){
            materials.push_back(Material {
                .reflectance = SolidTexture { .reflectance = Vector3 { 1.0, 1.0, 1.0 } },
                .ref_index = dielec->inEta,
                .exponent = dielec->outEta,
                .material_type = material_e::DielectricType
            });
        } else if (auto *plastic = std::get_if<ParsedPlastic>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&plastic->reflectance)){
                materials.push_back(Material {
                    .reflectance = SolidTexture { .reflectance = *tex },
                    .ref_index = plastic->eta,
                    .material_type = material_e::PlasticType,
                });
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&plastic->reflectance)){
                materials.push_back(Material {
                    .reflectance = ImgTexture(
                        tex->filename,
                        tex->uscale,
                        tex->vscale,
                        tex->uoffset,
                        tex->voffset
                    ),
                    .ref_index = plastic->eta,
                    .material_type = material_e::PlasticType,
                });
            }
        } else if (auto *phong = std::get_if<ParsedPhong>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&phong->reflectance)){
                materials.push_back(Material {
                    .reflectance = SolidTexture { .reflectance = *tex },
                    .exponent = phong->exponent,
                    .material_type = material_e::PhongType,
                });
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&phong->reflectance)){
                materials.push_back(Material {
                    .reflectance = ImgTexture(
                        tex->filename,
                        tex->uscale,
                        tex->vscale,
                        tex->uoffset,
                        tex->voffset
                    ),
                    .exponent = phong->exponent,
                    .material_type = material_e::PhongType,
                });
            }
        } else if (auto *blinnphong = std::get_if<ParsedBlinnPhong>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&blinnphong->reflectance)){
                materials.push_back(Material {
                    .reflectance = SolidTexture { .reflectance = *tex },
                    .exponent = blinnphong->exponent,
                    .material_type = material_e::BlinnPhongType,
                });
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&blinnphong->reflectance)){
                materials.push_back(Material {
                    .reflectance = ImgTexture(
                        tex->filename,
                        tex->uscale,
                        tex->vscale,
                        tex->uoffset,
                        tex->voffset
                    ),
                    .exponent=blinnphong->exponent,
                    .material_type = material_e::BlinnPhongType,
                });
            }
        } else if (auto *blinnphongmicrofacet = std::get_if<ParsedBlinnPhongMicrofacet>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&blinnphongmicrofacet->reflectance)){
                materials.push_back(Material {
                    .reflectance = SolidTexture { .reflectance = *tex },
                    .exponent=blinnphongmicrofacet->exponent,
                    .material_type = material_e::BlinnPhongMicrofacetType,
                });
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&blinnphongmicrofacet->reflectance)){
                materials.push_back(Material {
                    .reflectance = ImgTexture(
                        tex->filename,
                        tex->uscale,
                        tex->vscale,
                        tex->uoffset,
                        tex->voffset
                    ),
                    .exponent=blinnphongmicrofacet->exponent,
                    .material_type = material_e::BlinnPhongMicrofacetType,
                });
            }
        } else {
            assert(false);
        }
    }
    // Copy the lights
    for (const ParsedLight &parsed_light : scene.lights) {
        // We assume all lights are point lights for now.
        if (auto *point_light = std::get_if<ParsedPointLight>(&parsed_light)){
            lights.push_back(PointLight{.position= point_light->position , .intensity= point_light->intensity});
        } else if (auto *area_light = std::get_if<ParsedDiffuseAreaLight>(&parsed_light)){
            for (int s: arealight_map.at(area_light->shape_id)){
                lights.push_back(AreaLight {
                    .shape_idx = s,
                    .radiance = area_light->radiance,
                });
            }
        } else {
            assert(false);
        }
    }
}
