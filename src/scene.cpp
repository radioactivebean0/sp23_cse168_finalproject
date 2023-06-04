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
            shapes.push_back( Sphere{
                .center = sph->position,
                .material_id = sph->material_id,
                .radius = sph->radius,
                .area_light_id = sph->area_light_id
            });
            if (std::get<Sphere>(shapes.back()).area_light_id != -1){
                arealight_map.emplace(std::make_pair(i,std::vector<int>{(int)(shapes.size()-1)}));
            }
        } else if (auto *parsed_mesh = std::get_if<ParsedTriangleMesh>(&parsed_shape)) {
            meshes[tri_mesh_count] = TriangleMesh{.material_id= parsed_mesh->material_id,
                                                  .positions= parsed_mesh->positions,
                                                  .indices= parsed_mesh->indices, 
                                                  .normals = parsed_mesh->normals,
                                                  .uvs = parsed_mesh->uvs,
                                                  .area_light_id= parsed_mesh->area_light_id};
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
                materials.push_back(Material{.material_type = material_e::DiffuseType, .reflectance = SolidTexture{.reflectance = *tex}});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&diffuse->reflectance)){
                materials.push_back(Material{.material_type = material_e::DiffuseType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset)});
            }
        } else if (auto *mirror = std::get_if<ParsedMirror>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&mirror->reflectance)){
                materials.push_back(Material{.material_type = material_e::MirrorType, .reflectance = SolidTexture{.reflectance = *tex}});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&mirror->reflectance)){
                materials.push_back(Material{.material_type = material_e::MirrorType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset)});
            }
        } else if (auto *dielec = std::get_if<ParsedDielectric>(&parsed_mat)){
            materials.push_back(Material{.material_type = material_e::DielectricType, .reflectance = SolidTexture{.reflectance = Vector3{1.0,1.0,1.0}}, .ref_index = dielec->inEta, .exponent = dielec->outEta});
        } else if (auto *plastic = std::get_if<ParsedPlastic>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&plastic->reflectance)){
                materials.push_back(Material{.material_type = material_e::PlasticType, .reflectance = SolidTexture{.reflectance = *tex}, .ref_index = plastic->eta});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&plastic->reflectance)){
                materials.push_back(Material{.material_type = material_e::PlasticType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .ref_index = plastic->eta});
            }
        } else if (auto *phong = std::get_if<ParsedPhong>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&phong->reflectance)){
                materials.push_back(Material{.material_type = material_e::PhongType, .reflectance = SolidTexture{.reflectance = *tex}, .exponent=phong->exponent});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&phong->reflectance)){
                materials.push_back(Material{.material_type = material_e::PhongType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .exponent=phong->exponent});
            }
        } else if (auto *blinnphong = std::get_if<ParsedBlinnPhong>(&parsed_mat)){
            if (auto *tex = std::get_if<Vector3>(&blinnphong->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongType, .reflectance = SolidTexture{.reflectance = *tex}, .exponent=blinnphong->exponent});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&blinnphong->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .exponent=blinnphong->exponent});
            }
        } else if (auto *blinnphongmicrofacet = std::get_if<ParsedBlinnPhongMicrofacet>(&parsed_mat)) {
            if (auto *tex = std::get_if<Vector3>(&blinnphongmicrofacet->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongMicrofacetType, .reflectance = SolidTexture{.reflectance = *tex}, .exponent=blinnphongmicrofacet->exponent});
            } else if (auto *tex = std::get_if<ParsedImageTexture>(&blinnphongmicrofacet->reflectance)){
                materials.push_back(Material{.material_type = material_e::BlinnPhongMicrofacetType, .reflectance = ImgTexture(tex->filename, tex->uscale, tex->vscale, tex->uoffset, tex->voffset), .exponent=blinnphongmicrofacet->exponent});
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
                lights.push_back(AreaLight{.radiance = area_light->radiance, .shape_idx = s});
            }
        } else {
            assert(false);
        }
    }
}
