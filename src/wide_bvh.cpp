#include "wide_bvh.h"
#include "3rdparty/pcg.h"
#include <immintrin.h>

pcg32_state wide_bvh_state = init_pcg32();

void set_surrounding_wide_aabb(std::shared_ptr<WideAABB> &parent, const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b){
    parent->a = Vector3{
        min(a->a.x, b->a.x),
        min(a->a.y, b->a.y),
        min(a->a.z, b->a.z)
    };
    parent->b = Vector3{
        max(a->b.x, b->b.x),
        max(a->b.y, b->b.y),
        max(a->b.z, b->b.z)
    };
    parent->shapes = std::vector<Shape*>();
    parent->shapes.insert(parent->shapes.end(), a->shapes.begin(), a->shapes.end());
    parent->shapes.insert(parent->shapes.end(), b->shapes.begin(), b->shapes.end());
    return;
}
void set_surrounding_wide_aabb(std::shared_ptr<WideAABB> &parent, const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b, const std::shared_ptr<WideAABB> &c){
    parent->a = min(min(a->a, b->a), c->a);
    parent->b = max(max(a->b, b->b), c->b);
    parent->shapes = std::vector<Shape*>();
    parent->shapes.insert(parent->shapes.end(), a->shapes.begin(), a->shapes.end());
    parent->shapes.insert(parent->shapes.end(), b->shapes.begin(), b->shapes.end());
    parent->shapes.insert(parent->shapes.end(), c->shapes.begin(), c->shapes.end());
    return;
}
void set_surrounding_wide_aabb(std::shared_ptr<WideAABB> &parent, const std::shared_ptr<WideAABB> &a, const std::shared_ptr<WideAABB> &b, const std::shared_ptr<WideAABB> &c, const std::shared_ptr<WideAABB> &d){
    parent->a = min(min(min(a->a, b->a), c->a), d->a);
    parent->b = max(max(max(a->b, b->b), c->b), d->b);
    parent->shapes = std::vector<Shape*>();
    parent->shapes.insert(parent->shapes.end(), a->shapes.begin(), a->shapes.end());
    parent->shapes.insert(parent->shapes.end(), b->shapes.begin(), b->shapes.end());
    parent->shapes.insert(parent->shapes.end(), c->shapes.begin(), c->shapes.end());
    parent->shapes.insert(parent->shapes.end(), d->shapes.begin(), d->shapes.end());
    return;
}

std::shared_ptr<WideAABB> build_wide_bvh_helper(std::vector<std::shared_ptr<WideAABB>> &boxes, size_t start, size_t end){
    // base case
    if (end-start == 1) {
        return boxes.at(start);
    } 
    std::shared_ptr<WideAABB> node = std::make_shared<WideAABB>();
    uint32_t axis = 0;

    Real min_x = (*std::min_element(boxes.begin()+start, boxes.begin()+end, wide_aabb_x_compare))->a.x;
    Real min_y = (*std::min_element(boxes.begin()+start, boxes.begin()+end, wide_aabb_y_compare))->a.y;
    Real min_z = (*std::min_element(boxes.begin()+start, boxes.begin()+end, wide_aabb_z_compare))->a.z;
    Real max_x = (*std::max_element(boxes.begin()+start, boxes.begin()+end, wide_aabb_x_compare_max))->b.x;
    Real max_y = (*std::max_element(boxes.begin()+start, boxes.begin()+end, wide_aabb_y_compare_max))->b.y;
    Real max_z = (*std::max_element(boxes.begin()+start, boxes.begin()+end, wide_aabb_z_compare_max))->b.z;
    Real z_diff = max_z - min_z;
    Real y_diff = max_y - min_y;
    Real x_diff = max_x - min_x;
    if (z_diff > y_diff && z_diff > x_diff){
        axis = 2;
    } else if (y_diff > z_diff && y_diff > x_diff){
        axis = 1;
    } else {
        axis = 0;
    }

    auto comparator = (axis == 0) ? wide_aabb_x_compare
                    : (axis == 1) ? wide_aabb_y_compare
                                  : wide_aabb_z_compare;
    // other base cases
    if (end - start == 2){
        if (comparator(boxes.at(start), boxes.at(start+1))) {
            node->c1 = boxes.at(start);
            node->c2 = boxes.at(start+1);
        } else {
            node->c1 = boxes.at(start+1);
            node->c2 = boxes.at(start);
        }
        set_surrounding_wide_aabb(node, node->c1, node->c2);
        return node;
    } else if (end - start == 3){
        std::sort(boxes.begin()+start, boxes.begin()+end, comparator);
        node->c1 = boxes.at(start);
        node->c2 = boxes.at(start+1);
        node->c3 = boxes.at(start+2);
        set_surrounding_wide_aabb(node, node->c1, node->c2, node->c3);
        return node;
    } else if (end - start == 4){
        std::sort(boxes.begin()+start, boxes.begin()+end, comparator);
        node->c1 = boxes.at(start);
        node->c2 = boxes.at(start+1);
        node->c3 = boxes.at(start+2);
        node->c4 = boxes.at(start+3);
        set_surrounding_wide_aabb(node, node->c1, node->c2, node->c3, node->c4);
        return node;
    } else {
        //sort subvector on chosen axis
        std::sort(boxes.begin()+start, boxes.begin()+end, comparator);
        // split apart and recurse
        size_t split_size = (end-start)/4;
        node->c1 = std::shared_ptr<WideAABB>(build_wide_bvh_helper(boxes, start, start+split_size));
        node->c2 = std::shared_ptr<WideAABB>(build_wide_bvh_helper(boxes, start+split_size, start+split_size*2));
        node->c3 = std::shared_ptr<WideAABB>(build_wide_bvh_helper(boxes, start+split_size*2, start+split_size*3));
        node->c4 = std::shared_ptr<WideAABB>(build_wide_bvh_helper(boxes, start+split_size*3, end));
        set_surrounding_wide_aabb(node, node->c1, node->c2, node->c3, node->c4);
        return node;
    }
}

std::shared_ptr<WideAABB> build_wide_bvh(std::vector<Shape> &shapes){
    std::vector<std::shared_ptr<WideAABB>> boxes;
    Vector3 p0, p1, p2;

    for (int shape = 0; shape < shapes.size(); ++shape){
        std::shared_ptr<WideAABB> box = std::make_shared<WideAABB>();
        if (auto *sph = std::get_if<Sphere>(&shapes.at(shape))){
            box->a = sph->center - sph->radius;
            box->b = sph->center + sph->radius;
        } else if (auto *tri = std::get_if<Triangle>(&shapes.at(shape))) {
            p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
            p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
            p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
            box->a = Vector3{min(p2.x, min(p0.x, p1.x)), min(p2.y, min(p0.y, p1.y)), min(p2.z, min(p0.z, p1.z))};
            box->b = Vector3{max(p2.x, max(p0.x, p1.x)), max(p2.y, max(p0.y, p1.y)), max(p2.z, max(p0.z, p1.z))};
        } else {
            assert(false);
        }
        box->shapes = std::vector<Shape*>{&shapes.at(shape)};
        boxes.push_back(box);
    }
    return std::shared_ptr<WideAABB>(build_wide_bvh_helper(boxes, 0, boxes.size()));
}

bool hit_box(const WideAABB &box, const Vector3 &ray, const Vector3 &lookfrom, Real t_min, Real t_max){
    for (int d = 0; d < 3; d++) {
        Real invD = 1.0 / ray[d];
        auto t0 = (box.a[d] - lookfrom[d]) * invD;
        auto t1 = (box.b[d] - lookfrom[d]) * invD;
        if (invD < 0.0)
            std::swap(t0, t1);
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if (t_max < t_min){
            return false;
        }
    }
    return true;
}


int hit_wide_box_simd(const WideAABB &box, const Vector3 &ray, const Vector3 &lookfrom, Real t_min, Real t_max){
    __m256d t0, t1, m256_lookfrom, invD, vect_min, vect_max, hitacc, temp;
    vect_min = _mm256_set1_pd(t_min);
    vect_max = _mm256_set1_pd(t_max);
    hitacc = _mm256_setzero_pd();
    int res;
    for (int d = 0; d < 3; d++) {
        Real real_invD = 1.0 / ray[d];
        m256_lookfrom = _mm256_set1_pd(lookfrom[d]);
        invD = _mm256_set1_pd(real_invD);
        t0 = _mm256_set_pd(box.c1->a[d], box.c2->a[d], box.c3->a[d], box.c4->a[d]);
        t1 = _mm256_set_pd(box.c1->b[d], box.c2->b[d], box.c3->b[d], box.c4->b[d]);
        t0 = _mm256_sub_pd(t0, m256_lookfrom);
        t1 = _mm256_sub_pd(t1, m256_lookfrom);
        t0 = _mm256_mul_pd(t0, invD);
        t1 = _mm256_mul_pd(t1, invD);
        if (real_invD < 0.0){
            std::swap(t0, t1);
        }
        vect_min = _mm256_max_pd(t0, vect_min);
        vect_max = _mm256_min_pd(t1, vect_max);
        temp = _mm256_cmp_pd(vect_max, vect_min, 0x1);
        hitacc = _mm256_or_pd(temp, hitacc);
        res = _mm256_movemask_pd(hitacc);
        if (res == 15){
            return res;
        }
    }
    return res;
}

bool hit_bvh(const WideAABB &bvh, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, Shape **hit_shape){
    // base case
    if (bvh.shapes.size()==1){
       // std::cout << "hit base" << std::endl;
        // Real temp_t_val;
        // Shape* temp_shape = bvh.shapes.at(0);
        // if (auto *sph = std::get_if<Sphere>(temp_shape)){
        //     temp_t_val = hit_sphere(*sph, ray, ray_origin);
        //     if (temp_t_val > 0.0){
        //         if ((t_val < 0.0 || temp_t_val < t_val)&& distance(ray_origin, ray_origin + ray * temp_t_val)> eps){
        //             *hit_shape = temp_shape;
        //             t_val = temp_t_val;
        //             return true;
        //         } else {
        //             return false;
        //         }
        //     } else {
        //         return false;
        //     }
        // } else if (auto *tri = std::get_if<Triangle>(temp_shape)){
        //     Vector3 p0 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).x);
        //     Vector3 p1 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).y);
        //     Vector3 p2 = tri->mesh->positions.at(tri->mesh->indices.at(tri->face_index).z);
        //     if ( get_tri_intersect(p0, p1, p2, ray, ray_origin, eps, b_coord, temp_t_val)){
        //         if ((t_val < 0.0 || temp_t_val < t_val)&& distance(ray_origin, ray_origin + ray * temp_t_val)> eps){
        //             *hit_shape = temp_shape;
        //             t_val = temp_t_val;
        //             return true;
        //         } else {
        //             return false;
        //         }
        //     } else {
        //         return false;
        //     }
        // } else {
        //     assert(false);
        // }
        if (hit_box(bvh, ray, ray_origin, eps, infinity<Real>() )){
            *hit_shape = bvh.shapes.at(0);
            return true;
        } else {
            return false;
        }
    } else {
        bool c1_hit = false, c2_hit = false, c3_hit = false, c4_hit = false;
        //if (with_simd){
        int hits = hit_wide_box_simd(bvh, ray, ray_origin, eps, infinity<Real>());
        if (hits == 15){
            return false;
        }
        if ((hits & 0b1000) == 0b1000){// c1 hit
            c1_hit = hit_bvh(*bvh.c1, ray, ray_origin, eps, hit_shape);
        }
        if ((hits & 0b0100) == 0b0100){
            c2_hit = hit_bvh(*bvh.c2, ray, ray_origin, eps, hit_shape);
        }
        if ((hits & 0b0010) == 0b0010){
            c3_hit = hit_bvh(*bvh.c3, ray, ray_origin, eps, hit_shape);
        }
        if ((hits & 0b0001) == 0b0001){
            c4_hit = hit_bvh(*bvh.c4, ray, ray_origin, eps, hit_shape);
        }
        return c1_hit || c2_hit || c3_hit || c4_hit;
        // if (hits == 0){
        //     l_hit = hit_bvh(*bvh.left, ray, ray_origin, eps, hit_shape, t_val, b_coord, with_simd);
        //     r_hit = hit_bvh(*bvh.right, ray, ray_origin, eps, hit_shape, t_val, b_coord, with_simd);
        //     return l_hit || r_hit;
        // } else if (hits==2){
        //     return hit_bvh(*bvh.right, ray, ray_origin, eps, hit_shape, t_val, b_coord, with_simd);
        // }
        // else if (hits==1){
        //     return hit_bvh(*bvh.left, ray, ray_origin, eps, hit_shape, t_val, b_coord, with_simd);
        // } else {
        //     return false;
        // }
        // } else {
        //     if (!hit_box(bvh, ray, ray_origin, eps, infinity<Real>())){
        //         return false;
        //     } else {
        //         l_hit = hit_bvh(*bvh.left, ray, ray_origin, eps, hit_shape, t_val, b_coord, with_simd);
        //         r_hit = hit_bvh(*bvh.right, ray, ray_origin, eps, hit_shape, t_val, b_coord, with_simd);
        //         return l_hit || r_hit;
        //     }
        // }
    }
}