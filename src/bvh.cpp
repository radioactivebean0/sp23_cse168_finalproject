#include "bvh.h"
#include "shapes.h"
#include "3rdparty/pcg.h"
#include <immintrin.h>
#include <emmintrin.h>

#define NUM_BINS 4

pcg32_state bvh_state = init_pcg32();

void set_surrounding_aabb(std::shared_ptr<AABB> &parent, const std::shared_ptr<AABB> &a, const std::shared_ptr<AABB> &b){
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

Real get_union_aabb_surface(std::vector<std::shared_ptr<AABB>> &boxes, size_t start, size_t end){
    __m256d min_vec = _mm256_set_pd(boxes.at(start)->a.x, boxes.at(start)->a.y, boxes.at(start)->a.z, 0.0);
    __m256d max_vec = _mm256_set_pd(boxes.at(start)->b.x,boxes.at(start)->b.y,boxes.at(start)->b.z,0.0);
    __m256d tempa, tempb;
    for (int i = start+1; i < end; i++){
        tempa = _mm256_set_pd(boxes[i]->a.x, boxes[i]->a.y, boxes[i]->a.z, 0.0);
        tempb = _mm256_set_pd(boxes[i]->b.x, boxes[i]->b.y, boxes[i]->b.z, 0.0);
        min_vec = _mm256_min_pd(tempa, min_vec);
        max_vec = _mm256_max_pd(tempb, max_vec);
    }
    Real diag[4];
    min_vec = _mm256_sub_pd(max_vec, min_vec);
    _mm256_store_pd(diag, min_vec);
    return 2 * (diag[3] * diag[2] + diag[3]*diag[1] + diag[2]*diag[1]);
}

void get_bounds_aabb_list(std::vector<std::shared_ptr<AABB>> &boxes, size_t start, size_t end, Vector3 &a, Vector3 &b){
    __m256d tempa, tempb;
    __m256d min_vec = _mm256_set_pd(boxes.at(start)->a.x, boxes.at(start)->a.y, boxes.at(start)->a.z, 0.0);
    __m256d max_vec = _mm256_set_pd(boxes.at(start)->b.x,boxes.at(start)->b.y,boxes.at(start)->b.z,0.0);
    for (int i = start+1; i < end; i++){
        tempa = _mm256_set_pd(boxes[i]->a.x, boxes[i]->a.y, boxes[i]->a.z, 0.0);
        tempb = _mm256_set_pd(boxes[i]->b.x, boxes[i]->b.y, boxes[i]->b.z, 0.0);
        min_vec = _mm256_min_pd(tempa, min_vec);
        max_vec = _mm256_max_pd(tempb, max_vec);
    }
    Real extract_a[4];
    Real extract_b[4];
    _mm256_store_pd(extract_a, min_vec);
    _mm256_store_pd(extract_b, max_vec);
    a = Vector3{extract_a[3],extract_a[2],extract_a[1]};
    b = Vector3{extract_b[3],extract_b[2],extract_b[1]};
}

Real surface_area_one(Vector3 a, Vector3 b){
    Vector3 diag = b - a;
    return 2 * (diag.x * diag.y + diag.x * diag.z + diag.y * diag.z);
}

Real surface_area_two(Vector3 a1, Vector3 b1, Vector3 a2, Vector3 b2){
    Vector3 diag = max(b2, b1) - min(a1,a2);
    return 2 * (diag.x * diag.y + diag.x * diag.z + diag.y * diag.z);
}

Real surface_area_three(Vector3 a1, Vector3 b1, Vector3 a2, Vector3 b2, Vector3 a3, Vector3 b3){
    Vector3 diag = max(max(b2, b1),b3) - min(min(a1,a2),a3);
    return 2 * (diag.x * diag.y + diag.x * diag.z + diag.y * diag.z);
}

uint32_t get_best_four_split(Vector3 a[], Vector3 b[], Real &sah, size_t chunksize){
    Real temp_sah;
    uint32_t splitind = 1;
    sah = (surface_area_one(a[0],b[0])*chunksize)+ (surface_area_three(a[1],b[1],a[2],b[2],a[3],b[3])*3*chunksize);
    temp_sah = (surface_area_two(a[0],b[0],a[1],b[1])+surface_area_two(a[2],b[2],a[3],b[3]))*(2*chunksize);
    if (temp_sah < sah){
        splitind = 2;
        sah = temp_sah;
    }
    temp_sah = (surface_area_one(a[3],b[3])*chunksize)+ (surface_area_three(a[0],b[0],a[1],b[1],a[2],b[2])*(3*chunksize));
    if (temp_sah < sah){
        splitind = 3;
        sah = temp_sah;
    }
    return splitind;
}

Real get_cost_ranges(std::vector<std::shared_ptr<AABB>> &boxes, size_t start, size_t mid, size_t end){
    return (get_union_aabb_surface(boxes, start, mid)*(mid-start))+(get_union_aabb_surface(boxes, mid, end)*(end-mid))/get_union_aabb_surface(boxes, start, end);
}

std::shared_ptr<AABB> build_sah_bvh_helper(std::vector<std::shared_ptr<AABB>> &boxes, size_t start, size_t end, int *bvhsize){
    // base cases
    if (end-start == 1) {
        return boxes.at(start);
    } 
    std::shared_ptr<AABB> node = std::make_shared<AABB>();
    uint32_t axis = 0;
    if (end - start ==2){
        Real min_x = (*std::min_element(boxes.begin()+start, boxes.begin()+end, aabb_x_compare))->a.x;
        Real min_y = (*std::min_element(boxes.begin()+start, boxes.begin()+end, aabb_y_compare))->a.y;
        Real min_z = (*std::min_element(boxes.begin()+start, boxes.begin()+end, aabb_z_compare))->a.z;
        Real max_x = (*std::max_element(boxes.begin()+start, boxes.begin()+end, aabb_x_compare_max))->b.x;
        Real max_y = (*std::max_element(boxes.begin()+start, boxes.begin()+end, aabb_y_compare_max))->b.y;
        Real max_z = (*std::max_element(boxes.begin()+start, boxes.begin()+end, aabb_z_compare_max))->b.z;
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
            auto comparator = (axis == 0) ? aabb_x_compare
                    : (axis == 1) ? aabb_y_compare
                                    : aabb_z_compare;
        if (comparator(boxes.at(start), boxes.at(start+1))) {
            node->left = boxes.at(start);
            node->right = boxes.at(start+1);
        } else {
            node->left = boxes.at(start+1);
            node->right = boxes.at(start);
        }
        set_surrounding_aabb(node, node->left, node->right);
        *bvhsize += 1;
        return node;
    }

    uint32_t split_ind = -1;
    Real best_sah = infinity<Real>();
    Real temp_sah;

    if (end-start == 3){
        // sort by x
        std::sort(boxes.begin()+start, boxes.begin()+end, aabb_x_compare);
        best_sah = get_cost_ranges(boxes, start, start+1, start+3);
        split_ind = 1;
        temp_sah = get_cost_ranges(boxes, start, start+2, start+3);
        if (temp_sah < best_sah){
            best_sah = temp_sah;
            split_ind = 2;
        }
        // sort by y
        std::sort(boxes.begin()+start, boxes.begin()+end, aabb_y_compare);
        temp_sah = get_cost_ranges(boxes, start, start+1, start+3);
        if (temp_sah < best_sah){
            split_ind = 1;
            axis = 1;
            best_sah = temp_sah;
        }
        temp_sah = get_cost_ranges(boxes, start, start+2, start+3);
        if (temp_sah < best_sah){
            split_ind = 2;
            axis = 1;
            best_sah = temp_sah;
        }
        // sort by z
        std::sort(boxes.begin()+start, boxes.begin()+end, aabb_z_compare);
        temp_sah = get_cost_ranges(boxes, start, start+1, start+3);
        if (temp_sah < best_sah){
            split_ind = 1;
            axis = 2;
            best_sah = temp_sah;
        }
        temp_sah = get_cost_ranges(boxes, start, start+2, start+3);;
        if (temp_sah < best_sah){
            split_ind = 2;
            axis = 2;
            best_sah = temp_sah;
        }
        split_ind = start + split_ind;
    } else { // split into even 4 on each axis and try
        int chunksize = (end - start)/NUM_BINS;
        Vector3 chunks_a[4];
        Vector3 chunks_b[4];
        int i;
        uint32_t temp_splitind;
        Real outerbox_surface = get_union_aabb_surface(boxes, start, end);
        for (int d = 0; d < 3; d++){
            auto comparator = (d == 0) ? aabb_x_compare
                            : (d == 1) ? aabb_y_compare
                                       : aabb_z_compare;
            for (i = 1; i < NUM_BINS; i++){
                std::nth_element(boxes.begin()+start, boxes.begin()+(chunksize*i), boxes.begin()+end, comparator);
            }
            i = 0;
            while (i < NUM_BINS-1){
                get_bounds_aabb_list(boxes, start+(chunksize*i), start+(i+1)*chunksize, chunks_a[i], chunks_b[i]);
                i++;
            }
            get_bounds_aabb_list(boxes, start+(chunksize)*i, end, chunks_a[i], chunks_b[i]);

            temp_splitind = get_best_four_split(chunks_a, chunks_b, temp_sah, chunksize);
            temp_sah = temp_sah /outerbox_surface;
            if (temp_sah < best_sah){
                axis = d;
                best_sah = temp_sah;
                split_ind = temp_splitind;
            }
        }
        split_ind = start + (chunksize) * split_ind;
    }
    auto comparator = (axis == 0) ? aabb_x_compare
                    : (axis == 1) ? aabb_y_compare
                                    : aabb_z_compare;
    // sort subvector on chosen axis
    std::nth_element(boxes.begin()+start, boxes.begin()+split_ind, boxes.begin()+end, comparator);
    // split apart and recurse
    node->left = build_sah_bvh_helper(boxes, start, split_ind, bvhsize);
    node->right = build_sah_bvh_helper(boxes, split_ind, end, bvhsize);
    set_surrounding_aabb(node, node->left, node->right);
    *bvhsize += 1;
    return node;
}

std::shared_ptr<AABB> build_bvh_helper(std::vector<std::shared_ptr<AABB>> &boxes, size_t start, size_t end, int split_method){
    // base cases
    if (end-start == 1) {
        return boxes.at(start);
    } 
    std::shared_ptr<AABB> node = std::make_shared<AABB>();
    uint32_t axis = 0;
    if (split_method == 0) { // random split
        axis = nextUInt(bvh_state, 3);
    } else if (split_method == 1){ // split axis with largest extent
        Real min_x = (*std::min_element(boxes.begin()+start, boxes.begin()+end, aabb_x_compare))->a.x;
        Real min_y = (*std::min_element(boxes.begin()+start, boxes.begin()+end, aabb_y_compare))->a.y;
        Real min_z = (*std::min_element(boxes.begin()+start, boxes.begin()+end, aabb_z_compare))->a.z;
        Real max_x = (*std::max_element(boxes.begin()+start, boxes.begin()+end, aabb_x_compare_max))->b.x;
        Real max_y = (*std::max_element(boxes.begin()+start, boxes.begin()+end, aabb_y_compare_max))->b.y;
        Real max_z = (*std::max_element(boxes.begin()+start, boxes.begin()+end, aabb_z_compare_max))->b.z;
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
    } else {
        assert(false);
    }
    // select random axis
    auto comparator = (axis == 0) ? aabb_x_compare
                    : (axis == 1) ? aabb_y_compare
                                  : aabb_z_compare;
    node->axis = axis;
    if (end - start == 2){
        if (comparator(boxes.at(start), boxes.at(start+1))) {
            node->left = boxes.at(start);
            node->right = boxes.at(start+1);
        } else {
            node->left = boxes.at(start+1);
            node->right = boxes.at(start);
        }
        set_surrounding_aabb(node, node->left, node->right);
        return node;
    } else {
        // sort subvector on chosen axis
        std::sort(boxes.begin()+start, boxes.begin()+end, comparator);
        // split apart and recurse
        size_t middle = start + (end-start)/2;
        node->left = build_bvh_helper(boxes, start, middle, split_method);
        node->right = build_bvh_helper(boxes, middle, end, split_method);
        set_surrounding_aabb(node, node->left, node->right);
        return node;
    }
}


AABB build_bvh(std::vector<Shape> &shapes, int split_method = 0, int *bvhsize = nullptr){
    std::vector<std::shared_ptr<AABB>> boxes;
    Vector3 p0, p1, p2;

    for (int shape = 0; shape < shapes.size(); ++shape){
        std::shared_ptr<AABB> box = std::make_shared<AABB>();
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
    if (split_method == 2){
        *bvhsize = boxes.size();
        return *build_sah_bvh_helper(boxes, 0, boxes.size(), bvhsize);
    }
    return *build_bvh_helper(boxes, 0, boxes.size(), split_method);
}
AABB build_bvh(std::vector<Shape> &shapes, int split_method = 0){
    int dummy = 0;
    return build_bvh(shapes, split_method, &dummy);
}
int flatten_bvh(AABB &box, compact_AABB *root, int *offset){
    compact_AABB *cnode = &root[*offset];
    cnode->a = box.a;
    cnode->b = box.b;
    int myoffs = (*offset)++;
    //std::cout << box.shapes.size() << std::endl;
    if (box.shapes.size()==1){
        cnode->shape = box.shapes.at(0);
        cnode->next = -1;
    } else {
        cnode->shape = nullptr;
        flatten_bvh(*box.left, root, offset);
        cnode->next = flatten_bvh(*box.right, root, offset);
    }
    cnode->axis = box.axis;
    return myoffs;
}

int hit_boxes_simd(const AABB &box, const Vector3 &ray, const Vector3 &lookfrom, Real t_min, Real t_max){
    __m128d t0, t1, m128_lookfrom, invD, vect_min, vect_max, hitacc, temp;
    vect_min = _mm_set_pd1(t_min);
    vect_max = _mm_set_pd1(t_max);
    hitacc = _mm_setzero_pd();
    int res;
    for (int d = 0; d < 3; d++) {
        Real real_invD = 1.0 / ray[d];
        m128_lookfrom = _mm_set_pd1(lookfrom[d]);
        invD = _mm_set_pd1(real_invD);
        t0 = _mm_set_pd(box.left->a[d], box.right->a[d]);
        t1 = _mm_set_pd(box.left->b[d], box.right->b[d]);
        t0 = _mm_sub_pd(t0, m128_lookfrom);
        t1 = _mm_sub_pd(t1, m128_lookfrom);
        t0 = _mm_mul_pd(t0, invD);
        t1 = _mm_mul_pd(t1, invD);
        if (real_invD < 0.0){
            std::swap(t0, t1);
        }
        vect_min = _mm_max_pd(t0, vect_min);
        vect_max = _mm_min_pd(t1, vect_max);
        temp = _mm_cmplt_pd(vect_max, vect_min);
        hitacc = _mm_or_pd(temp, hitacc);
        res = _mm_movemask_pd(hitacc);
        if (res == 3){
            return res;
        }
    }
    return res;
}

bool hit_box(const AABB &box, const Vector3 &ray, const Vector3 &lookfrom, Real t_min, Real t_max){
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

bool hit_bvh_old(const AABB &bvh, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, Real t_min, Real t_max, Shape **hs, Real &t, Vector2 &uv){
    if (bvh.shapes.size()==1){
        Real temp_t_val;
        temp_t_val = hit_shape(bvh.shapes.at(0), ray, ray_origin, eps, uv);
        if (temp_t_val > 0.0){
            if ((t < 0.0 || temp_t_val < t)&& distance(ray_origin, ray_origin + ray * temp_t_val)> eps){
                *hs = bvh.shapes.at(0);
                t = temp_t_val;
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        bool l_hit = false, r_hit = false;
        if (!hit_box(bvh, ray, ray_origin, t_min, t_max)){
            return false;
        } else {
            l_hit = hit_bvh_old(*bvh.left, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
            r_hit = hit_bvh_old(*bvh.right, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
            return l_hit || r_hit;
        }
    }
}

bool hit_bvh(const AABB &bvh, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, Real t_min, Real t_max,  Shape **hs, Real &t, Vector2 &uv, const bool with_simd){
    //std::cout << "hitbvh" << std::endl;
    if (!with_simd){
        return hit_bvh_old(bvh, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
    }
    // base case
    if (bvh.shapes.size()==1){
        Real temp_t_val;
        Vector2 temp_uv;
        temp_t_val = hit_shape(bvh.shapes.at(0), ray, ray_origin, eps, temp_uv);
        if (temp_t_val > 0.0){
            if ((t < 0.0 || temp_t_val < t)&& distance(ray_origin, ray_origin + ray * temp_t_val)> eps){
                *hs = bvh.shapes.at(0);
                t = temp_t_val;
                uv = temp_uv;
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        bool l_hit = false, r_hit = false;
        int hits = hit_boxes_simd(bvh, ray, ray_origin, t_min, t_max);
        if (hits == 3){
            return false;
        }
        if (hits == 0){
            l_hit = hit_bvh(*bvh.left, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
            r_hit = hit_bvh(*bvh.right, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
            return l_hit || r_hit;
        } else if (hits==2){
            return hit_bvh(*bvh.right, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
        }
        else if (hits==1){
            return hit_bvh(*bvh.left, ray, ray_origin, eps, t_min, t_max, hs, t, uv);
        } else {
            return false;
        }
    }
}

bool hit_bvh(const AABB &bvh, const Vector3 &ray, const Vector3 &ray_origin, const Real eps,  Real t_min, Real t_max,  Shape **hs,Real &t, Vector2 &uv){
    return hit_bvh(bvh, ray, ray_origin, eps, t_min, t_max, hs, t, uv, true);
}

inline int hit_cboxes_simd(const Vector3 &la, const Vector3 &lb, const Vector3 &ra, const Vector3 &rb,
        const __m128d lookfromx, const __m128d lookfromy, const __m128d lookfromz,
        const __m128d invdx, const __m128d invdy, const __m128d invdz, const Vector3 invD,
        __m128d vect_min, __m128d vect_max){
    __m128d hitacc = _mm_setzero_pd();
    __m128d t0, t1, temp;
    int res;
    t0 = _mm_set_pd(la[0], ra[0]);
    t1 = _mm_set_pd(lb[0], rb[0]);
    t0 = _mm_sub_pd(t0, lookfromx);
    t1 = _mm_sub_pd(t1, lookfromx);
    t0 = _mm_mul_pd(t0, invdx);
    t1 = _mm_mul_pd(t1, invdx);
    if (invD[0] < 0.0){
        std::swap(t0, t1);
    }
    vect_min = _mm_max_pd(t0, vect_min);
    vect_max = _mm_min_pd(t1, vect_max);
    temp = _mm_cmplt_pd(vect_max, vect_min);
    hitacc = _mm_or_pd(temp, hitacc);
    res = _mm_movemask_pd(hitacc);
    if (res == 3){
        return res;
    }
    t0 = _mm_set_pd(la[1], ra[1]);
    t1 = _mm_set_pd(lb[1], rb[1]);
    t0 = _mm_sub_pd(t0, lookfromy);
    t1 = _mm_sub_pd(t1, lookfromy);
    t0 = _mm_mul_pd(t0, invdy);
    t1 = _mm_mul_pd(t1, invdy);
    if (invD[1] < 0.0){
        std::swap(t0, t1);
    }
    vect_min = _mm_max_pd(t0, vect_min);
    vect_max = _mm_min_pd(t1, vect_max);
    temp = _mm_cmplt_pd(vect_max, vect_min);
    hitacc = _mm_or_pd(temp, hitacc);
    res = _mm_movemask_pd(hitacc);
    if (res == 3){
        return res;
    }
    t0 = _mm_set_pd(la[2], ra[2]);
    t1 = _mm_set_pd(lb[2], rb[2]);
    t0 = _mm_sub_pd(t0, lookfromz);
    t1 = _mm_sub_pd(t1, lookfromz);
    t0 = _mm_mul_pd(t0, invdz);
    t1 = _mm_mul_pd(t1, invdz);
    if (invD[2] < 0.0){
        std::swap(t0, t1);
    }
    vect_min = _mm_max_pd(t0, vect_min);
    vect_max = _mm_min_pd(t1, vect_max);
    temp = _mm_cmplt_pd(vect_max, vect_min);
    hitacc = _mm_or_pd(temp, hitacc);
    res = _mm_movemask_pd(hitacc);
    return res;
}

bool hit_cbvh(const compact_AABB *cbvh, const Vector3 &ray, const Vector3 &ray_origin, const Real eps, Real t_min, Real t_max, Shape **hs,  Real &t, Vector2 &uv){
    bool hit = false;
    __m128d vect_min = _mm_set_pd1(t_min);
    __m128d vect_max = _mm_set_pd1(t_max);
    const __m128d lookfromx = _mm_set_pd1(ray_origin.x);
    const __m128d lookfromy = _mm_set_pd1(ray_origin.y);
    const __m128d lookfromz = _mm_set_pd1(ray_origin.z);
    const Vector3 invD = 1.0/ray;
    const __m128d invdx = _mm_set_pd1(invD.x);
    const __m128d invdy = _mm_set_pd1(invD.y);
    const __m128d invdz = _mm_set_pd1(invD.z);
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[128];
    while (true){
        const compact_AABB *node = &cbvh[currentNodeIndex];
        if (node->next ==-1){
            Real temp_t_val;
            Vector2 temp_uv;
            temp_t_val = hit_shape(node->shape, ray, ray_origin, eps, temp_uv);
            if (temp_t_val > 0.0){
                if ((t < 0.0 || temp_t_val < t)&& distance(ray_origin, ray_origin + ray * temp_t_val)> eps){
                    *hs = node->shape;
                    t = temp_t_val;
                    uv = temp_uv;
                    hit = true;
                }
            }
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        } else {
            int hits = hit_cboxes_simd(cbvh[currentNodeIndex+1].a, cbvh[currentNodeIndex+1].b, cbvh[node->next].a,cbvh[node->next].b,
                                        lookfromx, lookfromy, lookfromz, invdx, invdy, invdz, invD, vect_min, vect_max);
            if (hits == 3){
                //std::cout << "hit none" << std::endl;
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            }
            else if (hits == 0){
                //std::cout << "hit both" << std::endl;
                nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                currentNodeIndex = node->next;
            } else if (hits==2){
                //std::cout << "hit left" << std::endl;
                currentNodeIndex = node->next;
            } else if (hits==1){
                //std::cout << "hit right" << std::endl;
                currentNodeIndex +=1;
            } else {
                assert(false);
            }
        }
    }
    return hit;
}
