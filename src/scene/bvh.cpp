#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.
  BBox bbox;
  auto tmp_l = start;
  auto tmp_r = end;
  for (auto p = start; p != end; ++p) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
  }
  auto *node = new BVHNode(bbox);
  const int length = static_cast<int>(end - start);
  if (length > max_leaf_size) {
    // determine axes
    const double x_center = bbox.centroid().x;
    const double y_center = bbox.centroid().y;
    const double z_center = bbox.centroid().z;
    BBox x_left, x_right, y_left, y_right, z_left, z_right;
    int x_left_count = 0;
    int x_right_count = 0;
    int y_left_count = 0;
    int y_right_count = 0;
    int z_left_count = 0;
    int z_right_count = 0;
    for (auto p = start; p != end; ++p) {
      const BBox bb = (*p)->get_bbox();
      const Vector3D center = bb.centroid();
      if (center.x <= x_center) {
        x_left.expand(bb);
        x_left_count++;
      } else {
        x_right.expand(bb);
        x_right_count++;
      }
      if (center.y <= y_center) {
        y_left.expand(bb);
        y_left_count++;
      } else {
        y_right.expand(bb);
        y_right_count++;
      }
      if (center.z <= z_center) {
        z_left.expand(bb);
        z_left_count++;
      } else {
        z_right.expand(bb);
        z_right_count++;
      }
    }
    double cost_x = (x_left.surface_area()* x_left_count
    + x_right.surface_area() * x_right_count) / bbox.surface_area();
    double cost_y = (y_left.surface_area() * y_left_count + y_right.surface_area() * y_right_count) / bbox.surface_area();
    double cost_z = (z_left.surface_area() * z_left_count + z_right.surface_area() * z_right_count) / bbox.surface_area();
    if (cost_x < cost_y && cost_x < cost_z) {
      while (tmp_l < tmp_r) {
        if ((*tmp_l)->get_bbox().centroid().x > x_center) {
          do {
            --tmp_r;
          } while (!((*tmp_r)->get_bbox().centroid().x < x_center) && tmp_l < tmp_r);
          if (tmp_l < tmp_r) {
            iter_swap(tmp_l, tmp_r);
          }else {break;}
        }else {
          ++ tmp_l;
        }
      }
      if (tmp_r == end || tmp_l == start) {
        tmp_l = start + (end - start) / 2;
      }
    }else if  (cost_y < cost_x && cost_y < cost_z) {
      while (tmp_l < tmp_r) {
        if ((*tmp_l)->get_bbox().centroid().y > y_center) {
          do {
            --tmp_r;
          } while (!((*tmp_r)->get_bbox().centroid().y < y_center) && tmp_l < tmp_r);
          if (tmp_l < tmp_r) {
            iter_swap(tmp_l, tmp_r);
          }else {break;}
        }else {
          ++ tmp_l;
        }
      }
      if (tmp_r == end || tmp_l == start) {
        tmp_l = start + (end - start) / 2;
      }
    }else {
      while (tmp_l < tmp_r) {
        if ((*tmp_l)->get_bbox().centroid().z > z_center) {
          do {
            --tmp_r;
          } while (!((*tmp_r)->get_bbox().centroid().z < z_center) && tmp_l < tmp_r);
          if (tmp_l < tmp_r) {
            iter_swap(tmp_l, tmp_r);
          }else {break;}
        }else {
          ++ tmp_l;
        }
      }
      if (tmp_r == end || tmp_l == start) {
        tmp_l = start + (end - start) / 2;
      }
    }
    node->l = construct_bvh(start, tmp_l, max_leaf_size);
    node->r = construct_bvh(tmp_l, end, max_leaf_size);
  } else {
    node->start = start;
    node->end = end;
  }
  return node;
}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  double t_begin = ray.min_t;
  double t_end = ray.max_t;
  if (node->bb.intersect(ray, t_begin, t_end)) {
    if (! node->isLeaf()) {return has_intersection(ray, node->l) || has_intersection(ray, node->r);}
    // leaf case
    for (auto p = node->start; p != node->end; ++p) {
      total_isects++;
      if ((*p)->has_intersection(ray))
        return true;
    }
  }
  return false;
}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  bool hit = false;
  double t_begin = ray.min_t;
  double t_end = ray.max_t;
  if (node->bb.intersect(ray, t_begin, t_end)) {
    if (! node->isLeaf()) {
      bool l_hit = intersect(ray, i, node->l);
      bool r_hit = intersect(ray, i, node->r);
      return l_hit | r_hit;
    }
    // leaf case
    for (auto p = node->start; p != node->end; ++p) {
      total_isects++;
      hit  = (*p)->intersect(ray, i) || hit;
    }
  }
  return hit;
}

} // namespace SceneObjects
} // namespace CGL
