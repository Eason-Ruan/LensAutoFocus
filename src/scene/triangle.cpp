#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.
  // use Moller Trumbore Algorithm
  const auto e1 = p2 - p1;
  const auto e2 = p3 - p1;
  const auto s = r.o - p1;
  const auto s1 = cross(r.d, e2);
  const auto s2 = cross(s, e1);
  const auto s1e1 = dot(s1, e1);
  if (s1e1 == 0.) {return false;}
  const auto t =  dot(s2, e2) / s1e1;
  if (t < r.min_t || t > r.max_t) {return false;}
  const auto b1 = dot(s1, s) / s1e1;
  if (b1 < 0.0 || b1 > 1.0) {return false;}
  const auto b2 =dot(s2, r.d) / s1e1;
  if (b2 < 0.0 || b2 > 1.0) {return false;}
  const auto b0 = 1. - b1 - b2;
  if (b0 < 0.0 || b0 > 1.0) {return false;}
  // proved to be intersected here
  r.max_t = t;
  return true;
}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  // use Moller Trumbore Algorithm
  const auto e1 = p2 - p1;
  const auto e2 = p3 - p1;
  const auto s = r.o - p1;
  const auto s1 = cross(r.d, e2);
  const auto s2 = cross(s, e1);
  const auto s1e1 = dot(s1, e1);
  if (s1e1 == 0.) {return false;}
  const auto t =  dot(s2, e2) / s1e1;
  if (t < r.min_t || t > r.max_t) {return false;}
  const auto b1 = dot(s1, s) / s1e1;
  if (b1 < 0.0 || b1 > 1.0) {return false;}
  const auto b2 =dot(s2, r.d) / s1e1;
  if (b2 < 0.0 || b2 > 1.0) {return false;}
  const auto b0 = 1. - b1 - b2;
  if (b0 < 0.0 || b0 > 1.0) {return false;}
  // proved to be intersected here
  r.max_t = t;
  isect->t = t;
  isect->n =  b0 * n1 + b1 * n2 + b2 * n3;
  isect->primitive = this;
  isect->bsdf = get_bsdf();
  return true;
}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
