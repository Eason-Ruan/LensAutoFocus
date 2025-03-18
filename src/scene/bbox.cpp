#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  // calculate the intersection with three slabs
  std::vector<double> t_min;
  std::vector<double> t_max;
  t_min.reserve(3);
  t_max.reserve(3);
  if (r.d.x == 0.) {
    t_min[0] = - std::numeric_limits<double>::infinity();
    t_max[0] = std::numeric_limits<double>::infinity();
  }else {
    t_min[0] = (min.x - r.o.x) / r.d.x;
    t_max[0] = (max.x - r.o.x) / r.d.x;
    if (t_min[0] > t_max[0]) {std::swap(t_min[0], t_max[0]);}
  }
  if (r.d.y == 0.) {
    t_min[1] = - std::numeric_limits<double>::infinity();
    t_max[1] = std::numeric_limits<double>::infinity();
  }else {
    t_min[1] = (min.y - r.o.y) / r.d.y;
    t_max[1] = (max.y - r.o.y) / r.d.y;
    if (t_min[1] > t_max[1]) {std::swap(t_min[1], t_max[1]);}
  }
  if (r.d.z == 0.) {
    t_min[2] = - std::numeric_limits<double>::infinity();
    t_max[2] = std::numeric_limits<double>::infinity();
  } else {
    t_min[2] = (min.z - r.o.z) / r.d.z;
    t_max[2] = (max.z - r.o.z) / r.d.z;
    if (t_min[2] > t_max[2]) {std::swap(t_min[2], t_max[2]);}
  }
  std::sort(t_min.begin(), t_min.end());
  std::sort(t_max.begin(), t_max.end());
  const double t_begin = std::max(t_min[2], t0);
  const double t_end = std::min(t_max[0], t1);
  if (t_begin > t_end){return false;}
  t0 = t_begin;
  t1 = t_end;
  return true;
}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
