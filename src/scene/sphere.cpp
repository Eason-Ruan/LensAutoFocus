#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.

  return true;

}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  const auto d = r.d;
  const auto o = r.o;
  const auto centerSphere  = this->o;
  const auto a = dot(d, d);
  const auto b = dot(2 * (o - centerSphere ), d);
  const auto c = dot(o - centerSphere, o - centerSphere) - r2;
  const auto delta = b * b - 4 * a * c;
  if (delta < 0.) {return false;}
  const auto t_small = (-b - sqrt(delta)) / (2 * a);
  if (t_small < r.min_t || t_small > r.max_t) {return false;}
  // here we confirm intersection
  r.max_t = t_small;
  return true;
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  const auto d = r.d;
  const auto o = r.o;
  const auto centerSphere  = this->o;
  const auto a = dot(d, d);
  const auto b = dot(2 * (o - centerSphere ), d);
  const auto c = dot(o - centerSphere, o - centerSphere) - r2;
  const auto delta = b * b - 4 * a * c;
  if (delta < 0.) {return false;}
  const auto t_small = (-b - sqrt(delta)) / (2 * a);
  if (t_small < r.min_t || t_small > r.max_t) {return false;}
  // here we confirm intersection
  r.max_t = t_small;
  i->t = t_small;
  i->n = normal(o + t_small * d);
  i->primitive = this;
  i->bsdf = get_bsdf();
  return true;
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
