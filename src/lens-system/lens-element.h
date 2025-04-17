#ifndef CGL_LENS_ELEMENT_H
#define CGL_LENS_ELEMENT_H

#include <cmath>
#include <memory>

#include "CGL/vector3D.h"
#include "CGL/ray.h"
#include "CGL/grid-data.h"
#include "lens-system/ior.h"

namespace CGL {

inline Vector3D align_normal(const Vector3D& v, const Vector3D& n) {
  return dot(v, n) < 0 ? n : -n;
}

class LensElement {
 public:
  unsigned int index;
  double curvature_radius;
  double aperture_radius;
  double thickness;
  double z;

  std::shared_ptr<IOREquation> ior_equation;
  bool is_stop;

  LensElement(unsigned int index, double aperture_radius, double thickness,
              double curvature_radius, const std::shared_ptr<IOREquation>& ior_equation,
              bool is_stop)
      : index(index),
        curvature_radius(curvature_radius),
        aperture_radius(aperture_radius),
        thickness(thickness),
        z(0),
        ior_equation(ior_equation),
        is_stop(is_stop) {}

  bool intersect(const Ray& ray, Hit& res) const {
    // If the element is an aperture or the curvature radius is too large, treat it as a plane
    if (is_stop || curvature_radius > 10000) {
      double t = -(ray.origin.z - z) / ray.direction.z;
      Vector3D hit_pos = ray.at(t);

      double r = hit_pos.x * hit_pos.x + hit_pos.y * hit_pos.y;
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hit_pos = ray.at(t);
      res.hit_normal = align_normal(ray.direction, Vector3D(0, 0, -1));
      return true;
    } else {
      Vector3D center(0, 0, z + curvature_radius);
      double b = dot(ray.origin - center, ray.direction);
      double c = (ray.origin - center).norm2() - curvature_radius * curvature_radius;
      double discriminant = b * b - c;
      if (discriminant < 0) return false;

      double t0 = -b - std::sqrt(discriminant);
      double t1 = -b + std::sqrt(discriminant);
      double t = curvature_radius * ray.direction.z > 0 ? t0 : t1;
      Vector3D hit_pos = ray.at(t);

      double r = hit_pos.x * hit_pos.x + hit_pos.y * hit_pos.y;
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hit_pos = hit_pos;
      res.hit_normal = align_normal(ray.direction, (hit_pos - center).unit());
      return true;
    }
  }

  double ior(double lambda) const { return ior_equation->ior(lambda); }

  GridData<Vector3D> sample_points(unsigned int N) const {
    GridData<Vector3D> ret(N, N);

    // Grid sampling
    for (unsigned int j = 0; j < N; ++j) {
      double v = (2.0 * (j + 0.5) - N) / N;
      double y = v * aperture_radius;
      for (unsigned int i = 0; i < N; ++i) {
        double u = (2.0 * (i + 0.5) - N) / N;
        double x = u * aperture_radius;
        ret.set(j, i, Vector3D(x, y, z));
      }
    }

    return ret;
  }
};

}  // namespace CGL

#endif  // CGL_LENS_ELEMENT_H