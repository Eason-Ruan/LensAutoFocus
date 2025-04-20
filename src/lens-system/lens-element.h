#ifndef CGL_LENS_ELEMENT_H
#define CGL_LENS_ELEMENT_H

#include <cmath>
#include <memory>

#include "grid_data.h"
#include "CGL/vector3D.h"
#include "ray.h"
#include "lens-system/ior.h"

namespace CGL {

inline Vector3D align_normal(const Vector3D& v, const Vector3D& n) {
  return dot(v, n) < 0 ? n : -n;
}
  struct Hit {
    Real t;
    Vector3D hitPos;
    Vector3D hitNormal;

    Hit(): t(0) {} ;
  };
class LensElement {
 public:
  unsigned int index;
  double curvature_radius;
  double aperture_radius;
  double thickness;
  double z;

  std::shared_ptr<IOREquation> ior_equation;
  bool is_stop;

  LensElement(const unsigned int index, const  double aperture_radius, const double thickness,
              const double curvature_radius, const std::shared_ptr<IOREquation>& ior_equation,
              const bool is_stop)
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
      double t = -(ray.o.z - z) / ray.d.z;
      const Vector3D hit_pos = ray.at_time(t);
      const double r = hit_pos.x * hit_pos.x + hit_pos.y * hit_pos.y;
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hitPos = ray.at_time(t);
      res.hitNormal = align_normal(ray.d, Vector3D(0, 0, -1));
      return true;
    } else {
      const Vector3D center(0, 0, z + curvature_radius);
      const double b = dot(ray.o - center, ray.d);
      const double c = (ray.o - center).norm2() - curvature_radius * curvature_radius;
      const double discriminant = b * b - c;
      if (discriminant < 0) return false;

      const double t0 = -b - std::sqrt(discriminant);
      const double t1 = -b + std::sqrt(discriminant);
      const double t = curvature_radius * ray.d.z > 0 ? t0 : t1;
      const Vector3D hit_pos = ray.at_time(t);

      const double r = hit_pos.x * hit_pos.x + hit_pos.y * hit_pos.y;
      if (r > aperture_radius * aperture_radius) return false;

      res.t = t;
      res.hitPos = hit_pos;
      res.hitNormal = align_normal(ray.d, (hit_pos - center).unit());
      return true;
    }
  }

  double ior(const double lambda) const { return ior_equation->ior(lambda); }

  GridData<Vector3D> sample_points(const unsigned int N) const {
    GridData<Vector3D> ret(N, N);

    // Grid sampling
    for (unsigned int j = 0; j < N; ++j) {
      const double v = (2.0 * (j + 0.5) - N) / N;
      const double y = v * aperture_radius;
      for (unsigned int i = 0; i < N; ++i) {
        const double u = (2.0 * (i + 0.5) - N) / N;
        const double x = u * aperture_radius;
        ret.set(j, i, Vector3D(x, y, z));
      }
    }

    return ret;
  }
};

}  // namespace CGL

#endif  // CGL_LENS_ELEMENT_H