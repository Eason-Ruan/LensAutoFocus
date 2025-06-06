#ifndef CGL_RAY_H
#define CGL_RAY_H

#include "CGL/CGL.h"
#include "CGL/vector3D.h"
#include "CGL/vector4D.h"
#include "CGL/matrix4x4.h"

#define PART 5

#define PART_1 (PART >= 1)
#define PART_2 (PART >= 2)
#define PART_3 (PART >= 3)
#define PART_4 (PART >= 4)
#define PART_5 (PART >= 5)

namespace CGL {


struct Ray {
  size_t depth;  ///< depth of the Ray
  Vector3D o;  ///< origin
  Vector3D d;  ///< direction
    double lambda;
    double lambda_pdf;
    double ray_pdf;
  mutable double min_t; ///< treat the ray as a segment (ray "begin" at min_t)
  mutable double max_t; ///< treat the ray as a segment (ray "ends" at max_t)

  Vector3D inv_d;  ///< component wise inverse

  Ray() {}

  /**
   * Constructor.
   * Create a ray instance with given origin and direction.
   * \param o origin of the ray
   * \param d direction of the ray
   * \param depth depth of the ray
   */
    Ray(const Vector3D o, const Vector3D d, int depth = 0)
        : o(o), d(d), min_t(0.0), max_t(INF_D), depth(depth), lambda(0) {
    inv_d = 1.0 / d;
  }

  /**
   * Constructor.
   * Create a ray instance with given origin and direction.
   * \param o origin of the ray
   * \param d direction of the ray
   * \param max_t max t value for the ray (if it's actually a segment)
   * \param depth depth of the ray
   */
    Ray(const Vector3D o, const Vector3D d, double max_t, int depth = 0)
        : o(o), d(d), min_t(0.0), max_t(max_t), depth(depth), lambda(0) {
    inv_d = 1.0 / d;
  }


  /**
   * Returns the point t * |d| along the ray.
   */
  inline Vector3D at_time(double t) const { return o + t * d; }

  /**
   * Returns the result of transforming the ray by the given transformation
   * matrix.
   */
  Ray transform_by(const Matrix4x4& t) const {
    const Vector4D& newO = t * Vector4D(o, 1.0);
    return Ray((newO / newO.w).to3D(), (t * Vector4D(d, 0.0)).to3D());
  }
};

// structure used for logging rays for subsequent visualization
struct LoggedRay {

    LoggedRay(const Ray& r, double hit_t)
        : o(r.o), d(r.d), hit_t(hit_t) {}

    Vector3D o;
    Vector3D d;
    double hit_t;
};

}  // namespace CGL

#endif  // CGL_RAY_H
