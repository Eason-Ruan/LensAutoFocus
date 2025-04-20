#ifndef LENS_SYSTEM_H
#define LENS_SYSTEM_H

#include <array>
#include <vector>
#include <cmath>

#include "CGL/bound2.h"
#include "spectral-ray.h"
#include "grid_data.h"
#include "lens-system/lens-element.h"
#include "lens-sampler/sampler.h"
namespace CGL {

// Reflect a vector
inline Vector3D reflect(const Vector3D& v, const Vector3D& n) {
  return -v + 2 * dot(v, n) * n;
}

// Compute Fresnel coefficient
inline double fresnel(const Vector3D& wo, const Vector3D& n, const double n1, const double n2) {
  const double f0 = std::pow((n1 - n2) / (n1 + n2), 2.0);
  return f0 + (1.0 - f0) * std::pow(1.0 - dot(wo, n), 5.0);
}

// Compute refracted vector
inline bool refract(const Vector3D& wi, Vector3D& wt, const Vector3D& n, const double ior1, const double ior2) {
  const double eta = ior1 / ior2;
  const double cos_theta_i = dot(wi, n);
  const double sin2_theta_i = std::max(0.0, 1.0 - cos_theta_i * cos_theta_i);
  const double sin2_theta_t = eta * eta * sin2_theta_i;
  if (sin2_theta_t >= 1.0) return false;
  const double cos_theta_t = std::sqrt(1.0 - sin2_theta_t);
  wt = eta * (-wi) + (eta * cos_theta_i - cos_theta_t) * n;
  return true;
}

// Rotate a point around the origin in 2D
inline Vector2D rotate_2d(const Vector2D& p, double theta) {
  return Vector2D(p.x * std::cos(theta) - p.y * std::sin(theta),
                  p.x * std::sin(theta) + p.y * std::cos(theta));
}

// Paraxial ray for ray tracing
struct ParaxialRay {
  double u;  // Angle with the z-axis
  double h;  // Height of the ray

  ParaxialRay() : u(0), h(0) {}
  ParaxialRay(const double _u, const double _h) : u(_u), h(_h) {}
};

class LensSystem {
 public:
  std::vector<LensElement> elements;
  unsigned int aperture_index;

  double width;
  double height;  // Size of the image sensor

  double system_length;  // Length of the lens system

  double horizontal_fov;  // Horizontal FOV [radians]
  double vertical_fov;    // Vertical FOV [radians]
  double diagonal_fov;    // Diagonal FOV [radians]

  double object_focal_z;       // Object-side focal point position
  double object_principal_z;   // Object-side principal point position
  double object_focal_length;  // Object-side focal length
  double image_focal_z;        // Image-side focal point position
  double image_principal_z;    // Image-side principal point position
  double image_focal_length;   // Image-side focal length

  static constexpr unsigned int num_exit_pupil_bounds = 64;
  static constexpr unsigned int num_exit_pupil_bounds_samples = 1024;
  std::vector<Bounds2> exit_pupil_bounds;

  LensSystem(const std::string& filename, double _width, double _height);

  // Load lens JSON
  bool load_json(const std::string& filename);

  // Compute effective focal length
  double effective_focal_length() const;
  // Compute front focal length
  double front_focal_length() const;
  // Compute back focal length
  double back_focal_length() const;

  // Raytrace
  bool raytrace(const Ray& ray_in, Ray& ray_out, bool reflection = false,
                Prl2::Sampler* sampler = nullptr) const;
  // Raytrace and return raytraced path
  std::vector<Ray> raytrace_path(const Ray& ray_in) const;
  // Raytrace many rays
  GridData<std::pair<bool, Ray>> raytrace_n(const GridData<Ray>& rays_in,
                                            bool reflection = false,
                                            Prl2::Sampler* sampler = nullptr) const;

  // Paraxial raytrace
  std::vector<ParaxialRay> raytrace_paraxial(const ParaxialRay& ray_in,
                                             int start = 0, int end = -1,
                                             double lambda = 550.0) const;

  // Compute principal, focal points, and focal length
  void compute_cardinal_points();

  // Focus lens at z = focus_z
  bool focus(double focus_z);

  // Compute exit pupil bounds at a given point on the film
  Bounds2 compute_exit_pupil_bound(const Vector2D& p) const;
  // Compute exit pupil bounds
  bool compute_exit_pupil_bounds();

  // Compute exit pupil
  std::pair<GridData<double>, std::array<double, 4>> compute_exit_pupil(
      const Vector2D& p_film, unsigned int n_grids = 512) const;

  // Compute spot diagram
  std::vector<Vector3D> compute_spot_diagram(const Vector3D& origin,
                                             unsigned int n_grids) const;

  // Compute geometric PSF
  std::pair<GridData<double>, std::array<double, 4>> compute_geometric_psf(
      const Vector3D& origin, unsigned int n_rays = 512,
      unsigned int n_grids = 512) const;

  // Compute primary ray
  bool compute_primary_ray(const Vector3D& origin, Ray& primary_ray,
                           unsigned int n_grids = 512) const;

  // Sample ray going from the image sensor to object space
  bool sample_ray(double u, double v, double lambda, Prl2::Sampler& sampler, Ray& ray_out,
                  double& pdf, bool reflection = false) const;
};

}  // namespace CGL

#endif  // LENS_SYSTEM_H