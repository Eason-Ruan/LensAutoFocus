#include "lens-system/lens-system.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"
using JSON = nlohmann::json;

#include "CGL/parallel.h"

namespace CGL {

LensSystem::LensSystem(const std::string& filename, const std::shared_ptr<Film>& film)
    : film(film) {
  // Load JSON
  if (!load_json(filename)) exit(EXIT_FAILURE);

  // Compute aperture index
  for (size_t i = 0; i < elements.size(); ++i) {
    if (elements[i].is_stop) {
      aperture_index = i;
      break;
    }
  }

  // Compute system length and z
  double length = 0;
  for (auto itr = elements.rbegin(); itr != elements.rend(); ++itr) {
    length += itr->thickness;
    itr->z = -length;
  }
  system_length = length;

  // Compute cardinal points
  compute_cardinal_points();

  // Compute FOV
  horizontal_fov = 2.0 * std::atan2(film->width_length, 2.0 * image_focal_length);
  vertical_fov = 2.0 * std::atan2(film->height_length, 2.0 * image_focal_length);
  diagonal_fov = 2.0 * std::atan2(film->diagonal_length, 2.0 * image_focal_length);
}

bool LensSystem::load_json(const std::string& filename) {
  // Open file
  std::ifstream stream(filename);
  if (!stream) {
    std::cerr << "Failed to open " << filename << std::endl;
    return false;
  }

  // Parse JSON
  JSON json;
  stream >> json;
  stream.close();

  // Push lens elements
  for (const auto& [key, value] : json.items()) {
    // Required
    const unsigned int index = value["index"].get<unsigned int>();
    const double curvature_radius = value["curvature_radius"].get<double>();
    const double thickness = value["thickness"].get<double>();
    const double aperture_radius = 0.5 * value["aperture_diameter"].get<double>();

    // Optional
    bool is_stop = value.value("is_stop", false);

    std::optional<double> nd = value.contains("nd") ? value["nd"].get<double>() : std::nullopt;
    std::optional<double> nD = value.contains("nD") ? value["nD"].get<double>() : std::nullopt;
    std::optional<double> nF = value.contains("nF") ? value["nF"].get<double>() : std::nullopt;

    // Select IOR equation
    std::shared_ptr<IOREquation> ior_equation;
    if (nd) {
      ior_equation = std::make_shared<ConstantIOR>(nd.value());
    } else if (nD && nF) {
      ior_equation = std::make_shared<CauchyEquation>(fit_cauchy(nD.value(), nF.value()));
    } else {
      std::cerr << "Failed to create IOR equation for this lens element" << std::endl;
      return false;
    }

    // Make lens element
    elements.emplace_back(index, aperture_radius, thickness, curvature_radius, ior_equation, is_stop);
  }

  // Sort lens elements by index
  std::sort(elements.begin(), elements.end(),
            [](const LensElement& x1, const LensElement& x2) {
              return x1.index < x2.index;
            });

  return true;
}

double LensSystem::effective_focal_length() const { return image_focal_length; }
double LensSystem::front_focal_length() const {
  return std::abs(elements.front().z - object_focal_z);
}
double LensSystem::back_focal_length() const {
  return std::abs(image_focal_z - elements.back().z);
}

bool LensSystem::raytrace(const Ray& ray_in, Ray& ray_out, bool reflection, Sampler* sampler) const {
  int element_index = ray_in.direction.z > 0 ? -1 : elements.size();
  const int initial_element_index = element_index;

  Ray ray = ray_in;
  double ior = 1.0;

  while (true) {
    // Update element index
    element_index += ray.direction.z > 0 ? 1 : -1;
    if (element_index < 0 || element_index >= elements.size()) break;
    const LensElement& element = elements[element_index];

    // Aperture
    if (element.is_stop) {
      Hit res;
      if (!element.intersect(ray, res)) return false;

      // Update ray
      ray.origin = res.hit_pos;

      // Update IOR
      ior = 1.0;
    }
    // Lens
    else {
      // Compute next element IOR
      double next_ior = 1.0;

      // Compute next element
      const int next_element_index = ray.direction.z > 0 ? element_index : element_index - 1;
      if (next_element_index >= 0) {
        const LensElement& next_element = elements[next_element_index];
        if (!next_element.is_stop) {
          next_ior = next_element.ior(ray.lambda);
        }
      }

      // Compute intersection with lens
      Hit res;
      if (!element.intersect(ray, res)) return false;

      // Refract and reflect
      Vector3D next_direction;
      if (reflection) {
        const double fr = fresnel(-ray.direction, res.hit_normal, ior, next_ior);
        if (sampler->get_next() < fr) {
          // Reflection
          next_direction = reflect(-ray.direction, res.hit_normal);
        } else {
          // Refract
          if (!refract(-ray.direction, next_direction, res.hit_normal, ior, next_ior)) {
            // Total reflection
            next_direction = reflect(-ray.direction, res.hit_normal);
          }
        }
      } else {
        if (!refract(-ray.direction, next_direction, res.hit_normal, ior, next_ior))
          return false;
      }

      // Set next ray
      ray.origin = res.hit_pos;
      ray.direction = next_direction.unit();

      // Update IOR
      ior = next_ior;
    }
  }

  // If ray exits from the same side
  if (element_index == initial_element_index) return false;

  ray_out = ray;

  return true;
}


std::vector<Ray> LensSystem::raytrace_path(const Ray& ray_in) const {
  std::vector<Ray> result;

  int element_index = ray_in.direction.z > 0 ? -1 : elements.size();
  const int initial_element_index = element_index;

  Ray ray = ray_in;
  double ior = 1.0;
  while (true) {
    // Push ray origin
    result.push_back(ray);

    // Update element index
    element_index += ray.direction.z > 0 ? 1 : -1;
    if (element_index < 0 || element_index >= elements.size()) break;
    const LensElement& element = elements[element_index];

    // Aperture
    if (element.is_stop) {
      Hit res;
      if (!element.intersect(ray, res)) break;

      // Update ray
      ray.origin = res.hit_pos;

      // Update IOR
      ior = 1.0;
    }
    // Lens
    else {
      // Compute next element IOR
      double next_ior = 1.0;

      // Compute next element
      const int next_element_index =
          ray.direction.z > 0 ? element_index : element_index - 1;
      if (next_element_index >= 0) {
        const LensElement& next_element = elements[next_element_index];
        if (!next_element.is_stop) {
          next_ior = next_element.ior(ray.lambda);
        }
      }

      // Compute intersection with lens
      Hit res;
      if (!element.intersect(ray, res)) break;

      // Refract
      Vector3D next_direction;
      if (!refract(-ray.direction, next_direction, res.hit_normal, ior, next_ior))
        break;

      // Set next ray
      ray.origin = res.hit_pos;
      ray.direction = next_direction.unit();

      // Update IOR
      ior = next_ior;
    }
  }

  return result;
}

std::vector<ParaxialRay> LensSystem::raytrace_paraxial(const ParaxialRay& ray_in,
                                                       int start, int end,
                                                       double lambda) const {
  std::vector<ParaxialRay> result;
  result.push_back(ray_in);

  // Compute start and end indices
  int start_index = start == -1 ? elements.size() - 1 : start;
  int end_index = end == -1 ? elements.size() - 1 : end;

  // Invalid index case
  if (start_index < 0 || start_index >= elements.size() || end_index < 0 ||
      end_index >= elements.size()) {
    std::cerr << "Invalid start or end index" << std::endl;
    return result;
  }

  // Paraxial raytrace
  double ior;
  double ior_prev = (start_index == 0 || start_index == elements.size() - 1)
                        ? 1.0
                        : elements[start_index].ior(lambda);
  double u_prev = ray_in.u;
  double h_prev = ray_in.h;

  if (start_index < end_index) {
    for (int i = start_index; i <= end_index; ++i) {
      const auto& element = elements[i];

      // Compute curvature radius
      double r = element.is_stop ? 1e9 : element.curvature_radius;

      // Compute IOR of element
      ior = element.ior(lambda);

      // Compute thickness
      double thickness = (i != elements.size() - 1) ? element.thickness : 0;

      // Compute paraxial ray
      const double u = ior_prev / ior * u_prev + (ior - ior_prev) / (ior * r) * h_prev;
      const double h = h_prev - thickness * u;

      // Save paraxial ray
      result.push_back(ParaxialRay(u, h));

      // Update
      ior_prev = ior;
      u_prev = u;
      h_prev = h;
    }
  } else {
    for (int i = start_index; i >= end_index; --i) {
      // Compute curvature radius
      double r = elements[i].is_stop ? -1e9 : -elements[i].curvature_radius;

      // Compute IOR of element
      ior = (i > 0) ? elements[i - 1].ior(lambda) : 1.0;

      // Compute thickness
      double thickness = (i > 0) ? elements[i - 1].thickness : 0;

      // Compute paraxial ray
      const double u = ior_prev / ior * u_prev + (ior - ior_prev) / (ior * r) * h_prev;
      const double h = h_prev - thickness * u;

      // Save paraxial ray
      result.push_back(ParaxialRay(u, h));

      // Update
      ior_prev = ior;
      u_prev = u;
      h_prev = h;
    }
  }

  return result;
}

void LensSystem::compute_cardinal_points() {
  // Compute image focal point
  // Paraxial raytrace with (u, h) = (0, 1)
  auto result = raytrace_paraxial(ParaxialRay(0, 1));
  image_focal_z = elements.back().z + result.back().h / result.back().u;

  // Compute principal point
  image_principal_z = elements.back().z +
                      (result.back().h - result.front().h) / result.back().u;

  // Compute image focal length
  image_focal_length = image_focal_z - image_principal_z;

  // Compute object focal point
  // Paraxial reverse raytrace with (u, h) = (0, 1)
  result = raytrace_paraxial(ParaxialRay(0, 1), -1, 0);
  object_focal_z = elements.front().z - result.back().h / result.back().u;

  // Compute object principal point
  object_principal_z = elements.front().z -
                       (result.back().h - result.front().h) / result.back().u;

  // Compute object focal length
  object_focal_length = -(object_focal_z - object_principal_z);
}

bool LensSystem::focus(double focus_z) {
  const double delta =
      0.5 * (object_principal_z - focus_z + image_principal_z -
             std::sqrt((object_principal_z - focus_z - image_principal_z) *
                       (object_principal_z - focus_z - 4 * image_focal_length -
                        image_principal_z)));

  // Move lens elements
  for (auto& element : elements) {
    element.z -= delta;
  }

  // Recompute cardinal points
  compute_cardinal_points();

  return true;
}

Bounds2 LensSystem::compute_exit_pupil_bound(const Vector2D& p) const {
  Bounds2 bounds;

  const auto& last_element = elements.back();
  Ray ray_out;
  for (int i = 0; i < num_exit_pupil_bounds_samples; ++i) {
    for (int j = 0; j < num_exit_pupil_bounds_samples; ++j) {
      // Sample point on last element surface
      const double u =
          2.0 * static_cast<double>(i) / num_exit_pupil_bounds_samples - 1.0;
      const double v =
          2.0 * static_cast<double>(j) / num_exit_pupil_bounds_samples - 1.0;
      const Vector3D sample_point =
          Vector3D(last_element.aperture_radius * u,
                   last_element.aperture_radius * v, last_element.z);

      // Make ray
      const Vector3D origin(p.x, p.y, 0);
      const Ray ray_in(origin, (sample_point - origin).unit());

      // Raytrace
      if (!raytrace(ray_in, ray_out)) continue;

      // Extend bounding box
      bounds = extend_bounds(bounds, Vector2D(sample_point.x, sample_point.y));
    }
  }

  return bounds;
}

bool LensSystem::compute_exit_pupil_bounds() {
  exit_pupil_bounds.resize(num_exit_pupil_bounds);

  Parallel parallel;

  parallel.parallel_for_1d(
      [&](unsigned int idx) {
        const double r = static_cast<double>(idx) / num_exit_pupil_bounds * 0.5 *
                         film->diagonal_length;
        exit_pupil_bounds[idx] = compute_exit_pupil_bound(Vector2D(r, 0));

        std::cout << "Finished " << idx << "th computation of exit pupil bounds"
                  << std::endl;
        std::cout << exit_pupil_bounds[idx] << std::endl;
      },
      16, num_exit_pupil_bounds);

  return true;
}

bool LensSystem::sample_ray(double u, double v, double lambda, Sampler& sampler,
                            Ray& ray_out, double& pdf, bool reflection) const {
  // Compute position on film
  const Vector2D p = film->compute_position(u, v);

  // Choose exit pupil bound
  const double r = p.norm();
  const unsigned int exit_pupil_bounds_index =
      r / (0.5 * film->diagonal_length) * num_exit_pupil_bounds;
  const Bounds2& exit_pupil_bound = exit_pupil_bounds[exit_pupil_bounds_index];
  if (!exit_pupil_bound.is_valid()) return false;

  // Sample point on exit pupil bound
  double pdf_area;
  Vector2D p_bound = exit_pupil_bound.sample_point(sampler, pdf_area);

  // Rotate sampled point
  if (r > 0) {
    const double theta = std::atan2(v, u);
    p_bound = rotate_2d(p_bound, theta);
  }

  // Make input ray
  const Vector3D origin = Vector3D(p.x, p.y, 0);
  const Vector3D p_bound_3d = Vector3D(p_bound.x, p_bound.y, elements.back().z);
  const Vector3D direction = (p_bound_3d - origin).unit();
  const Ray ray_in(origin, direction, lambda);

  // Convert area pdf to solid angle pdf
  const double l = (p_bound_3d - origin).norm();
  pdf = l * l / std::abs(dot(direction, Vector3D(0, 0, -1))) * pdf_area;

  // Raytrace
  if (!raytrace(ray_in, ray_out, reflection, &sampler)) return false;

  return true;
}

GridData<std::pair<bool, Ray>> LensSystem::raytrace_n(
  const GridData<Ray>& rays_in, bool reflection, Sampler* sampler) const {
GridData<std::pair<bool, Ray>> result(rays_in.nrows, rays_in.ncols);

Parallel parallel;
parallel.parallel_for_2d(
    [&](unsigned int i, unsigned int j) {
      bool traced;
      Ray ray_out;
      traced = raytrace(rays_in.get(i, j), ray_out, reflection, sampler);

      result.set(i, j, {traced, ray_out});
    },
    16, 16, rays_in.nrows, rays_in.ncols);

return result;
}

std::pair<GridData<double>, std::array<double, 4>> LensSystem::compute_exit_pupil(
  const Vector2D& p_film, unsigned int n_grids) const {
const auto& last_element = elements.back();

// Compute extents
std::array<double, 4> extents;
extents[0] = -last_element.aperture_radius;
extents[1] = last_element.aperture_radius;
extents[2] = -last_element.aperture_radius;
extents[3] = last_element.aperture_radius;

// Compute grids
const GridData<Vector3D> grids = elements.back().sample_points(n_grids);

// Make rays
GridData<Ray> rays_in(n_grids, n_grids);
for (unsigned int i = 0; i < n_grids; ++i) {
  for (unsigned int j = 0; j < n_grids; ++j) {
    const Vector3D p_film_3d = Vector3D(p_film.x, p_film.y, 0);
    rays_in.set(i, j, Ray(p_film_3d, (grids.get(i, j) - p_film_3d).unit()));
  }
}

// Raytrace
const auto result = raytrace_n(rays_in);

// Represent exit pupil as 0, 1
GridData<double> exit_pupil(result.nrows, result.ncols);
for (unsigned int i = 0; i < result.nrows; ++i) {
  for (unsigned int j = 0; j < result.ncols; ++j) {
    exit_pupil.set(i, j, result.get(i, j).first);
  }
}

return {exit_pupil, extents};
}

bool LensSystem::compute_primary_ray(const Vector3D& origin, Ray& primary_ray,
                                   unsigned int n_grids) const {
// Compute grids
const GridData<Vector3D> grids = elements.front().sample_points(n_grids);

// Make rays
GridData<Ray> rays_in(n_grids, n_grids);
for (unsigned int i = 0; i < n_grids; ++i) {
  for (unsigned int j = 0; j < n_grids; ++j) {
    rays_in.set(i, j, Ray(origin, (grids.get(i, j) - origin).unit()));
  }
}

// Raytrace
const auto result = raytrace_n(rays_in);

// Compute entrance pupil center
unsigned int n_average = 0;
Vector3D entrance_pupil_center;
for (unsigned int i = 0; i < n_grids; ++i) {
  for (unsigned int j = 0; j < n_grids; ++j) {
    if (result.get(i, j).first) {
      entrance_pupil_center += grids.get(i, j);
      n_average++;
    }
  }
}
if (n_average == 0) {
  std::cerr << "Failed to compute primary ray" << std::endl;
  return false;
}
entrance_pupil_center /= n_average;

primary_ray = Ray(origin, (entrance_pupil_center - origin).unit());

return true;
}

std::vector<Vector3D> LensSystem::compute_spot_diagram(const Vector3D& origin,
                                                     unsigned int n_grids) const {
std::vector<Vector3D> result;

// Compute grids
const GridData<Vector3D> grids = elements.front().sample_points(n_grids);

// Make rays
GridData<Ray> rays_in(n_grids, n_grids);
for (unsigned int i = 0; i < n_grids; ++i) {
  for (unsigned int j = 0; j < n_grids; ++j) {
    rays_in.set(i, j, Ray(origin, (grids.get(i, j) - origin).unit()));
  }
}

// Raytrace
const auto raytrace_result = raytrace_n(rays_in);

// Compute intersect position at Gaussian plane
for (unsigned int i = 0; i < n_grids; ++i) {
  for (unsigned int j = 0; j < n_grids; ++j) {
    if (raytrace_result.get(i, j).first) {
      const Ray& ray = raytrace_result.get(i, j).second;
      const double t = -(ray.origin.z - image_focal_z) / ray.direction.z;
      const Vector3D p_film = ray.at(t);
      result.push_back(p_film);
    }
  }
}

return result;
}

std::pair<GridData<double>, std::array<double, 4>> LensSystem::compute_geometric_psf(
  const Vector3D& origin, unsigned int n_rays, unsigned int n_grids) const {
// Compute grids
const GridData<Vector3D> grids = elements.front().sample_points(n_rays);

// Make rays
GridData<Ray> rays_in(n_rays, n_rays);
for (unsigned int i = 0; i < n_rays; ++i) {
  for (unsigned int j = 0; j < n_rays; ++j) {
    rays_in.set(i, j, Ray(origin, (grids.get(i, j) - origin).unit()));
  }
}

// Raytrace
const auto raytrace_result = raytrace_n(rays_in);

// Compute intersect position at Gaussian plane
std::vector<Vector3D> points;
for (unsigned int i = 0; i < n_rays; ++i) {
  for (unsigned int j = 0; j < n_rays; ++j) {
    if (raytrace_result.get(i, j).first) {
      const Ray& ray = raytrace_result.get(i, j).second;
      const double t = -(ray.origin.z - image_focal_z) / ray.direction.z;
      const Vector3D p_film = ray.at(t);
      points.push_back(p_film);
    }
  }
}

// Compute primary ray
Ray primary_ray;
if (!compute_primary_ray(origin, primary_ray)) {
  std::cerr << "Failed to compute primary ray" << std::endl;
  std::exit(EXIT_FAILURE);
}

// Raytrace
Ray primary_ray_out;
if (!raytrace(primary_ray, primary_ray_out)) {
  std::cerr << "Failed to raytrace primary ray" << std::endl;
  std::exit(EXIT_FAILURE);
}

// Compute intersect position at Gaussian plane
const Vector3D p_film_primary =
    primary_ray_out.at(-(primary_ray_out.origin.z - image_focal_z) /
                       primary_ray_out.direction.z);

// Compute mean and variance
const Vector3D p_mean =
    std::accumulate(points.begin(), points.end(), Vector3D()) / points.size();
const Vector3D p_var =
    std::accumulate(points.begin(), points.end(), Vector3D(),
                    [&](const Vector3D& v1, const Vector3D& v2) {
                      return v1 + (v2 - p_mean) * (v2 - p_mean);
                    }) /
    points.size();

// Compute extent
const double x_std = std::sqrt(p_var.x);
const double y_std = std::sqrt(p_var.y);
const double grid_width = x_std > y_std ? 3 * x_std : 3 * y_std;
const double xmin = p_film_primary.x - grid_width;
const double xmax = p_film_primary.x + grid_width;
const double ymin = p_film_primary.y - grid_width;
const double ymax = p_film_primary.y + grid_width;
const std::array<double, 4> extent = {xmin, xmax, ymin, ymax};

// Initialize PSF grid
GridData<double> psf(n_grids, n_grids);
for (unsigned int i = 0; i < n_grids; ++i) {
  for (unsigned int j = 0; j < n_grids; ++j) {
    psf.set(i, j, 0);
  }
}

// Compute PSF by grid
for (const auto& p : points) {
  // Compute grid index
  const unsigned int i = (p.x - xmin) / (xmax - xmin) * n_grids;
  const unsigned int j = (p.y - ymin) / (ymax - ymin) * n_grids;

  // Accumulate
  if (i < n_grids && j < n_grids) {
    psf.set(j, i, psf.get(j, i) + 1);
  }
}

// Normalize
const double max_value = *std::max_element(psf.data.begin(), psf.data.end());
for (auto& value : psf.data) {
  value /= max_value;
}

return {psf, extent};
}
}  // namespace CGL