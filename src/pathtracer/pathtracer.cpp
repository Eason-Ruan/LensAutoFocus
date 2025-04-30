#include "pathtracer.h"

#include "type.h"
#include "CGL/spectrum.h"
#include "scene/light.h"
#include "scene/sphere.h"
#include "misc.h"


using namespace CGL::SceneObjects;

namespace CGL {
  Real getSpectrum(const Vector3D& radiance,const Ray& ray) {
    return RGB2Spectrum(radiance).sample(ray.lambda);
  }

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  const int num_samples = scene->lights.size() * ns_area_light;
  auto L_out = Vector3D(0, 0, 0);

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading
  for (int i = 0; i < num_samples; i++) {
    Vector3D w_in, L_in;
    double pdf;
    Intersection bounce_intersect;
    Vector3D sample_f = isect.bsdf->sample_f(w_out, &w_in, &pdf);
    auto r_in = Ray(hit_p, o2w * w_in);
    r_in.min_t = EPS_F;
    if (!bvh->intersect(r_in, &bounce_intersect)) {
      continue;
    }
    L_in = zero_bounce_radiance(r_in, bounce_intersect);
    L_out += sample_f * L_in * cos_theta(w_in) / (pdf * num_samples);
  }
  return L_out;
}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  auto L_out = Vector3D(0, 0, 0);
  for (const auto l: scene->lights) {
    if (l->is_delta_light()) {
      Vector3D w_in_world, w_in_o, L_in;
      double pdf, dis;
      L_in = l->sample_L(hit_p, &w_in_world, &dis, &pdf);
      if (dot(w_in_world, isect.n) < - EPS_F) {
        continue;
      }
      w_in_o = w2o * w_in_world;
      auto r_in = Ray(hit_p, w_in_world);
      r_in.min_t = EPS_F;
      Intersection bounce_intersect;
      bvh->intersect(r_in, &bounce_intersect);
      double bounce_dis = bounce_intersect.t * w_in_world.norm();
      if (bounce_dis - dis > EPS_F || dis - bounce_dis > EPS_F) {
        continue;
      }
      Vector3D f = isect.bsdf->f(w_out, w_in_o);
      L_out += f * L_in * cos_theta(w_in_o) / pdf;
    }else {
      for (int i = 0; i < ns_area_light; i++) {
        Vector3D w_in_world, w_in_o, L_in;
        double pdf, dis;
        L_in = l->sample_L(hit_p, &w_in_world, &dis, &pdf);
        if (dot(w_in_world, isect.n) < - EPS_F) {
          continue;
        }
        w_in_o = w2o * w_in_world;
        auto r_in = Ray(hit_p, w_in_world);
        r_in.min_t = EPS_F;
        Intersection bounce_intersect;
        bvh->intersect(r_in, &bounce_intersect);
        double bounce_dis = bounce_intersect.t * w_in_world.norm();
        if (bounce_dis - dis > EPS_F || dis - bounce_dis > EPS_F) {
          continue;
        }
        Vector3D f = isect.bsdf->f(w_out, w_in_o);
        L_out += f * L_in * cos_theta(w_in_o) / (pdf * static_cast<double>(ns_area_light));
      }
    }
  }
  return L_out;
}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light
  return isect.bsdf->get_emission();


}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
  return direct_hemisphere_sample ? estimate_direct_lighting_hemisphere(r, isect) : estimate_direct_lighting_importance(r, isect);

  return Vector3D(1.0);


}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out = isAccumBounces ? one_bounce_radiance(r, isect) : Vector3D(0, 0, 0);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.
  if (coin_flip(0.4) || r.depth >= max_ray_depth - 1) {
    if (isAccumBounces) {
      return one_bounce_radiance(r, isect);
    }
    return L_out;
  }
  Vector3D w_in;
  double pdf;
  Intersection bounce_intersect;
  Vector3D sample_f = isect.bsdf->sample_f(w_out, &w_in, &pdf);
  auto r_in = Ray(hit_p, o2w * w_in, static_cast<int>(r.depth + 1));
  r_in.min_t = EPS_F;
  if (!bvh->intersect(r_in, &bounce_intersect)) {
    return isAccumBounces ? L_out : Vector3D(0, 0, 0);
  }
  L_out += sample_f * at_least_one_bounce_radiance(r_in, bounce_intersect) * cos_theta(w_in) / (pdf * 0.4);
  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  if (!bvh->intersect(r, &isect))
    return Vector3D(0, 0, 0);
  return zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect);
  // TODO (Part 3): Return the direct illumination.

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct

  return L_out;
}

void PathTracer::raytrace_pixel(const size_t x, const size_t y, const bool is_focusing) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.
  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"
  const int num_samples = static_cast<int>(ns_aa);          // total samples to evaluate
  const auto origin = Vector2D(static_cast<int>(x), static_cast<int>(y)); // bottom left corner of the pixel
  auto estSample = Vector3D(0.0, 0.0, 0.0);
  float s1 = 0.0;
  float s2 = 0.0;
  int sampled_num = 0;
  for (int i = 0; i < num_samples; i ++) {
    const auto samplePoint = origin + gridSampler->get_sample();
    Ray ray_out = camera->generate_ray(samplePoint.x / static_cast<double>(sampleBuffer.w),
                           samplePoint.y / static_cast<double>(sampleBuffer.h));
    if (ray_out.d == Vector3D(0, 0, 0)) {
      continue; // skip if Ray failed to pass lens,
    }
    Vector3D sample = est_radiance_global_illumination(ray_out);
    const Real cos  = std::abs(dot(ray_out.d, camera->look_direction()));
    if (spectrumSampling) {
      sample = srgbToLinear(sample / 255.);
      const Real radiance = getSpectrum(sample, ray_out) * cos / (ray_out.ray_pdf * ray_out.lambda_pdf);
      const int index = static_cast<int>((ray_out.lambda - 380) / 5);
      Vector3D XYZ;
      if (index >= 0 && index <= SPD::color_matching_func_samples - 1) {
        XYZ.x = radiance * SPD::color_matching_func_x[index];
        XYZ.y = radiance * SPD::color_matching_func_y[index];
        XYZ.z = radiance * SPD::color_matching_func_z[index];
      }
      sample = XYZ;
    } else {
      sample = sample * cos * gain / (ray_out.ray_pdf);
    }
    estSample += sample;
    sampled_num ++;
    const float sample_illum = sample.illum();
    s1 += sample_illum;
    s2 += sample_illum * sample_illum;
    if (sampled_num % samplesPerBatch == 0 && sampled_num != 0) {
      if (1.96 * sqrt((s2 - s1 * s1 / static_cast<float>(sampled_num)) / (static_cast<float>(sampled_num) - 1)) / sqrt(sampled_num) <= maxTolerance * s1 / sampled_num) {
        break;
      }
    }
  }
  if (spectrumSampling) {
    estSample /= sampled_num;
    estSample = XYZ2RGB(estSample);
    estSample = linearToSrgb(estSample);
  } else {
    estSample /= sampled_num;
  }
  sampleBuffer.update_pixel(estSample, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = sampled_num;
}


void PathTracer::focus(double delta_distance) {
  // TODO Redesign the autofocus function
}

} // namespace CGL
