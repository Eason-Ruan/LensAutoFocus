#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include "core/spectrum.h"

namespace CGL {

// Construct SPD from non-uniform spectral power distribution (SPD)
// Interpolate non-uniform SPD to create uniform SPD
SPD::SPD(const std::vector<double>& _lambda, const std::vector<double>& _phi) {
  assert(_lambda.size() == _phi.size());  // Ensure wavelength and power arrays match

  // If only one sample
  if (_lambda.size() == 1) {
    // Initialize other wavelengths' power to 0
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] = 0;
    }
    add_phi(_lambda[0], _phi[0]);  // Add power for the single sample
  }
  // If multiple samples
  else {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      // Compute uniform wavelength
      const double lambda_value = LAMBDA_MIN + LAMBDA_INTERVAL * i;

      // If wavelength is out of range, set power to 0
      if (lambda_value < _lambda.front() || lambda_value > _lambda.back()) {
        phi[i] = 0;
      } else if (lambda_value == _lambda.front()) {
        phi[i] = _phi.front();  // If wavelength equals minimum, use corresponding power
      } else if (lambda_value == _lambda.back()) {
        phi[i] = _phi.back();  // If wavelength equals maximum, use corresponding power
      } else {
        // Interpolate power for non-uniform SPD
        const size_t lambda1_index = static_cast<size_t>(
            std::lower_bound(_lambda.begin(), _lambda.end(), lambda_value) -
            _lambda.begin());  // Find index of the first wavelength >= current wavelength
        const size_t lambda0_index = lambda1_index - 1;  // Previous index
        assert(lambda0_index != lambda1_index);

        // Compute interpolation factor
        const double t = (lambda_value - _lambda[lambda0_index]) /
                         (_lambda[lambda1_index] - _lambda[lambda0_index]);
        assert(t >= 0 && t <= 1);

        // Interpolate power
        const double interpolated_phi =
            (1.0 - t) * _phi[lambda0_index] + t * _phi[lambda1_index];
        phi[i] = interpolated_phi;
      }
    }
  }
}

// Add power to a specific wavelength
void SPD::add_phi(const double& _lambda, const double& _phi) {
  // If wavelength is out of range, do nothing
  if (_lambda < LAMBDA_MIN || _lambda >= LAMBDA_MAX) {
    return;
  }

  // Compute index for the wavelength
  const size_t lambda_index = (_lambda - LAMBDA_MIN) / LAMBDA_INTERVAL;

  // Distribute power to adjacent wavelengths
  const double lambda0 =
      LAMBDA_MIN + lambda_index * LAMBDA_INTERVAL;  // Left wavelength
  const double t =
      (_lambda - lambda0) / LAMBDA_INTERVAL;  // Position within the interval
  phi[lambda_index] += (1.0 - t) * _phi;  // Assign to left wavelength
  phi[lambda_index + 1] += t * _phi;      // Assign to right wavelength
}

// Sample power at a specific wavelength
double SPD::sample(const double& l) const {
  assert(l >= LAMBDA_MIN && l < LAMBDA_MAX);  // Ensure wavelength is within range

  // Compute index for the wavelength
  const size_t lambda_index = (l - LAMBDA_MIN) / LAMBDA_INTERVAL;

  // Interpolate power
  if (lambda_index == 0 || lambda_index == LAMBDA_SAMPLES - 1) {
    return phi[lambda_index];  // If at boundary, return power directly
  } else {
    const double lambda_nearest = LAMBDA_MIN + lambda_index * LAMBDA_INTERVAL;
    const double t = (l - lambda_nearest) / LAMBDA_INTERVAL;  // Interpolation factor
    assert(t >= 0 && t <= 1);
    return (1.0 - t) * phi[lambda_index] + t * phi[lambda_index + 1];
  }
}

// Convert SPD to XYZ color space
XYZ SPD::to_xyz() const {
  XYZ xyz;

  for (std::size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
    const double lambda_value = LAMBDA_MIN + LAMBDA_INTERVAL * i;  // Current wavelength
    const double phi_value = phi[i];  // Power at current wavelength

    // Skip if power is 0
    if (phi_value == 0) continue;

    // Compute index for color matching function
    const int index = (lambda_value - 380) / 5;
    assert(index >= 0 && index < color_matching_func_samples);

    // Interpolate color matching function
    double cmf_x, cmf_y, cmf_z;
    if (index != color_matching_func_samples - 1) {
      const double cmf_lambda = 5 * index + 380;
      const double t = (lambda_value - cmf_lambda) / 5;
      assert(t >= 0 && t <= 1);

      cmf_x = (1.0 - t) * color_matching_func_x[index] +
              t * color_matching_func_x[index + 1];
      cmf_y = (1.0 - t) * color_matching_func_y[index] +
              t * color_matching_func_y[index + 1];
      cmf_z = (1.0 - t) * color_matching_func_z[index] +
              t * color_matching_func_z[index + 1];
    } else {
      cmf_x = color_matching_func_x[index];
      cmf_y = color_matching_func_y[index];
      cmf_z = color_matching_func_z[index];
    }

    // Compute XYZ values using trapezoidal approximation
    xyz[0] += cmf_x * phi_value;
    xyz[1] += cmf_y * phi_value;
    xyz[2] += cmf_z * phi_value;
  }

  return xyz;
}

// Convert RGB color to spectral representation
SPD rgb_to_spectrum(const RGB& rgb) {
  // Sampled wavelengths
  static const std::vector<double> sampled_lambda = {
      380, 417.7, 455.55, 493.33, 531.11, 568.88, 606.66, 644.44, 682.22, 720};

  // Define base spectra
  static const SPD white_spectrum =
      SPD(sampled_lambda, {1, 1, 0.9999, 0.9993, 0.9992, 0.9998, 1, 1, 1, 1});
  static const SPD cyan_spectrum =
      SPD(sampled_lambda,
          {0.9710, 0.9426, 1.0007, 1.0007, 1.0007, 1.0007, 0.1564, 0, 0, 0});
  static const SPD magenta_spectrum = SPD(
      sampled_lambda, {1, 1, 0.9685, 0.2229, 0, 0.0458, 0.8369, 1, 1, 0.9959});
  static const SPD yellow_spectrum =
      SPD(sampled_lambda,
          {0.0001, 0, 0.1088, 0.6651, 1, 1, 0.9996, 0.9586, 0.9685, 0.9840});
  static const SPD red_spectrum =
      SPD(sampled_lambda,
          {0.1012, 0.0515, 0, 0, 0, 0, 0.8325, 1.0149, 1.0149, 1.0149});
  static const SPD green_spectrum = SPD(
      sampled_lambda, {0, 0, 0.0273, 0.7937, 1, 0.9418, 0.1719, 0, 0, 0.0025});
  static const SPD blue_spectrum =
      SPD(sampled_lambda,
          {1, 1, 0.8916, 0.3323, 0, 0, 0.0003, 0.0369, 0.0483, 0.0496});

  SPD result;
  const double red = rgb.x();   // Extract red component
  const double green = rgb.y(); // Extract green component
  const double blue = rgb.z();  // Extract blue component

  // Combine base spectra based on RGB components
  if (red <= green && red <= blue) {
    result += red * white_spectrum;
    if (green <= blue) {
      result += (green - red) * cyan_spectrum;
      result += (blue - green) * blue_spectrum;
    } else {
      result += (blue - red) * cyan_spectrum;
      result += (green - blue) * green_spectrum;
    }
  } else if (green <= red && green <= blue) {
    result += green * white_spectrum;
    if (red <= blue) {
      result += (red - green) * magenta_spectrum;
      result += (blue - red) * blue_spectrum;
    } else {
      result += (blue - green) * magenta_spectrum;
      result += (red - blue) * red_spectrum;
    }
  } else {
    result += blue * white_spectrum;
    if (red <= green) {
      result += (red - blue) * yellow_spectrum;
      result += (green - red) * green_spectrum;
    } else {
      result += (green - blue) * yellow_spectrum;
      result += (red - green) * red_spectrum;
    }
  }

  return result;
}

}  // namespace CGL