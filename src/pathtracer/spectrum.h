#ifndef CGL_SPECTRUM_H
#define CGL_SPECTRUM_H

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"

namespace CGL {

using XYZ = Vector3D;
using RGB = Vector3D;

// Convert XYZ to sRGB (D65)
// Reference: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
inline RGB XYZ_to_RGB(const XYZ& xyz) {
  return RGB(
      3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
      -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
      0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z);
}

// Represents a Spectral Power Distribution (SPD) sampled at equal intervals
class SPD {
 public:
  // Wavelength range stored in SPD
  static constexpr double LAMBDA_MIN = 380.0;
  static constexpr double LAMBDA_MAX = 780.0;

  // Number of wavelength samples
  static constexpr size_t LAMBDA_SAMPLES = 80;

  // Interval between sampled wavelengths
  static constexpr double LAMBDA_INTERVAL =
      (LAMBDA_MAX - LAMBDA_MIN) / LAMBDA_SAMPLES;

  std::array<double, LAMBDA_SAMPLES> phi;  // Radiant flux

  // Default constructor initializes to zero
  SPD() { phi.fill(0.0); }

  // Initialize with a constant value
  explicit SPD(double value) { phi.fill(value); }

  // Construct SPD from arbitrary wavelength and flux samples
  SPD(const std::vector<double>& wavelengths, const std::vector<double>& fluxes);

  // Access radiant flux at a specific index
  double operator[](size_t i) const {
    assert(i < LAMBDA_SAMPLES);
    return phi[i];
  }

  // Clear all values
  void clear() { phi.fill(0.0); }

  // Add radiant flux at a specific wavelength
  void add_flux(double wavelength, double flux);

  // Check if SPD is black (all zero)
  bool is_black() const {
    return std::all_of(phi.begin(), phi.end(), [](double v) { return v == 0.0; });
  }

  // Linearly interpolate radiant flux at a specific wavelength
  double sample(double wavelength) const;

  // Convert SPD to XYZ color space
  XYZ to_XYZ() const;

  // Convert SPD to sRGB color space
  RGB to_RGB() const {
    XYZ xyz = to_XYZ();
    RGB rgb = XYZ_to_RGB(xyz);
    return clamp(rgb, Vector3D(0), Vector3D(INFINITY));
  }

  // Arithmetic operations
  SPD& operator+=(const SPD& other);
  SPD& operator-=(const SPD& other);
  SPD& operator*=(const SPD& other);
  SPD& operator/=(const SPD& other);
  SPD& operator+=(double scalar);
  SPD& operator-=(double scalar);
  SPD& operator*=(double scalar);
  SPD& operator/=(double scalar);

  // CIE 1931 color matching functions
  static constexpr int COLOR_MATCHING_FUNC_SAMPLES = 85;
  static constexpr double COLOR_MATCHING_FUNC_X[COLOR_MATCHING_FUNC_SAMPLES];
  static constexpr double COLOR_MATCHING_FUNC_Y[COLOR_MATCHING_FUNC_SAMPLES];
  static constexpr double COLOR_MATCHING_FUNC_Z[COLOR_MATCHING_FUNC_SAMPLES];
};

// Arithmetic operators for SPD
SPD operator+(const SPD& lhs, const SPD& rhs);
SPD operator-(const SPD& lhs, const SPD& rhs);
SPD operator*(const SPD& lhs, const SPD& rhs);
SPD operator/(const SPD& lhs, const SPD& rhs);
SPD operator+(const SPD& spd, double scalar);
SPD operator-(const SPD& spd, double scalar);
SPD operator*(const SPD& spd, double scalar);
SPD operator/(const SPD& spd, double scalar);
SPD operator+(double scalar, const SPD& spd);
SPD operator-(double scalar, const SPD& spd);
SPD operator*(double scalar, const SPD& spd);
SPD operator/(double scalar, const SPD& spd);

// Normalize SPD
SPD normalize(const SPD& spd);

// Output SPD to stream
std::ostream& operator<<(std::ostream& os, const SPD& spd);

// Convert RGB to SPD (Smits 2001)
SPD RGB_to_Spectrum(const RGB& rgb);

// CIE 1931 Color Matching Functions
// Reference: http://cvrl.ucl.ac.uk/cmfs.htm
static constexpr int COLOR_MATCHING_FUNC_SAMPLES = 85;

static constexpr double COLOR_MATCHING_FUNC_X[COLOR_MATCHING_FUNC_SAMPLES] = {
    0.001368, 0.002236, 0.004243, 0.007650, 0.014310, 0.023190, 0.043510, 0.077630,
    0.134380, 0.214770, 0.283900, 0.328500, 0.348280, 0.348060, 0.336200, 0.318700,
    0.290800, 0.251100, 0.195360, 0.142100, 0.095640, 0.057950, 0.032010, 0.014700,
    0.004900, 0.002400, 0.009300, 0.029100, 0.063270, 0.109600, 0.165500, 0.225750,
    0.290400, 0.359700, 0.433450, 0.512050, 0.594500, 0.678400, 0.762100, 0.842500,
    0.916300, 0.978600, 1.026300, 1.056700, 1.062200, 1.045600, 1.002600, 0.938400,
    0.854450, 0.751400, 0.642400, 0.541900, 0.447900, 0.360800, 0.283500, 0.218700,
    0.164900, 0.121200, 0.087400, 0.063600, 0.046770, 0.032900, 0.022700, 0.015840,
    0.011359, 0.008111, 0.005790, 0.004109, 0.002899, 0.002049, 0.001440, 0.001000,
    0.000690, 0.000476, 0.000332, 0.000235, 0.000166, 0.000117, 0.000083, 0.000059,
    0.000042
};

static constexpr double COLOR_MATCHING_FUNC_Y[COLOR_MATCHING_FUNC_SAMPLES] = {
    0.000039, 0.000064, 0.000120, 0.000217, 0.000396, 0.000640, 0.001210, 0.002180,
    0.004000, 0.007300, 0.011600, 0.016840, 0.023000, 0.029800, 0.038000, 0.048000,
    0.060000, 0.073900, 0.090980, 0.112600, 0.139020, 0.169300, 0.208020, 0.258600,
    0.323000, 0.407300, 0.503000, 0.608200, 0.710000, 0.793200, 0.862000, 0.914850,
    0.954000, 0.980300, 0.994950, 1.000000, 0.995000, 0.978600, 0.952000, 0.915400,
    0.870000, 0.816300, 0.757000, 0.694900, 0.631000, 0.566800, 0.503000, 0.441200,
    0.381000, 0.321000, 0.265000, 0.217000, 0.175000, 0.138200, 0.107000, 0.081600,
    0.061000, 0.044580, 0.032000, 0.023200, 0.017000, 0.011920, 0.008210, 0.005723,
    0.004102, 0.002929, 0.002091, 0.001484, 0.001047, 0.000740, 0.000520, 0.000361,
    0.000249, 0.000172, 0.000120, 0.000085, 0.000060, 0.000042, 0.000030, 0.000021,
    0.000015
};

static constexpr double COLOR_MATCHING_FUNC_Z[COLOR_MATCHING_FUNC_SAMPLES] = {
    0.006450, 0.010550, 0.020050, 0.036210, 0.067850, 0.110200, 0.207400, 0.371300,
    0.645600, 1.039050, 1.385600, 1.622960, 1.747060, 1.782600, 1.772110, 1.744100,
    1.669200, 1.528100, 1.287640, 1.041900, 0.812950, 0.616200, 0.465180, 0.353300,
    0.272000, 0.212300, 0.158200, 0.111700, 0.078250, 0.057250, 0.042160, 0.029840,
    0.020300, 0.013400, 0.008750, 0.005750, 0.003900, 0.002750, 0.002100, 0.001800,
    0.001650, 0.001400, 0.001100, 0.001000, 0.000800, 0.000600, 0.000340, 0.000240,
    0.000190, 0.000100, 0.000050, 0.000030, 0.000020, 0.000010, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000
};

}  // namespace CGL

#endif  // CGL_SPECTRUM_H