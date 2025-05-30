//
// Created by 阮立心 on 25-4-14.
//

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "type.h"
#include <vector3D.h>

namespace CGL {
  /* ---------- 常量，来自 IEC 61966-2-1 ---------- */
  constexpr float kDecodeT  = 0.04045f;   // sRGB→线性 分段阈值 :contentReference[oaicite:0]{index=0}
  constexpr float kEncodeT  = 0.0031308f; // 线性→sRGB 分段阈值 :contentReference[oaicite:1]{index=1}
  constexpr float kInvGamma = 1.0f / 2.4f;
  constexpr float kA = 1.055f;
  constexpr float kB = 0.055f;
  constexpr float kS = 12.92f;

  inline double srgbToLinearCh(const double c)                 // sRGB→线性
  {
    return (c <= kDecodeT) ? c / kS
                           : std::pow((c + kB) / kA, 2.4f);
  }

  inline double linearToSrgbCh(const double c)                 // 线性→sRGB
  {
    return (c <= kEncodeT) ? c * kS
                           : kA * std::pow(c, kInvGamma) - kB;
  }

  inline Vector3D srgbToLinear(const Vector3D& srgb)             // in:[0,1] sRGB
  {
    return { srgbToLinearCh(srgb.r),
             srgbToLinearCh(srgb.g),
             srgbToLinearCh(srgb.b) };
  }

  inline Vector3D linearToSrgb(const Vector3D& lin)              // in:[0,1] 线性
  {
    return { linearToSrgbCh(lin.r),
             linearToSrgbCh(lin.g),
             linearToSrgbCh(lin.b) };
  }

using XYZ = Vector3D;
using RGB = Vector3D;

// Convert XYZ to sRGB color space
// XYZ to sRGB(D65)
// http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
inline RGB XYZ2RGB(const XYZ& xyz) {
  return RGB(
      3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
      -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
      0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z);
}

// Represents a SPD (Spectral Power Distribution) sampled at equal intervals
// Holds a series of wavelength and radiant flux samples
// Wavelengths are stored in [nm]
// SPDs with peaks narrower than the wavelength division width may not be represented accurately
// Originally wanted to store samples directly to represent non-equidistant SPDs, but data size was too large
class SPD {
 public:
  // Wavelength range stored in SPD
  static constexpr Real LAMBDA_MIN = 380;
  static constexpr Real LAMBDA_MAX = 780;

  // Number of wavelength divisions
  static constexpr size_t LAMBDA_SAMPLES = 80;

  // Divided wavelength width
  static constexpr Real LAMBDA_INTERVAL =
      (LAMBDA_MAX - LAMBDA_MIN) / LAMBDA_SAMPLES;

  std::array<Real, LAMBDA_SAMPLES> phi;  // Radiant flux

  // Initialize with 0
  SPD() : phi(){
    phi.fill(0);
  }

  // Initialize with a specific value
  explicit SPD(const Real& v) : phi(){
    phi.fill(v);
  }

  // Construct equidistant SPD from arbitrary wavelength and radiant flux sampling series
  // Assumes wavelengths and corresponding radiant flux are arranged in ascending order
  SPD(const std::vector<Real>& _lambda, const std::vector<Real>& _phi);

  // Returns the i-th radiant flux
  Real operator[](const size_t i) const {
    assert(i < SPD::LAMBDA_SAMPLES);
    return phi[i];
  }

  // Clear the SPD
  void clear() {
    phi.fill(0);
  }

  // Add spectral radiant flux
  void addPhi(const Real& _lambda, const Real& _phi);

  // Returns whether it's black
  bool isBlack() const {
    return std::all_of(phi.begin(), phi.end(), [](const Real& x) {return x != 0.0f;});
  }

  // Returns linearly interpolated radiant flux at specified wavelength
  // l : wavelength [nm]
  Real sample(const Real& l) const;

  // Convert to XYZ color space
  XYZ toXYZ() const;

  // sRGB color space conversion
  // XYZ to sRGB(D65)
  // http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  RGB toRGB() const {
    return clamp(XYZ2RGB(this->toXYZ()), Vector3D(0), Vector3D(INF_D));
  }

  // Operations
  SPD& operator+=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] += spd.phi[i];
    }
    return *this;
  }
  SPD& operator+=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] += k;
    }
    return *this;
  }
  SPD& operator-=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] -= spd.phi[i];
    }
    return *this;
  }
  SPD& operator-=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] -= k;
    }
    return *this;
  }
  SPD& operator*=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] *= spd.phi[i];
    }
    return *this;
  }
  SPD& operator*=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] *= k;
    }
    return *this;
  }
  SPD& operator/=(const SPD& spd) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] /= spd.phi[i];
    }
    return *this;
  }
  SPD& operator/=(const Real& k) {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] /= k;
    }
    return *this;
  }

  // Color matching function (CIE1931)
  // http://cvrl.ucl.ac.uk/cmfs.htm
  static constexpr int color_matching_func_samples = 85;
  static constexpr Real color_matching_func_x[color_matching_func_samples] = {
      0.001368000000f, 0.002236000000f, 0.004243000000f, 0.007650000000f,
      0.014310000000f, 0.023190000000f, 0.043510000000f, 0.077630000000f,
      0.134380000000f, 0.214770000000f, 0.283900000000f, 0.328500000000f,
      0.348280000000f, 0.348060000000f, 0.336200000000f, 0.318700000000f,
      0.290800000000f, 0.251100000000f, 0.195360000000f, 0.142100000000f,
      0.095640000000f, 0.057950010000f, 0.032010000000f, 0.014700000000f,
      0.004900000000f, 0.002400000000f, 0.009300000000f, 0.029100000000f,
      0.063270000000f, 0.109600000000f, 0.165500000000f, 0.225749900000f,
      0.290400000000f, 0.359700000000f, 0.433449900000f, 0.512050100000f,
      0.594500000000f, 0.678400000000f, 0.762100000000f, 0.842500000000f,
      0.916300000000f, 0.978600000000f, 1.026300000000f, 1.056700000000f,
      1.062200000000f, 1.045600000000f, 1.002600000000f, 0.938400000000f,
      0.854449900000f, 0.751400000000f, 0.642400000000f, 0.541900000000f,
      0.447900000000f, 0.360800000000f, 0.283500000000f, 0.218700000000f,
      0.164900000000f, 0.121200000000f, 0.087400000000f, 0.063600000000f,
      0.046770000000f, 0.032900000000f, 0.022700000000f, 0.015840000000f,
      0.011359160000f, 0.008110916000f, 0.005790346000f, 0.004109457000f,
      0.002899327000f, 0.002049190000f, 0.001439971000f, 0.000999949300f,
      0.000690078600f, 0.000476021300f, 0.000332301100f, 0.000234826100f,
      0.000166150500f, 0.000117413000f, 0.000083075270f, 0.000058706520f,
      0.000041509940f};
  static constexpr Real color_matching_func_y[color_matching_func_samples] = {
      0.000039000000f, 0.000064000000f, 0.000120000000f, 0.000217000000f,
      0.000396000000f, 0.000640000000f, 0.001210000000f, 0.002180000000f,
      0.004000000000f, 0.007300000000f, 0.011600000000f, 0.016840000000f,
      0.023000000000f, 0.029800000000f, 0.038000000000f, 0.048000000000f,
      0.060000000000f, 0.073900000000f, 0.090980000000f, 0.112600000000f,
      0.139020000000f, 0.169300000000f, 0.208020000000f, 0.258600000000f,
      0.323000000000f, 0.407300000000f, 0.503000000000f, 0.608200000000f,
      0.710000000000f, 0.793200000000f, 0.862000000000f, 0.914850100000f,
      0.954000000000f, 0.980300000000f, 0.994950100000f, 1.000000000000f,
      0.995000000000f, 0.978600000000f, 0.952000000000f, 0.915400000000f,
      0.870000000000f, 0.816300000000f, 0.757000000000f, 0.694900000000f,
      0.631000000000f, 0.566800000000f, 0.503000000000f, 0.441200000000f,
      0.381000000000f, 0.321000000000f, 0.265000000000f, 0.217000000000f,
      0.175000000000f, 0.138200000000f, 0.107000000000f, 0.081600000000f,
      0.061000000000f, 0.044580000000f, 0.032000000000f, 0.023200000000f,
      0.017000000000f, 0.011920000000f, 0.008210000000f, 0.005723000000f,
      0.004102000000f, 0.002929000000f, 0.002091000000f, 0.001484000000f,
      0.001047000000f, 0.000740000000f, 0.000520000000f, 0.000361100000f,
      0.000249200000f, 0.000171900000f, 0.000120000000f, 0.000084800000f,
      0.000060000000f, 0.000042400000f, 0.000030000000f, 0.000021200000f,
      0.000014990000f};
  static constexpr Real color_matching_func_z[color_matching_func_samples] = {
      0.006450001000f, 0.010549990000f, 0.020050010000f, 0.036210000000f,
      0.067850010000f, 0.110200000000f, 0.207400000000f, 0.371300000000f,
      0.645600000000f, 1.039050100000f, 1.385600000000f, 1.622960000000f,
      1.747060000000f, 1.782600000000f, 1.772110000000f, 1.744100000000f,
      1.669200000000f, 1.528100000000f, 1.287640000000f, 1.041900000000f,
      0.812950100000f, 0.616200000000f, 0.465180000000f, 0.353300000000f,
      0.272000000000f, 0.212300000000f, 0.158200000000f, 0.111700000000f,
      0.078249990000f, 0.057250010000f, 0.042160000000f, 0.029840000000f,
      0.020300000000f, 0.013400000000f, 0.008749999000f, 0.005749999000f,
      0.003900000000f, 0.002749999000f, 0.002100000000f, 0.001800000000f,
      0.001650001000f, 0.001400000000f, 0.001100000000f, 0.001000000000f,
      0.000800000000f, 0.000600000000f, 0.000340000000f, 0.000240000000f,
      0.000190000000f, 0.000100000000f, 0.000049999990f, 0.000030000000f,
      0.000020000000f, 0.000010000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f, 0.000000000000f, 0.000000000000f, 0.000000000000f,
      0.000000000000f};
};

// Operations between SPDs
// Performs operations on each element
inline SPD operator+(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] + spd2.phi[i];
  }
  return ret;
}
inline SPD operator-(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] - spd2.phi[i];
  }
  return ret;
}
inline SPD operator*(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] * spd2.phi[i];
  }
  return ret;
}
inline SPD operator/(const SPD& spd1, const SPD& spd2) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd1.phi[i] / spd2.phi[i];
  }
  return ret;
}

// Operations between SPD and Real
inline SPD operator+(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] + k;
  }
  return ret;
}
inline SPD operator+(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] + k;
  }
  return ret;
}
inline SPD operator-(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] - k;
  }
  return ret;
}
inline SPD operator-(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = k - spd.phi[i];
  }
  return ret;
}
inline SPD operator*(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] * k;
  }
  return ret;
}
inline SPD operator*(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] * k;
  }
  return ret;
}
inline SPD operator/(const SPD& spd, const Real& k) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = spd.phi[i] / k;
  }
  return ret;
}
inline SPD operator/(const Real& k, const SPD& spd) {
  SPD ret;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    ret.phi[i] = k / spd.phi[i];
  }
  return ret;
}

// Normalization
inline SPD normalize(const SPD& spd) {
  const auto m = std::max_element(spd.phi.begin(), spd.phi.end());
  return spd / *m;
}

// SPD output
inline std::ostream& operator<<(std::ostream& stream, const SPD& spd) {
  stream << std::setw(12) << "lambda" << std::setw(12) << "phi" << std::endl;
  for (size_t i = 0; i < SPD::LAMBDA_SAMPLES; ++i) {
    const Real lambda = SPD::LAMBDA_MIN + i * SPD::LAMBDA_INTERVAL;
    stream << std::setw(12) << lambda << std::setw(12) << spd.phi[i]
           << std::endl;
  }
  return stream;
}

// An RGB to Spectrum Conversion for Reflectances, Smits(2001)
SPD RGB2Spectrum(const RGB& rgb);

}  // namespace Prl2

#endif //SPECTRUM_H
