#ifndef CGL_IOR_H
#define CGL_IOR_H

#include <cmath>
#include "CGL/type.h"

namespace CGL {

class IOREquation {
 public:
    virtual ~IOREquation() = default;

    virtual double ior(double lambda) const = 0;
};

class ConstantIOR final : public IOREquation {
 public:
  double ior_value;

  explicit ConstantIOR(const double ior_value) : ior_value(ior_value) {}

  double ior(double lambda) const override { return ior_value; }
};

class CauchyEquation final : public IOREquation {
 public:
  double A;
  double B;

  CauchyEquation() : A(0.0), B(0.0) {}
  CauchyEquation(const double A, const double B) : A(A), B(B) {}

  double ior(const double lambda) const override {
    const double mu = 1e-3 * lambda;
    return A + B / (mu * mu);
  }
};

inline CauchyEquation fit_cauchy(const double nD, const double nF) {
  constexpr double lambdaD = 0.5893;
  constexpr double lambdaF = 0.4861;

  constexpr double alpha1 = 1 / (lambdaD * lambdaD);
  constexpr double alpha2 = 1 / (lambdaF * lambdaF);
  constexpr double det = 1 / (alpha2 - alpha1);

  return CauchyEquation(det * (alpha2 * nD - alpha1 * nF), det * (nF - nD));
}

class SellmeierCoefficient final: public IOREquation {
 public:
  double alpha;
  double B[3];
  double C[3];

  SellmeierCoefficient() : alpha(0.0), B{0.0, 0.0, 0.0}, C{0.0, 0.0, 0.0} {}
  SellmeierCoefficient(const double ior550, const double b0, const double b1, const double b2, const double c0, const double c1, const double c2) {
    B[0] = b0;
    B[1] = b1;
    B[2] = b2;

    C[0] = c0;
    C[1] = c1;
    C[2] = c2;

    constexpr double lambda2 = 0.550 * 0.550;
    alpha = ior550 - std::sqrt(1 + B[0] * lambda2 / (lambda2 - C[0]) +
                               B[1] * lambda2 / (lambda2 - C[1]) +
                               B[2] * lambda2 / (lambda2 - C[2]));
  }

  double ior(const double lambda) const override {
    const double l = lambda * 1e-3;
    const double l2 = l * l;
    return std::sqrt(1 + B[0] * l2 / (l2 - C[0]) + B[1] * l2 / (l2 - C[1]) +
                     B[2] * l2 / (l2 - C[2])) +
           alpha;
  }
};

}  // namespace CGL

#endif  // CGL_IOR_H