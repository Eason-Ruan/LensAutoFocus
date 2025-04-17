#ifndef CGL_IOR_H
#define CGL_IOR_H

#include <cmath>
#include "CGL/types.h"

namespace CGL {

class IOREquation {
 public:
  virtual double ior(double lambda) const = 0;
};

class ConstantIOR : public IOREquation {
 public:
  double ior_value;

  explicit ConstantIOR(double ior_value) : ior_value(ior_value) {}

  double ior(double lambda) const override { return ior_value; }
};

class CauchyEquation : public IOREquation {
 public:
  double A;
  double B;

  CauchyEquation() : A(0.0), B(0.0) {}
  CauchyEquation(double A, double B) : A(A), B(B) {}

  double ior(double lambda) const override {
    const double mu = 1e-3 * lambda;
    return A + B / (mu * mu);
  }
};

inline CauchyEquation fit_cauchy(double nD, double nF) {
  constexpr double lambdaD = 0.5893;
  constexpr double lambdaF = 0.4861;

  const double alpha1 = 1 / (lambdaD * lambdaD);
  const double alpha2 = 1 / (lambdaF * lambdaF);
  const double det = 1 / (alpha2 - alpha1);

  return CauchyEquation(det * (alpha2 * nD - alpha1 * nF), det * (nF - nD));
}

class SellmeierCoefficient : public IOREquation {
 public:
  double alpha;
  double B[3];
  double C[3];

  SellmeierCoefficient() : alpha(0.0), B{0.0, 0.0, 0.0}, C{0.0, 0.0, 0.0} {}
  SellmeierCoefficient(double ior550, double b0, double b1, double b2, double c0, double c1, double c2) {
    B[0] = b0;
    B[1] = b1;
    B[2] = b2;

    C[0] = c0;
    C[1] = c1;
    C[2] = c2;

    const double lambda2 = 0.550 * 0.550;
    alpha = ior550 - std::sqrt(1 + B[0] * lambda2 / (lambda2 - C[0]) +
                               B[1] * lambda2 / (lambda2 - C[1]) +
                               B[2] * lambda2 / (lambda2 - C[2]));
  }

  double ior(double lambda) const override {
    const double l = lambda * 1e-3;
    const double l2 = l * l;
    return std::sqrt(1 + B[0] * l2 / (l2 - C[0]) + B[1] * l2 / (l2 - C[1]) +
                     B[2] * l2 / (l2 - C[2])) +
           alpha;
  }
};

}  // namespace CGL

#endif  // CGL_IOR_H