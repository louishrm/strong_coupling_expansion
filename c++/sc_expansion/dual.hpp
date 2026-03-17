#pragma once

#include <cmath>

struct Dual {

  double value;
  double derivative;

  Dual(double v = 0.0, double d = 0.0) : value(v), derivative(d) {}

  Dual operator+(const Dual &other) const { return Dual(value + other.value, derivative + other.derivative); }
  Dual operator-(const Dual &other) const { return Dual(value - other.value, derivative - other.derivative); }
  Dual operator*(const Dual &other) const { return Dual(value * other.value, value * other.derivative + derivative * other.value); }
  Dual operator/(const Dual &other) const {
    double v2 = other.value * other.value;
    return Dual(value / other.value, (derivative * other.value - value * other.derivative) / v2);
  }

  Dual operator-() const { return Dual(-value, -derivative); }
};

inline Dual exp(const Dual &xi) {
  double exp_val = std::exp(xi.value);
  return Dual(exp_val, exp_val * xi.derivative);
}

inline Dual pow(const Dual &xi, double p) {

  double val = std::pow(xi.value, p);
  double der = p * std::pow(xi.value, p - 1) * xi.derivative;
  return Dual(val, der);
}

inline Dual sqrt(const Dual &xi) {
  double sqrt_val = std::sqrt(xi.value);
  double der      = 0.5 * xi.derivative / sqrt_val;
  return Dual(sqrt_val, der);
}