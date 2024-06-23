#pragma once
#include <cmath>

inline double alpha1(double Beta, double tStep) { return 1. / Beta / std::pow(tStep, 2); }
inline double alpha2(double Beta, double tStep) { return 1. / Beta / tStep; }
inline double alpha3(double Beta) { return (1 - 2 * Beta) / 2 / Beta; }
inline double alpha4(double Beta, double Gamma, double tStep) { return Gamma / Beta / tStep; }
inline double alpha5(double Beta, double Gamma) { return 1 - Gamma / Beta; }
inline double alpha6(double Beta, double Gamma, double tStep) { return (1 - Gamma / Beta / 2) * tStep; }