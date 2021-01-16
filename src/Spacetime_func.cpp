#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <Spacetime_func.hpp>

namespace cartdg
{

std::vector<double> Constant_func::operator()(std::vector<double> pos, double time)
{
  return value;
}

Isentropic_vortex::Isentropic_vortex(std::vector<double> state) : freestream(state) {}

std::vector<double> Isentropic_vortex::operator()(std::vector<double> pos, double time)
{
  double beta = 5;
  double r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
  double u = beta/(2*M_PI)*std::exp((1 - r*r)/2)*-pos[1];
  double v = beta/(2*M_PI)*std::exp((1 - r*r)/2)* pos[0];
  double T = 1. - (heat_rat - 1.)*beta*beta/(8*M_PI*M_PI)*std::exp(1. - r*r);
  double rho = std::pow(T, 1./(heat_rat - 1.));
  return std::vector<double> {rho*u, rho*v, rho, T/(heat_rat - 1.) + 0.5*rho*(u*u + v*v)};
}

}
