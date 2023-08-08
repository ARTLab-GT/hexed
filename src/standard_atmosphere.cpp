#include <standard_atmosphere.hpp>
#include <constants.hpp>

namespace hexed
{

// altitude vs temperature
std::vector<std::array<double, 2>> alt_temp {
  {    0., 288.15},
  {11000., 216.65},
  {20000., 216.65},
  {32000., 228.65},
  {47000., 270.65},
  {51000., 270.65},
  {71000., 214.65},
  {80000., 196.65},
};

std::array<double, 2> standard_atmosphere(double alt_geom, double temp_offset)
{
  return standard_atmosphere_geopot(constants::earth_radius*alt_geom/(constants::earth_radius + alt_geom), temp_offset);
}

std::array<double, 2> standard_atmosphere_geopot(double alt_geopot, double temp_offset)
{
  double pres = constants::atmosphere;
  double temp;
  const double epsilon = 1e-10;
  HEXED_ASSERT(alt_geopot > alt_temp[0][0] - epsilon, "altitude is below model limit");
  HEXED_ASSERT(alt_geopot < alt_temp.back()[0] + epsilon, "altitude is above model limit");
  for (unsigned i_row = 0; alt_temp[i_row][0] - epsilon < alt_geopot && i_row < alt_temp.size() - 1; ++i_row) {
    double base_temp = alt_temp[i_row][1];
    double lapse = (alt_temp[i_row + 1][1] - base_temp)/(alt_temp[i_row + 1][0] - alt_temp[i_row][0]);
    double h_diff = std::min(alt_temp[i_row + 1][0], alt_geopot) - alt_temp[i_row][0];
    temp = base_temp + lapse*h_diff;
    if (std::abs(lapse) < epsilon) {
      pres *= std::exp(-constants::std_grav*h_diff/constants::specific_gas_air/base_temp);
    } else {
      pres *= std::pow(temp/base_temp, -constants::std_grav/lapse/constants::specific_gas_air);
    }
  }
  temp += temp_offset;
  double dens = pres/constants::specific_gas_air/temp;
  return {dens, pres};
}

}
