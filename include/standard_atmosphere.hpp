#ifndef HEXED_STANDARD_ATMOSPHERE_HPP_
#define HEXED_STANDARD_ATMOSPHERE_HPP_

#include <array>

namespace hexed
{

/*! \brief computes the [ICAO Standard atmosphere](http://www.aviationchief.com/uploads/9/2/0/9/92098238/icao_doc_7488_-_manual_of_icao_standard_atmosphere_-_3rd_edition_-_1994.pdf)
 * \details Valid from 0 to 80km.
 * \param alt_geom Geometric altitude (as opposed to geopotential altitude)
 * \param temp_offset Temperature will be incremented by `temp_offset` relative to the standard atmosphere without changing the pressure.
 * \returns density and pressure (in that order)
 * \see \ref units
 */
std::array<double, 2> standard_atmosphere(double alt_geom, double temp_offset = 0);

}
#endif
