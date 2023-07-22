#ifndef HEXED_CONSTANTS_HPP_
#define HEXED_CONSTANTS_HPP_

/*! \file constants.hpp
 * Definition of mathematical and physical constants as `const` global variables
 * (in terms of [SI base units](https://en.wikipedia.org/wiki/SI_base_unit) when possible).
 * Note the doc descriptions with the mathematical symbols for the variables are mostly just
 * so the variables can be linked to in Doxygen but maybe they will be useful to someone...
 */

#include "math.hpp"

namespace hexed
{

//! \name mathematical constants
//!\{
const double pi = M_PI; //!< \brief \f$ \pi \f$
const double degree = pi/360;//!< \brief \f$ ^{\circ} \f$
//!\}

//! \name SI definitions (wikipedia)
//!\{
const double light_speed = 299792458; //!< \brief \f$ c \f$
const double plank = 6.62607015e-34; //!< \brief \f$ h \f$
const double boltzmann = 1.380649e-23; //!< \brief \f$ k \f$
const double avogadro = 6.02214076e23; //!< \brief \f$ mol \f$
//!\}

//! \name empirical values
//!\{
const double mol_mass_air = 28.9647e-3; //!< \brief \f$ \mu_{air} \f$ (from engineering toolbox)
const double earth_radius = 6356766; //!< \f$ \f$
//!\}

#define P(base, exp) math::pow((base), (exp))
//! \name derived quantities
//!\{
const double universal_gas = boltzmann*avogadro; //!< \brief \f$ R_u \f$
const double specific_gas_air = universal_gas/mol_mass_air; //!< \brief \f$ R_{air} \f$
const double stefan_boltzmann = 2*P(pi, 5)*P(boltzmann, 4)/(15*P(light_speed, 2)*P(plank, 3)); //!< \brief \f$ \sigma \f$
//!\}

//! \name unit definitions
//! Conversion factors from various units to standard SI units, based on the official exact definitions.
//! \see \ref units
//!\{
const double meter = 1; //!< \brief m
const double kilogram = 1; //!< \brief kg
const double second = 1; //!< \brief s
const double kelvin = 1; //!< \brief K
const double std_grav = 9.80665; //!< \brief \f$ g_0 \f$ [standard gravity](https://en.wikipedia.org/wiki/Standard_gravity)
const double foot = 0.3048; //!< \brief [ft](https://en.wikipedia.org/wiki/International_yard_and_pound)
const double yard = 3*foot; //!< \brief [yd](https://en.wikipedia.org/wiki/International_yard_and_pound)
const double inch = foot/12; //!< \brief in
const double mile = 5280*foot; //!< \brief mi
const double nautical_mile = 1852; //!< nmi
const double pound_mass = 0.45359237; //!< \brief [lb_m](https://en.wikipedia.org/wiki/International_yard_and_pound)
//! \brief [slug](https://en.wikipedia.org/wiki/Slug_(unit))
//! \todo I'm not sure if this is an exact definition or an approximate conversion. Find an official answer somewhere...
const double slug = 32.1740*pound_mass;
const double minute = 60*second; //!< \brief min
const double hour = 60*minute; //!< \brief hr
const double rankine = 9./5.; //!< \brief \f$ ^{\circ} \f$ R
const double pound_force = slug*foot; //!< lb_f
const double knot = nautical_mile/hour; //!< \brief kn
const double zero_celsius = 273.15; //!< \brief 0 \f$ ^{\circ} \f$ [C](https://en.wikipedia.org/wiki/Celsius) (not technically a unit but still important)
const double atmosphere = 101325; //!< \brief [atm](https://en.wikipedia.org/wiki/Standard_atmosphere_(unit))
//!\}
#undef P

}
#endif
