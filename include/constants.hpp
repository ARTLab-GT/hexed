#ifndef HEXED_CONSTANTS_HPP_
#define HEXED_CONSTANTS_HPP_

#include "math.hpp"

/*! \brief Mathematical and physical constants.
 * \details Definition of mathematical and physical constants, including unit conversions, as `const` global variables.
 * If you're in Python, you can import with `from hexed.cpp.constants import *`.
 * \note
 * - Whenever I say here something is "exact", I really mean it is exact up to the
 *   [inherent limitations](https://docs.python.org/3/tutorial/floatingpoint.html) of floating-point representations.
 * - the doc descriptions with the mathematical symbols for the variables are mostly just
 *   so the variables can be linked to in Doxygen but maybe they will be useful to someone...
 */
namespace hexed::constants
{

//! \name mathematical constants
//!\{
const double pi = M_PI; //!< \brief \f$ \pi \f$ \details obviously a computed value, not a definition, but should be exact to machine epsilon
const double degree = pi/360;//!< \brief \f$ ^{\circ} \f$
//!\}

/*! \name SI definitions
 * \brief Exact numerical values of physical constants based on the
 * [2019 redefinition](https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units#Redefinition)
 * of [SI base units](https://en.wikipedia.org/wiki/SI_base_unit).
 * \details I expect only a few of these constants will be useful for CFD, but might as well include them all for completeness.
 */
//!\{
const double caesium133_freq = 9192631770; //!< \brief \f$ \Delta \nu_{Cs} \f$
const double light_speed = 299792458; //!< \brief \f$ c \f$
const double plank = 6.62607015e-34; //!< \brief \f$ h \f$
const double electron_charge = 1.602176634e-19; //!< \brief \f$ e \f$
const double boltzmann = 1.380649e-23; //!< \brief \f$ k \f$
const double avogadro = 6.02214076e23; //!< \brief \f$ mol \f$
const double lum_eff_540e12 = 683; //!< \brief \f$ K_cd \f$
//!\}

//! \name empirical values
//!\{
const double mol_mass_air = 28.9647e-3; //!< \brief \f$ \mu_{air} \f$ (from engineering toolbox)
const double earth_radius = 6356766; //!< \brief \f$ \f$
//!\}

//! \name derived quantities
//!\{
const double universal_gas = boltzmann*avogadro; //!< \brief \f$ R_u \f$
const double specific_gas_air = universal_gas/mol_mass_air; //!< \brief \f$ R_{air} \f$
const double stefan_boltzmann = 2*math::pow(pi, 5)*math::pow(boltzmann, 4)/(15*math::pow(light_speed, 2)*math::pow(plank, 3)); //!< \brief \f$ \sigma \f$
//!\}

/*! \name unit definitions
 * Conversion factors from various units to standard SI units. Since hexed works \ref units "exclusively in SI units",
 * these are the values of each of these units _in m, kg, s, K_.
 * E.g. If you have a wing with a chord of 3 feet, then you should tell hexed it's chord is `3*foot`.
 * If hexed told you your drag force is 10 (implying newtons) and your (ill advised) chief engineer wants it in pounds,
 * the number they're looking for is `10/pound_force` (note that `pound_force != pound_mass`).
 * All quantities here are numerically exact according to some form of official standard.
 */
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
const double pound_mass = 0.45359237; //!< \brief [\f$ lb_m \f$](https://en.wikipedia.org/wiki/International_yard_and_pound)
/*! \brief [slug](https://en.wikipedia.org/wiki/Slug_(unit))
 * \details This definition may be backward from the perspective of slug, ft, s proponents, but
 * \f[ slug \cdot ft/s^2 = lb_f = g_0 lb_m \f]
 * is an exact equation for the slug since the pound mass has been exactly defined by the
 * [international yard and pound](https://en.wikipedia.org/wiki/International_yard_and_pound).
 */
const double slug = pound_mass*std_grav/foot*second*second;
const double minute = 60*second; //!< \brief min
const double hour = 60*minute; //!< \brief hr
const double rankine = 9./5.; //!< \brief \f$ ^{\circ} \f$ R
const double pound_force = std_grav*pound_mass; //!< \brief \f$ lb_f \f$
const double knot = nautical_mile/hour; //!< \brief kn
const double zero_celsius = 273.15; //!< \brief 0 \f$ ^{\circ} \f$ [C](https://en.wikipedia.org/wiki/Celsius) (not technically a unit but still important)
const double atmosphere = 101325; //!< \brief [atm](https://en.wikipedia.org/wiki/Standard_atmosphere_(unit))
//!\}

}
#endif
