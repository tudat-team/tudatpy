// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include <boost/python.hpp>
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "conversion.h"

namespace bp = boost::python;
namespace tbc = tudat::celestial_body_constants;

BOOST_PYTHON_MODULE (celestial) {

    bp::scope().attr("__doc__") = " EARTH_EQUATORIAL_RADIUS"
                                  " EARTH_FLATTENING_FACTOR"
                                  " EARTH_GEODESY_NORMALIZED_J2"
                                  " SUN_GRAVITATIONAL_PARAMETER"
                                  " MERCURY_GRAVITATIONAL_PARAMETER"
                                  " VENUS_GRAVITATIONAL_PARAMETER"
                                  " EARTH_GRAVITATIONAL_PARAMETER"
                                  " MOON_GRAVITATIONAL_PARAMETER"
                                  " MARS_GRAVITATIONAL_PARAMETER"
                                  " JUPITER_GRAVITATIONAL_PARAMETER"
                                  " SATURN_GRAVITATIONAL_PARAMETER"
                                  " URANUS_GRAVITATIONAL_PARAMETER"
                                  " NEPTUNE_GRAVITATIONAL_PARAMETER"
                                  " PLUTO_GRAVITATIONAL_PARAMETER";

    // Definition of celestial body constants.
    bp::scope().attr("EARTH_EQUATORIAL_RADIUS") = tbc::EARTH_EQUATORIAL_RADIUS;
    bp::scope().attr("EARTH_FLATTENING_FACTOR") = tbc::EARTH_FLATTENING_FACTOR;
    bp::scope().attr("EARTH_GEODESY_NORMALIZED_J2") = tbc::EARTH_GEODESY_NORMALIZED_J2;
    bp::scope().attr("SUN_GRAVITATIONAL_PARAMETER") = tbc::SUN_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("MERCURY_GRAVITATIONAL_PARAMETER") = tbc::MERCURY_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("VENUS_GRAVITATIONAL_PARAMETER") = tbc::VENUS_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("EARTH_GRAVITATIONAL_PARAMETER") = tbc::EARTH_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("MOON_GRAVITATIONAL_PARAMETER") = tbc::MOON_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("MARS_GRAVITATIONAL_PARAMETER") = tbc::MARS_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("JUPITER_GRAVITATIONAL_PARAMETER") = tbc::JUPITER_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("SATURN_GRAVITATIONAL_PARAMETER") = tbc::SATURN_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("URANUS_GRAVITATIONAL_PARAMETER") = tbc::URANUS_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("NEPTUNE_GRAVITATIONAL_PARAMETER") = tbc::NEPTUNE_GRAVITATIONAL_PARAMETER;
    bp::scope().attr("PLUTO_GRAVITATIONAL_PARAMETER") = tbc::PLUTO_GRAVITATIONAL_PARAMETER;

    // Definition of contained maps.
    bp::scope().attr("planet_names") = toPythonDict(tbc::planetNames);
    bp::scope().attr("planet_id_numbers") = toPythonDict(tbc::planetIdNumbers);

}

#pragma GCC diagnostic pop