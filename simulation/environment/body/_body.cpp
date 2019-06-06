// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wpragmas"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
//#include <boost/numpy.hpp>
//#include <Eigen/Core>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
//#include "util/conversion.h"

namespace bp = boost::python;
namespace bpn = boost::python::numpy;
namespace tss = tudat::simulation_setup;

BOOST_PYTHON_MODULE (_body) {

        bp::class_<tss::Body>("Body", bp::init<Eigen::Vector6d>())
                .def("get_ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame)
        ;


//        class_<Ctor>("Ctor", init<std::string>())
//        .def(init<double, double>())
//        .def("greet", &Ctor::greet)
//        .def("set", &Ctor::set)
//        ;

//        // Definition of celestial body constants.
//        bp::scope().attr("EARTH_EQUATORIAL_RADIUS") = tbc::EARTH_EQUATORIAL_RADIUS;
//        bp::scope().attr("EARTH_FLATTENING_FACTOR") = tbc::EARTH_FLATTENING_FACTOR;
//        bp::scope().attr("EARTH_GEODESY_NORMALIZED_J2") = tbc::EARTH_GEODESY_NORMALIZED_J2;
//        bp::scope().attr("SUN_GRAVITATIONAL_PARAMETER") = tbc::SUN_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("MERCURY_GRAVITATIONAL_PARAMETER") = tbc::MERCURY_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("VENUS_GRAVITATIONAL_PARAMETER") = tbc::VENUS_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("EARTH_GRAVITATIONAL_PARAMETER") = tbc::EARTH_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("MOON_GRAVITATIONAL_PARAMETER") = tbc::MOON_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("MARS_GRAVITATIONAL_PARAMETER") = tbc::MARS_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("JUPITER_GRAVITATIONAL_PARAMETER") = tbc::JUPITER_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("SATURN_GRAVITATIONAL_PARAMETER") = tbc::SATURN_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("URANUS_GRAVITATIONAL_PARAMETER") = tbc::URANUS_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("NEPTUNE_GRAVITATIONAL_PARAMETER") = tbc::NEPTUNE_GRAVITATIONAL_PARAMETER;
//        bp::scope().attr("PLUTO_GRAVITATIONAL_PARAMETER") = tbc::PLUTO_GRAVITATIONAL_PARAMETER;
//
//        // Definition of contained maps.
//        bp::scope().attr("planet_names") = toPythonDict(tbc::planetNames);
//        bp::scope().attr("planet_id_numbers") = toPythonDict(tbc::planetIdNumbers);

}

#pragma GCC diagnostic pop