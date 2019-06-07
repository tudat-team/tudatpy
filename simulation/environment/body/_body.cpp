// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wpragmas"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
//#include <Eigen/Core>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "util/conversion.h"

namespace bp = boost::python;
namespace bpn = boost::python::numpy;
namespace tss = tudat::simulation_setup;

//bp::init<Eigen::Vector6d>()
Eigen::Vector6d numpyToEigen6d(bpn::ndarray array) {
    Eigen::Vector6d result;
    result << array[1], array[2], array[3], array[4], array[5];
    return result;
};


struct numpyToEigen6d {
    static void *convertible(bpn::ndarray *obj_ptr) {
        if (!obj_ptr->size() == 6) return 0;
        return obj_ptr;
    }

    static void construct(
            bpn::ndarray *obj_ptr,
            boost::python::converter::rvalue_from_python_stage1_data *data
    ) {
        const Eigen::Vector6d result;
        result <<
               obj_ptr->operator[](0),
                obj_ptr->operator[](1),
                obj_ptr->operator[](2),
                obj_ptr->operator[](3),
                obj_ptr->operator[](4),
                obj_ptr->operator[](5);

    }

};


BOOST_PYTHON_MODULE (_body) {

    Py_Initialize();
    bpn::initialize();

    bp::object Body(
            bp::class_<tss::Body>("Body", bp::init<Eigen::Vector6d>())
                    .def("get_ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame)
    );

//    object body_obj = Body(numpyToEigen6d())

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