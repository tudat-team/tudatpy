// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wpragmas"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <iostream>
#include <Eigen/Core>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "conversion.h"
#include "eigen_numpy.h"

namespace bp = boost::python;
namespace bpn = boost::python::numpy;
namespace tss = tudat::simulation_setup;

BOOST_PYTHON_MODULE (_body) {

    Py_Initialize();
    SetupEigenConverters();

    bp::object Body(
            bp::class_<tss::Body>("Body", bp::init<Eigen::Vector6d const &>())
                    .def("get_ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame)
    );

}

#pragma GCC diagnostic pop