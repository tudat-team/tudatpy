//
// Created by ggarrett on 11-6-19.
//

// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wpragmas"

#include <boost/python/numpy.hpp>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

namespace bp = boost::python;
namespace bpn = boost::python::numpy;
namespace tsi = tudat::spice_interface;

void loadStandardSpiceKernels1() { tsi::loadStandardSpiceKernels(); }
void loadStandardSpiceKernels2(const std::vector<std::string> alternativeEphemerisKernels) {
    tsi::loadStandardSpiceKernels(alternativeEphemerisKernels);
}

BOOST_PYTHON_MODULE (_spice) {

    // TODO: See if BOOST_PYTHON_FUNCTION_OVERLOADS can replace line 29 and 30.
    // BOOST_PYTHON_FUNCTION_OVERLOADS(load_standard_spice_kernels_overloads, tsi::loadStandardSpiceKernels, 0,1)
    bp::def("load_standard_spice_kernels", loadStandardSpiceKernels1);
    bp::def("load_standard_spice_kernels", loadStandardSpiceKernels2);

}

#pragma GCC diagnostic pop