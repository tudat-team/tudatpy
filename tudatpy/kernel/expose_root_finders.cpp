//
// Created by ggarrett on 04-05-20.
//

#include "expose_root_finders.h"

namespace py = pybind11;
using namespace tudat::root_finders;

namespace tudatpy {
    void expose_root_finders(py::module &m) {

    	py::enum_<RootFinderType>(m, "RootFinderType")
    	        .value("bisection_root_finder", bisection_root_finder)
				.value("halley_root_finder", halley_root_finder)
				.value("newton_raphson_root_finder", newton_raphson_root_finder)
				.value("secant_root_finder", secant_root_finder)
				.export_values();

        py::class_<RootFinderSettings, std::shared_ptr<RootFinderSettings>>(m, "RootFinderSettings")
        	.def(py::init<
        	        const RootFinderType,
        	        const double,
        	        const unsigned int
        	        >(),
        	        py::arg("root_finder_type"),
        	        py::arg("termination_tolerance"),
        	        py::arg("maximum_number_of_iterations")
				);

    }
}
