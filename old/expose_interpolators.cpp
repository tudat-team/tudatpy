//
// Created by ggarrett on 29-04-20.
//

#include "expose_interpolators.h"

namespace py = pybind11;

namespace ti = tudat::interpolators;

namespace tudatpy {

    void expose_interpolators(py::module &m) {

        py::class_<ti::LagrangeInterpolatorSettings,
                std::shared_ptr<ti::LagrangeInterpolatorSettings>> lagrange_interpolator_settings(
                m,
                "LagrangeInterpolatorSettings",
                "PLACEHOLDER_CLASS");

    };

}
