//
// Created by ggarrett on 02-05-20.
//

#include "expose_ephemerides.h"

namespace py = pybind11;
using namespace tudat::ephemerides;

namespace tudatpy {

    void expose_ephemerides(py::module &m) {

        py::class_<RotationalEphemeris,
                std::shared_ptr<RotationalEphemeris>> RotationalEphemeris_(m, "RotationalEphemeris");

        m.def("transform_state_to_global_frame",
              &tudat::ephemerides::transformStateToGlobalFrame<double, double>,
              py::arg("state_in_local_frame"),
              py::arg("current_time"),
              py::arg("rotational_ephemeris")
        );

    };

}