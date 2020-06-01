//
// Created by ggarrett on 04-05-20.
//

#include "expose_root_finders.h"

namespace py = pybind11;
using namespace tudat::root_finders;

namespace tudatpy {
    void expose_root_finders(py::module &m) {

        py::class_<RootFinderSettings,
                std::shared_ptr<RootFinderSettings>> RootFinderSettings_(m, "RootFinderSettings");

    }
}
