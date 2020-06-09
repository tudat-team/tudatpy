//
// Created by ggarrett on 09-06-20.
//

#include "expose_io.h"

namespace py = pybind11;

using namespace tudat::input_output;

namespace tudatpy {

    void expose_io(py::module &m){

        m.def("get_tudat_data_path",
              &getTudatDataPath,
#ifdef TUDAT_INSTALL
              py::arg("state") = "install"
#else
              py::arg("state") = "build"
#endif
        );

        m.def("get_tudat_root_path", &getTudatRootPath);

        m.def("get_tudat_src_path", &get_tudat_src_dir_path);

    };

}
