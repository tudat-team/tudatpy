//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDAT_EXPOSE_IO_H
#define TUDAT_EXPOSE_IO_H

#include <pybind11/pybind11.h>

#include <tudat/io/basicInputOutput.h>

namespace py = pybind11;

namespace tudatpy {

    void expose_io(py::module &m);

}


#endif //TUDAT_EXPOSE_IO_H
