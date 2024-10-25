//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATPY_EXPOSE_IO_H
#define TUDATPY_EXPOSE_IO_H

#include <pybind11/pybind11.h>

#include <tudat/io/basicInputOutput.h>

namespace py = pybind11;

namespace tudatpy {

namespace data {

    void expose_data(py::module &m);

}
}


#endif //TUDATPY_EXPOSE_IO_H
