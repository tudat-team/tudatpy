//
// Created by ggarrett on 24-04-20.
//

#ifndef TUDATPY_DOCSTRINGS_H
#define TUDATPY_DOCSTRINGS_H

#include <string>

namespace tudatpy {

    // constants
    std::string constants_docstring();

    // spice_interface
    std::string load_standard_kernels_docstring();
    std::string clear_spice_kernels_docstring();

    // simulation_setup
    std::string body_settings_docstring();

}


#endif //TUDATPY_DOCSTRINGS_H
