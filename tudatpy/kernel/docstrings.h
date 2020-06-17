//
// Created by ggarrett on 24-04-20.
//

#ifndef TUDAT_DOSTRINGS_H
#define TUDAT_DOSTRINGS_H

#include <string>

namespace tudatpy {

    // constants
    std::string constants_docstring();

    // spice_interface
    std::string load_standard_spice_kernels_docstring();
    std::string clear_spice_kernels_docstring();

    // simulation_setup
    std::string body_settings_docstring();

}


#endif //TUDATBUNDLE_DOSTRINGS_H
