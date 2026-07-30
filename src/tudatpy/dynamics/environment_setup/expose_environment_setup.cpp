/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_environment_setup.h"

#include "aerodynamic_coefficients/expose_aerodynamic_coefficients.h"
#include "gravity_field_variation/expose_gravity_field_variation.h"
#include "radiation_pressure/expose_radiation_pressure.h"
#include "space_time/expose_space_time.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{

void expose_environment_setup_types( py::module& m )
{
    auto aerodynamic_coefficient_setup = m.def_submodule( "aerodynamic_coefficients" );
    aerodynamic_coefficients::expose_aerodynamic_coefficient_types( aerodynamic_coefficient_setup );

    auto radiation_pressure_setup = m.def_submodule( "radiation_pressure" );
    radiation_pressure::expose_radiation_pressure_types( radiation_pressure_setup );

    auto gravity_variation_setup = m.def_submodule( "gravity_field_variation" );
    gravity_field_variation::expose_gravity_field_variation_types( gravity_variation_setup );

    auto space_time_setup = m.def_submodule( "space_time" );
    space_time::expose_space_time_types( space_time_setup );
}

}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
