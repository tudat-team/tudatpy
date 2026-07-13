/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_environment.h"

#include "coma/expose_coma.h"
#include "expose_ilrs.h"
#include "expose_space_weather.h"
#include "missions/expose_missions.h"

namespace py = pybind11;

namespace tudatpy
{
namespace data_access
{
namespace environment
{

void expose_environment( py::module& m )
{
    auto coma_submodule = m.def_submodule( "coma" );
    coma::expose_coma( coma_submodule );

    auto missions_submodule = m.def_submodule( "missions" );
    missions::expose_missions( missions_submodule );

    auto ilrs_submodule = m.def_submodule( "ilrs" );
    ilrs::expose_ilrs( ilrs_submodule );

    auto space_weather_submodule = m.def_submodule( "space_weather" );
    space_weather::expose_space_weather( space_weather_submodule );
}

}  // namespace environment
}  // namespace data_access
}  // namespace tudatpy
