/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H
#define TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H

#include "json/src/json.hpp"

// #include <Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h>

namespace tudat
{

namespace json_interface
{

using json = nlohmann::json;

json emptyJson = json::parse( "{ }" );

} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H
