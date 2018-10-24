/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <getopt.h>

#include "Tudat/JsonInterface/jsonInterfaceVariational.h"

namespace tudat
{

namespace json_interface
{

template class JsonVariationalEquationsSimulationManager< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//template class JsonVariationalEquationsSimulationManager< Time, double >;
template class JsonVariationalEquationsSimulationManager< double, long double >;
//template class JsonVariationalEquationsSimulationManager< Time, long double >;
#endif

}

}
