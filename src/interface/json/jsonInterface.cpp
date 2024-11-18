/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <getopt.h>

#include "tudat/interface/json/jsonInterface.h"

namespace tudat
{

namespace json_interface
{

//template class JsonSimulationManager< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//template class JsonSimulationManager< double, long double >;
////template class JsonSimulationManager< Time, double >;
////template class JsonSimulationManager< Time, long double >;
#endif

}

}
