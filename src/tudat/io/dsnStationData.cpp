/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/dsnStationData.h"

namespace tudat
{

namespace input_output
{

std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( )
{
    return { { 10, { "DSS-12", "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26", "DSS-27" } },
             { 40, { "DSS-34", "DSS-35", "DSS-36", "DSS-42", "DSS-43", "DSS-45" } },
             { 60, { "DSS-53", "DSS-54", "DSS-55", "DSS-56", "DSS-61", "DSS-63", "DSS-65" } } };
}

}  // namespace input_output

}  // namespace tudat
