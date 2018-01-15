/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cmath>
#include <ctime>

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
{
namespace basic_mathematics
{

//! Get global random number generator.
GlobalRandomNumberGeneratorType& getGlobalRandomNumberGenerator( )
{
  static GlobalRandomNumberGeneratorType globalRandomNumberGenerator(
              static_cast< unsigned int >( std::time( 0 ) ) );
  return globalRandomNumberGenerator;
}

} // namespace basic_mathematics
} // namespace tudat
