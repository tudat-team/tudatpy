/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef ORBITDETERMINATIONTESTCASEUTILITIES_H
#define ORBITDETERMINATIONTESTCASEUTILITIES_H

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{


Eigen::VectorXd getDefaultInitialParameterPerturbation( );

}  // namespace unit_tests

}  // namespace tudat

#endif  // ORBITDETERMINATIONTESTCASEUTILITIES_H
