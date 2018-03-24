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

#include "Tudat/External/SofaInterface/earthOrientation.h"

namespace tudat
{

namespace unit_tests
{

Eigen::Matrix3d convertArrayToMatrix(
        const double array[3][3] );

Eigen::Matrix3d getSofaEarthOrientationExamples(
        const int calculationCase,
        const double dXInMas = 0.1750,
        const double dYInMas = -0.2259,
        const double xPInAs = 0.0349282,
        const double yPInAs = 0.48331639,
        const double ut1Correction = -0.072073685 );

}

}
