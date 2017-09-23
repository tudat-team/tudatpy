#include <iostream>

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
