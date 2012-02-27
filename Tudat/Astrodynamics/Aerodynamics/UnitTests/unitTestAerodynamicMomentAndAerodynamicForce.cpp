/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110622    F.M. Engelen      First creation of code.
 *      110822    D. Dirkx          Removed no longer necessary unit tests.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *
 *    References
 *
 */

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <limits>
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicForce.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicMoment.h"

//! Test implemntation of aerodynamic moment and aerodynamic force.
int main( )
{
    // Using declarations.
    using Eigen::Quaternion;
    using Eigen::AngleAxisd;
    using std::cerr;
    using std::endl;
    using std::fabs;
    using namespace tudat;

    bool isAerodynamicMomentBroken = false;
    bool isAerodynamicForceBroken = false;

    // Initialise objects.
    Eigen::Vector3d forceCoefficients;
    forceCoefficients( 0 ) = 1.1;
    forceCoefficients( 1 ) = 1.2;
    forceCoefficients( 2 ) = 1.3;

    Eigen::Vector3d momentCoefficients;
    momentCoefficients( 0 ) = 0.0;
    momentCoefficients( 1 ) = 1.0;
    momentCoefficients( 2 ) = 0.0;

    double dynamicPressure = 50.0;
    double referenceArea = 2.2;
    double referenceLength = 3.2;


    AerodynamicCoefficientInterface aerodynamicCoefficientInterface;
    aerodynamicCoefficientInterface.setCurrentForceCoefficients( forceCoefficients );
    aerodynamicCoefficientInterface.setCurrentMomentCoefficients( momentCoefficients );
    aerodynamicCoefficientInterface.setReferenceArea( referenceArea );
    aerodynamicCoefficientInterface.setReferenceLength( referenceLength );

    Eigen::Vector3d momentArm;
    momentArm( 0 ) = -2.0;
    momentArm( 1 ) = 0.0;
    momentArm( 2 ) = 0.0;

    AerodynamicForce aerodynamicForce( &aerodynamicCoefficientInterface );
    aerodynamicForce.setDynamicPressure( dynamicPressure );

    Eigen::Vector3d force;

    // Test 1. Test the force model without setting a transformation.

    // Calculated force.
    State dummyState;
    aerodynamicForce.computeForce( &dummyState );
    force = aerodynamicForce.getForce( );

    // Expected force.
    Eigen::Vector3d expectedForce = forceCoefficients * dynamicPressure * referenceArea;

    // Error in calculation.
    Eigen::Vector3d errorInForce;
    errorInForce( 0 ) = fabs( expectedForce( 0 ) - force( 0 ) );
    errorInForce( 1 ) = fabs( expectedForce( 1 ) - force( 1 ) );
    errorInForce( 2 ) = fabs( expectedForce( 2 ) - force( 2 ) );

    if ( std::numeric_limits< double >::epsilon( ) < errorInForce.sum( ) / 3.0 /
         ( dynamicPressure * referenceArea ) )
    {
        isAerodynamicForceBroken = true;
        cerr << "Test 1 of unitTestAerodynamicMomentAndAerodynamicForce failed" << endl;
    }


     // Test 2: Check the moment model without extra moment due to the force.
    AerodynamicMoment aerodynamicMoment( &aerodynamicCoefficientInterface );
    aerodynamicMoment.setDynamicPressure( dynamicPressure );

    Eigen::Vector3d moment;
    Eigen::Vector3d expectedMoment;
    Eigen::Vector3d errorInMoment;

    // Calculate moment.
    aerodynamicMoment.computeMoment( &dummyState );
    moment = aerodynamicMoment.getMoment( );

    // Calculate expected moment.
    expectedMoment = dynamicPressure * referenceArea * referenceLength * momentCoefficients;

    // Error in calculation.
    errorInMoment( 0 ) = fabs( expectedMoment( 0 ) - moment( 0 ) );
    errorInMoment( 1 ) = fabs( expectedMoment( 1 ) - moment( 1 ) );
    errorInMoment( 2 ) = fabs( expectedMoment( 2 ) - moment( 2 ) );

    if ( std::numeric_limits< double >::epsilon( ) < errorInMoment.sum( ) / 3.0 /
         ( dynamicPressure * referenceArea * referenceLength ) )
    {
        isAerodynamicMomentBroken = true;
        cerr << "Test 2 of unitTestAerodynamicMomentAndAerodynamicForce failed" << endl;
    }

    if ( isAerodynamicMomentBroken || isAerodynamicForceBroken )
    {
        cerr << "testAerodynamicMomentAndAerodynamicForce failed!" << std::endl;
    }

    return ( isAerodynamicMomentBroken || isAerodynamicForceBroken );
}
