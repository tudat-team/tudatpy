/*! \file unitTestAerodynamicMomentAndAerodynamicForce.cpp
*    This source file contains the unit test for aerodynamic forces model
*    included in Tudat.
*
*    Path              : /Astrodynamics/ForceModels/
*    Version           : 3
*    Check status      : Checked
*
*    Checker           : F. M. Engelen
*    Affiliation       : Delft University of Technology
*    E-mail address    : F.M.Engelen@student.tudelft.nl
*
*    Checker           : D. Dirkx
*    Affiliation       : Delft University of Technology
*    E-mail address    : D.Dirkx@tudelft.nl
*
*    Date created      : 22 June, 2011
*    Last modified     : 24 August, 2011
*
*    References
*
*    Notes
*
*    Copyright (c) 2010-2011 Delft University of Technology.
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
*/

// Include statements.
#include "Astrodynamics/ForceModels/unitTestAerodynamicMomentAndAerodynamicForce.h"
#include "Astrodynamics/ForceModels/Aerothermodynamics/aerodynamicCoefficientInterface.h"
#include "Astrodynamics/ForceModels/aerodynamicForce.h"
#include "Astrodynamics/MomentModels/aerodynamicMoment.h"
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/LinearAlgebra/linearAlgebra.h"
#include "Mathematics/unitConversions.h"

// Using statements.
using Eigen::Quaternion;
using Eigen::AngleAxisd;
using std::cerr;
using std::endl;
using std::fabs;
using mathematics::MACHINE_PRECISION_DOUBLES;

bool unit_tests::testAerodynamicMomentAndAerodynamicForce( )
{
    bool isAerodynamicMomentBroken = false;
    bool isAerodynamicForceBroken = false;

    // Initialise objects.
    Vector3d forceCoefficients;
    forceCoefficients(0) = 1.1;
    forceCoefficients(1) = 1.2;
    forceCoefficients(2) = 1.3;

    Vector3d momentCoefficients;
    momentCoefficients(0) = 0.0;
    momentCoefficients(1) = 1.0;
    momentCoefficients(2) = 0.0;

    double dynamicPressure = 50.0;
    double referenceArea = 2.2;
    double referenceLength = 3.2;

    AerodynamicCoefficientInterface aerodynamicCoefficientInterface;
    aerodynamicCoefficientInterface.setCurrentForceCoefficients( forceCoefficients );
    aerodynamicCoefficientInterface.setCurrentMomentCoefficients( momentCoefficients );
    aerodynamicCoefficientInterface.setReferenceArea( referenceArea );
    aerodynamicCoefficientInterface.setReferenceLength( referenceLength );

    Vector3d momentArm;
    momentArm(0) = -2.0;
    momentArm(1) = 0.0;
    momentArm(2) = 0.0;

    AerodynamicForce aerodynamicForce;
    aerodynamicForce.setDynamicPressure( dynamicPressure );
    aerodynamicForce.setAerodynamicCoefficientInterface( &aerodynamicCoefficientInterface );

    Vector3d force;

    // Test 1. Test the force model without setting a transformation.

    // Calculated force.
    State dummyState;
    aerodynamicForce.computeForce( &dummyState );
    force = aerodynamicForce.getForce( );

    // Expected force.
    Vector3d expectedForce = forceCoefficients * dynamicPressure * referenceArea;

    // Error in calculation.
    Vector3d errorInForce;
    errorInForce( 0 ) = fabs( expectedForce( 0 ) - force( 0 ) );
    errorInForce( 1 ) = fabs( expectedForce( 1 ) - force( 1 ) );
    errorInForce( 2 ) = fabs( expectedForce( 2 ) - force( 2 ) );

    if ( MACHINE_PRECISION_DOUBLES < errorInForce.sum( ) / 3.0 /
         ( dynamicPressure * referenceArea ) )
    {
        isAerodynamicForceBroken = true;
        cerr << "Test 1 of unitTestAerodynamicMomentAndAerodynamicForce failed" << endl;
    }

     // Test 2: Check the moment model without extra moment due to the force.
    AerodynamicMoment aerodynamicMoment;
    aerodynamicMoment.setDynamicPressure( dynamicPressure );
    aerodynamicMoment.setAerodynamicCoefficientInterface( &aerodynamicCoefficientInterface );

    Vector3d moment;
    Vector3d expectedMoment;
    Vector3d errorInMoment;

    // Calculate moment.
    aerodynamicMoment.computeMoment( &dummyState );
    moment = aerodynamicMoment.getMoment( );

    // Calculate expected moment.
    expectedMoment = dynamicPressure * referenceArea * referenceLength * momentCoefficients;

    // Error in calculation.
    errorInMoment( 0 ) = fabs( expectedMoment( 0 ) - moment( 0 ) );
    errorInMoment( 1 ) = fabs( expectedMoment( 1 ) - moment( 1 ) );
    errorInMoment( 2 ) = fabs( expectedMoment( 2 ) - moment( 2 ) );

    if ( MACHINE_PRECISION_DOUBLES < errorInMoment.sum( ) / 3.0 /
         ( dynamicPressure * referenceArea * referenceLength ) )
    {
        isAerodynamicMomentBroken = true;
        cerr << "Test 2 of unitTestAerodynamicMomentAndAerodynamicForce failed" << endl;
    }

    // Test 3: Check the moment model when the forcemodel is set and a constant arm is provided.
    aerodynamicMoment.setForceModel( &aerodynamicForce );
    aerodynamicMoment.setForceApplicationArm( momentArm );

    Vector3d intermediateExpectedForce =  ( forceCoefficients * dynamicPressure * referenceArea );

    // Calculate moment.
    aerodynamicMoment.computeMoment( &dummyState );
    moment = aerodynamicMoment.getMoment( );

    // Calculate expected moment.
    Vector3d expectedMomentDueToForce;
    expectedMomentDueToForce( 0 ) = 0.0;
    expectedMomentDueToForce( 1 ) = -1.0 * momentArm( 0 ) * intermediateExpectedForce( 2 );
    expectedMomentDueToForce( 2 ) = momentArm( 0 ) * intermediateExpectedForce( 1 );

    expectedMoment = dynamicPressure * referenceArea * referenceLength * momentCoefficients +
                     expectedMomentDueToForce ;

    // Error in calculation.
    errorInMoment( 0 ) = fabs( expectedMoment( 0 ) - moment( 0 ) );
    errorInMoment( 1 ) = fabs( expectedMoment( 1 ) - moment( 1 ) );
    errorInMoment( 2 ) = fabs( expectedMoment( 2 ) - moment( 2 ) );

    if ( MACHINE_PRECISION_DOUBLES < errorInMoment.sum( ) / 3.0 /
         ( dynamicPressure * referenceArea * referenceLength ) )
    {
        isAerodynamicMomentBroken = true;
        cerr << "Test 3 of unitTestAerodynamicMomentAndAerodynamicForce failed" << endl;
    }

    return ( isAerodynamicMomentBroken || isAerodynamicForceBroken );
}

// End of file.
