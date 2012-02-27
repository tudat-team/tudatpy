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
 *      120220    D. Dirkx          Moved moment due to force test to separate file after split
 *                                  of moment model base class.
 *
 *    References
 *
 */

// Include statements.
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <limits>
#include "Tudat/Astrodynamics/BasicAstrodynamics/momentDueToForceModel.h"

//! Test implementation of moment due to force.
int main( )
{
    using std::cerr;
    using std::endl;
    using std::fabs;
    using namespace tudat;

    bool isMomentDueToForceBroken = true;

    // Initialise objects.
    Eigen::Vector3d forceCoefficients;
    forceCoefficients( 0 ) = 1.1;
    forceCoefficients( 1 ) = 1.2;
    forceCoefficients( 2 ) = 1.3;

    Eigen::Vector3d momentArm;
    momentArm( 0 ) = -2.0;
    momentArm( 1 ) = 0.0;
    momentArm( 2 ) = 0.0;

    double dynamicPressure = 50.0;
    double referenceArea = 2.2;

    Eigen::Vector3d intermediateExpectedForce
            = forceCoefficients * dynamicPressure * referenceArea;

    // Calculate moment.
    MomentDueToForceModel momentDueToForce = MomentDueToForceModel();
    momentDueToForce.computeMoment( intermediateExpectedForce, momentArm );
    Eigen::Vector3d moment = momentDueToForce.getMomentDueToForce( );

    // Calculate expected moment.
    Eigen::Vector3d expectedMomentDueToForce;
    expectedMomentDueToForce( 0 ) = 0.0;
    expectedMomentDueToForce( 1 ) = -1.0 * momentArm( 0 ) * intermediateExpectedForce( 2 );
    expectedMomentDueToForce( 2 ) = momentArm( 0 ) * intermediateExpectedForce( 1 );

    // Error in calculation.
    Eigen::Vector3d errorInMoment;
    errorInMoment( 0 ) = fabs( ( expectedMomentDueToForce( 0 ) - moment( 0 ) ) / moment( 0 ) );
    errorInMoment( 1 ) = fabs( ( expectedMomentDueToForce( 1 ) - moment( 1 ) ) / moment( 1 ) );
    errorInMoment( 2 ) = fabs( ( expectedMomentDueToForce( 2 ) - moment( 2 ) ) / moment( 2 ) ) ;

    if ( errorInMoment( 0 ) > std::numeric_limits< double >::epsilon( ) ||
         errorInMoment( 1 ) > std::numeric_limits< double >::epsilon( ) ||
         errorInMoment( 2 ) > std::numeric_limits< double >::epsilon( ) )
    {
        isMomentDueToForceBroken = true;
        cerr << "Test momentDueToForce failed." << endl;
    }
}
