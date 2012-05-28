/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110622    F.M. Engelen      Creation of code.
 *      110822    D. Dirkx          Removed no longer necessary unit tests.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      120220    D. Dirkx          Moved moment due to force test to separate file after split
 *                                  of moment model base class.
 *
 *    References
 *
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/momentDueToForceModel.h"

//! Test implementation of moment due to force.
int main( )
{
    using std::cerr;
    using std::endl;
    using std::fabs;
    using namespace tudat;
    using namespace tudat::astrodynamics::moment_models;

    bool isMomentDueToForceBroken = false;

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
    errorInMoment( 2 ) = fabs( ( expectedMomentDueToForce( 2 ) - moment( 2 ) ) / moment( 2 ) );

    if ( errorInMoment( 0 ) > std::numeric_limits< double >::epsilon( ) ||
         errorInMoment( 1 ) > std::numeric_limits< double >::epsilon( ) ||
         errorInMoment( 2 ) > std::numeric_limits< double >::epsilon( ) )
    {
        isMomentDueToForceBroken = true;
        cerr << "Test momentDueToForce using class failed." << endl;
    }

    moment = computeMomentDueToForce( intermediateExpectedForce, momentArm );

    // Error in calculation.
    errorInMoment( 0 ) = fabs( ( expectedMomentDueToForce( 0 ) - moment( 0 ) ) / moment( 0 ) );
    errorInMoment( 1 ) = fabs( ( expectedMomentDueToForce( 1 ) - moment( 1 ) ) / moment( 1 ) );
    errorInMoment( 2 ) = fabs( ( expectedMomentDueToForce( 2 ) - moment( 2 ) ) / moment( 2 ) );

    if ( errorInMoment( 0 ) > std::numeric_limits< double >::epsilon( ) ||
         errorInMoment( 1 ) > std::numeric_limits< double >::epsilon( ) ||
         errorInMoment( 2 ) > std::numeric_limits< double >::epsilon( ) )
    {
        isMomentDueToForceBroken = true;
        cerr << "Test momentDueToForce using class failed." << endl;
    }

    return isMomentDueToForceBroken;
}
