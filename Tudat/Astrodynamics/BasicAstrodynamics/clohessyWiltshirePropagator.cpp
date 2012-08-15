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
 *      120810    E. Dekens         File created.
 *
 *    References
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications. Microcosm Press, 2001.
 */

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/clohessyWiltshirePropagator.h"

namespace tudat
{
namespace basic_astrodynamics
{

// Declare basic mathematical operations.
using std::sqrt;
using std::cos;
using std::sin;

//! Propagate linearized relative motion.
Eigen::VectorXd propagateClohessyWiltshire( const Eigen::VectorXd& initialState,
                                            const double propagationDuration,
                                            const double centralBodyGravitationalParameter,
                                            const double referenceOrbitRadius )
{
    // Initialize final state vector.
    Eigen::VectorXd finalState( 6 );

    // Calculate mean angular motion of reference orbit.
    const double meanAngularMotion =
            sqrt( centralBodyGravitationalParameter
                  / ( referenceOrbitRadius * referenceOrbitRadius * referenceOrbitRadius ) );

    // Calculate integration constants of Clohessy-Wiltshire equations.
    const double integrationConstant1 = 4.0 * initialState( 0 )
            + 2.0 * initialState( 4 ) / meanAngularMotion;
    const double integrationConstant2 = initialState( 1 )
            - 2.0 * initialState( 3 ) / meanAngularMotion;
    const double integrationConstant3 = 3.0 * initialState( 0 )
            + 2.0 * initialState( 4 ) / meanAngularMotion;
    const double integrationConstant4 = initialState( 3 ) / meanAngularMotion;
    const double integrationConstant5 = initialState( 2 );
    const double integrationConstant6 = initialState( 5 ) / meanAngularMotion;

    // Calculate dynamical terms of Clohessy-Wiltshire equations.
    const double cosineTerm = cos( meanAngularMotion * propagationDuration );
    const double sineTerm = sin( meanAngularMotion * propagationDuration );

    // Calculate Cartesian position components with Clohessy-Wilthsire equations as given by
    // by Vallado [2001, p. 382].
    finalState( 0 ) =  integrationConstant1 - integrationConstant3 * cosineTerm
            + integrationConstant4 * sineTerm;
    finalState( 1 ) = integrationConstant2
            - 1.5 * integrationConstant1 * meanAngularMotion * propagationDuration
            + 2.0 * integrationConstant3 * sineTerm
            + 2.0 * integrationConstant4 * cosineTerm;
    finalState( 2 ) = integrationConstant5 * cosineTerm + integrationConstant6 * sineTerm;

    // Calculate Cartesian velocity components with Clohessy-Wilthsire equations as given by
    // by Vallado [2001, p. 382].
    finalState( 3 ) = meanAngularMotion * integrationConstant3 * sineTerm
            + meanAngularMotion * integrationConstant4 * cosineTerm;
    finalState( 4 ) = - 1.5 * meanAngularMotion * integrationConstant1
            + 2.0 * meanAngularMotion * integrationConstant3 * cosineTerm
            - 2.0 * meanAngularMotion * integrationConstant4 * sineTerm;
    finalState( 5 ) = - meanAngularMotion * integrationConstant5 * sineTerm
            + meanAngularMotion * integrationConstant6 * cosineTerm;

    // Return final state.
    return finalState;
}

} // namespace basic_astrodynamics
} // namespace tudat
