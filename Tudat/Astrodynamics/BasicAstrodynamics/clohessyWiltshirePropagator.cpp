/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications. Microcosm Press, 2001.
 *
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
Eigen::Vector6d propagateClohessyWiltshire(
        const Eigen::Vector6d& initialState,
        const double propagationDuration,
        const double centralBodyGravitationalParameter,
        const double referenceOrbitRadius )
{
    // Initialize final state vector.
    Eigen::Vector6d finalState;

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
