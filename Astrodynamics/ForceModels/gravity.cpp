/*! \file gravity.cpp
 *    Source file that defines the gravity force model included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 17 September, 2010
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      YYMMDD    author              comment
 *      100917    K. Kumar            File created.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Include statements modified and
 *                                    fileheader updated.
 */

// Include statements.
#include "gravity.h"

//! Default constructor.
Gravity::Gravity()
{
}

//! Default destructor.
Gravity::~Gravity()
{
}

//! Set body for gravity field expansion.
void Gravity::setBody( CelestialBody& celestialBody )
{
    // Set celestial body.
    celestialBody_ = celestialBody;
}

//! Compute state derivatives for gravity field expansion.
void Gravity::computeStateDerivatives( Vector& stateVector,
                              Vector& stateDerivativeVector )
{
    // Declare local variables.
    // PLACEHOLDER
    double stateVectorNorm;
    double stateVectorNormCubed;

    // Compute norm of state vector.
    // PLACEHOLDER
    stateVectorNorm = sqrt( stateVector( 0 ) * stateVector( 0 )
                      + stateVector( 1 ) * stateVector( 1 )
                      + stateVector( 2 ) * stateVector( 2 ) );

    // Set state derivative vector to size of state vector and fill with zeros.
    // PLACEHOLDER
    stateDerivativeVector.setZero( stateVector.rows() );

    // Compute cube of norm of state vector.
    // PLACEHOLDER
    stateVectorNormCubed = stateVectorNorm * stateVectorNorm * stateVectorNorm;

    // Compute state derivatives.
    stateDerivativeVector( 0 ) = stateVector( 3 );
    stateDerivativeVector( 1 ) = stateVector( 4 );
    stateDerivativeVector( 2 ) = stateVector( 5 );
    // PLACEHOLDER
    stateDerivativeVector( 3 ) = -celestialBody_.getGravitationalParameter()
                                 * stateVector( 0 )
                                 / stateVectorNormCubed;
    // PLACEHOLDER
    stateDerivativeVector( 4 ) = -celestialBody_.getGravitationalParameter()
                                 * stateVector( 1 )
                                 / stateVectorNormCubed;
    // PLACEHOLDER
    stateDerivativeVector( 5 ) = -celestialBody_.getGravitationalParameter()
                                 * stateVector( 2 )
                                 / stateVectorNormCubed;
}

// End of file.
