/*! \file gravity.cpp
 *    Source file that defines the gravity force model included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 5
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
 *    Last modified     : 19 January, 2011
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
 *      110113    K. Kumar            Updated computeStateDerivatives().
 *      110119    K. Kumar            Changed computeStateDerivatives() to
 *                                    computeForce().
 */

// Include statements.
#include "gravity.h"

//! Default constructor.
Gravity::Gravity( )
{
    // Initialize variables.
    forcePerUnitMass_.setZero( 3 );
}

//! Default destructor.
Gravity::~Gravity( )
{
}

//! Set body for gravity field expansion.
void Gravity::setBody( CelestialBody* celestialBody )
{
    // Set celestial body.
    celestialBody_ = celestialBody;
}

//! Compute force per unit mass for gravity field expansion.
VectorXd& Gravity::computeForce( VectorXd& stateVector )
{
    // Set pointer to gravity field model to gravity field model stored in
    // CelestialBody object.
    pointerToGravityFieldModel_ = celestialBody_->getGravityFieldModel( );

    // Compute forces per unit mass using gradient of potential of gravity
    // expansion.
    forcePerUnitMass_ = pointerToGravityFieldModel_
                        ->getGradientOfPotential( stateVector.segment( 0, 3 ) );

    // Return computed forces per unit mass.
    return forcePerUnitMass_;
}

// End of file.
