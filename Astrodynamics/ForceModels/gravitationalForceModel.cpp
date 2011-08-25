/*! \file gravitationalForceModel.cpp
 *    Source file that defines the gravitational force model included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 17 September, 2010
 *    Last modified     : 15 August, 2011
 *
 *    References
 *
 *    Notes
 *    The gradient of potential which is supplied needs to be in cartesian coordinates,
 *    and not in spherical coordinates.
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
 *      YYMMDD    Author            Comment
 *      100917    K. Kumar          File created.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Include statements modified and
 *                                  fileheader updated.
 *      110113    K. Kumar          Updated computeStateDerivatives().
 *      110119    K. Kumar          Changed computeStateDerivatives() to
 *                                  computeForce().
 *      110202    K. Kumar          Updated code to make use of the State and
 *                                  CartesianPositionElements classes.
 *      110815    K. Kumar          Changed filename and class name; changed
 *                                  computeForce() function and added
 *                                  setMass() function.
 *      110815    F.M. Engelen      Updated comments.
 */

// Include statements.
#include "gravitationalForceModel.h"

//! Default constructor.
GravitationalForceModel::GravitationalForceModel( )
{
}

//! Default destructor.
GravitationalForceModel::~GravitationalForceModel( )
{
}

//! Set body subject to force.
void GravitationalForceModel::setBodySubjectToForce(
    Body* pointerToBodySubjectToForce )
{
    pointerToBodySubjectToForce_ = pointerToBodySubjectToForce;
}

//! Set body for gravity field expansion.
void GravitationalForceModel::setGravitationalBody( CelestialBody* pointerToCelestialBody )
{
    // Set celestial body.
    pointerToCelestialBody_ = pointerToCelestialBody;

    // Set pointer to gravity field model to gravity field model stored in
    // CelestialBody object.
    pointerToGravityFieldModel_ = pointerToCelestialBody_
            ->getGravityFieldModel( );
}

//! Compute force due to gravity field.
void GravitationalForceModel::computeForce( State* pointerToState )
{
    // Set Cartesian position elements state.
    cartesianPositionElements_.state = pointerToState->state.segment( 0, 3 );

    // Compute forces in Newtons using gradient of potential of gravity
    // expansion and the mass of the body subjected to the force.
    force_ = pointerToGravityFieldModel_->getGradientOfPotential(
                &cartesianPositionElements_ )
            * pointerToBodySubjectToForce_->getMass( );
}

// End of file.
