/*! \file gravitationalForceModel.h
 *    Header file that defines the gravitational force model included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 9
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
 *    Date created      : 16 September, 2010
 *    Last modified     : 15 August, 2011
 *
 *    References
 *
 *    Notes
 *      The mass and state of the body on which the force acts should be
 *      transferred to the Body class.
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
 *      100916    K. Kumar          File created.
 *      100916    K. Kumar          Filename modified.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor corrections to include statements and comments.
 *      110113    K. Kumar          Changed setBody() argument to pointer; added pointer to
 *                                  GravityFieldModel.
 *      110119    K. Kumar          Changed computeStateDerivatives() to computeForce().
 *      110202    K. Kumar          Updated code to make use of the State and
 *                                  CartesianPositionElements classes.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110815    K. Kumar          Changed filename and class name; changed computeForce()
 *                                  function and added setMass() function.
 */

#ifndef GRAVITATIONALFORCEMODEL_H
#define GRAVITATIONALFORCEMODEL_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/EnvironmentModels/gravityFieldModel.h"
#include "Astrodynamics/ForceModels/forceModel.h"

//! Gravitational force model class.
/*!
 * Class containing the gravitational force model.
 */
class GravitationalForceModel : public ForceModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GravitationalForceModel( ) : pointerToBodySubjectToForce_( NULL ),
        pointerToCelestialBody_( NULL ), pointerToGravityFieldModel_( NULL ) { }

    //! Set body subject to force.
    /*!
     * Sets body subject to force.
     * \param pointerToBodySubjectToForce Pointer to body subject to force.
     */
    void setBodySubjectToForce( Body* pointerToBodySubjectToForce )
    { pointerToBodySubjectToForce_ = pointerToBodySubjectToForce; }

    //! Set body for gravity field expansion.
    /*!
     * Sets the body for gravity field expansion.
     * \param pointerToCelestialBody Celestial body which is set.
     */
    void setGravitationalBody( CelestialBody* pointerToCelestialBody )
    {
        pointerToCelestialBody_ = pointerToCelestialBody;
        pointerToGravityFieldModel_ = pointerToCelestialBody_->getGravityFieldModel( );
    }

    //! Compute force due to gravity field.
    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class.
     */
    void computeForce( State* pointerToState )
    {
        cartesianPositionElements_.state = pointerToState->state.segment( 0, 3 );
        force_ = pointerToGravityFieldModel_->getGradientOfPotential( &cartesianPositionElements_ )
                * pointerToBodySubjectToForce_->getMass( );
    }

protected:

private:

    //! Body subject to force.
    /*!
     * Body subject to force.
     */
    Body* pointerToBodySubjectToForce_;

    //! Pointer to celestial body.
    /*!
     * Pointer to celestial body.
     */
    CelestialBody* pointerToCelestialBody_;

    //! Pointer to gravity field model.
    /*!
     * Pointer to gravity field model.
     */
    GravityFieldModel* pointerToGravityFieldModel_;

    //! Cartesian position elements.
    /*!
     * Cartesian position elements.
     */
    CartesianPositionElements cartesianPositionElements_;
};

#endif // GRAVITATIONALFORCEMODEL_H

// End of file.
