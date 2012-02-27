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
 *      100916    K. Kumar          File created.
 *      100916    K. Kumar          Filename modified.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor corrections to include statements and comments.
 *      110113    K. Kumar          Changed setBody( ) argument to pointer; added pointer to
 *                                  GravityFieldModel.
 *      110119    K. Kumar          Changed computeStateDerivatives( ) to computeForce( ).
 *      110202    K. Kumar          Updated code to make use of the State and
 *                                  CartesianPositionElements classes.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110815    K. Kumar          Changed filename and class name; changed computeForce( )
 *                                  function and added setMass( ) function.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// The mass and state of the body on which the force acts should be
// transferred to the Body class.
// 

#ifndef TUDAT_GRAVITATIONAL_FORCE_MODEL_H
#define TUDAT_GRAVITATIONAL_FORCE_MODEL_H

#include "Tudat/Astrodynamics/Bodies/body.h"
#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h"

namespace tudat
{

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
     * \param pointerToBodySubjectToForce Pointer to the body that on which the force is acting.
     * \param pointerToGravityFieldModel Pointer to gravity field model object that is causing
     * the gravitational force.
     */
    GravitationalForceModel( Body* pointerToBodySubjectToForce,
                             GravityFieldModel* pointerToGravityFieldModel ):
            pointerToBodySubjectToForce_( pointerToBodySubjectToForce ),
            pointerToGravityFieldModel_(  pointerToGravityFieldModel ) { }

    //! Get pointer to the body that on which the force is acting.
    /*!
     * Returns the pointer to the body that on which the force is acting.
     * \return Pointer to the body that on which the force is acting.
     */
    Body* getPointerToBodySubjectToForce( ) { return  pointerToBodySubjectToForce_; }

    //! Get pointer to the gravitational field model.
    /*!
     * Returns the pointer to the gravitational field model.
     * \return Pointer to the gravitational field model.
     */
    GravityFieldModel* getPointerToGravityFieldModel( )  { return pointerToGravityFieldModel_; }

    //! Compute force due to gravity field.
    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class.
     * \param time Time (or other independent variable).
     */
    void computeForce( State* pointerToState, double time );

protected:

private:

    //! Body subject to force.
    /*!
     * Body subject to force.
     */
    Body* pointerToBodySubjectToForce_;

    //! Pointer to gravity field model.
    /*!
     * Pointer to gravity field model.
     */
    GravityFieldModel* pointerToGravityFieldModel_;
};

} // namespace tudat

#endif // TUDAT_GRAVITATIONAL_FORCE_MODEL_H
