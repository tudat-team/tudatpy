/*! \file gravityFieldModel.h
 *    Header file that defines the gravity field model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 9
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 10 November, 2010
 *    Last modified     : 10 March, 2011
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
 *      101110    K. Kumar          File created.
 *      101116    K. Kumar          Changed filename and class name.
 *      101215    K. Kumar          Added virtual functions and missing
 *                                  Doxygen comments.
 *      101216    K. Kumar          Added set/get functions for origin.
 *      110107    K. Kumar          Removed reference radius get/set virtual
 *                                  functions.
 *      110113    K. Kumar          Updated arguments of get-functions with
 *                                  const.
 *      110202    K. Kumar          Updated code to use the
 *                                  CartesianPositionElements class.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient
 *                                  tensor.
 */

#ifndef GRAVITYFIELDMODEL_H
#define GRAVITYFIELDMODEL_H

// Include statements.
#include "Astrodynamics/EnvironmentModels/environmentModel.h"
#include "Astrodynamics/States/cartesianPositionElements.h"

//! GravityFieldModel class.
/*!
 * Gravity field model class included in Tudat.
 */
class GravityFieldModel : public EnvironmentModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GravityFieldModel( ): gravitationalParameter_( -0.0 ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~GravityFieldModel( ){ }

    //! Set the gravitational parameter.
    /*!
     * Define the gravitational parameter in meter^3 per second^2.
     * \param gravitationalParameter
     */
    void setGravitationalParameter( const double& gravitationalParameter )
        { gravitationalParameter_ = gravitationalParameter; }

    //! Set origin of gravity field.
    /*!
     * Set origin of gravity field.
     * \param pointerToPositionOfOrigin Position of origin given as a pointer to a
     *          CartesianPositionElements object.
     */
    void setOrigin( CartesianPositionElements* pointerToPositionOfOrigin )
        { positionOfOrigin_ = *pointerToPositionOfOrigin; }

    //! Get the gravitational parameter.
    /*!
     * Return the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    double getGravitationalParameter( ){ return gravitationalParameter_; }

    //! Get origin of gravity field.
    /*!
     * Get origin of gravity field.
     * \return Position of origin given as a pointer to a
     *          CartesianPositionElements object.
     */
    CartesianPositionElements* getOrigin( ){ return &positionOfOrigin_; }

    //! Get the potential.
    /*!
     * Returns the potential for the gravity field selected.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Potential.
     */
    virtual double getPotential( CartesianPositionElements*
                                 pointerToPosition ) = 0;

    //! Get the gradient of the potential.
    /*!
     * Returns the gradient of the potential for the gravity field selected.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Gradient of potential.
     */
    virtual Vector3d getGradientOfPotential(
            CartesianPositionElements* pointerToPosition ) = 0;

    //! Get gradient tensor of the potential.
    /*!
     * Returns the gradient tensor of the potential for the gravity field
     * selected.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Gradient tensor of potential.
     */
    virtual Matrix3d getGradientTensorOfPotential( CartesianPositionElements*
                                                   pointerToPosition ) = 0;

protected:

    //! Gravitational parameter.
    /*!
     * The gravitational parameter in meter^3 per second^2.
     */
    double gravitationalParameter_;

    //! Origin of gravity field.
    /*!
     * Origin of gravity field given as a CartesianPositionElements object.
     */
    CartesianPositionElements positionOfOrigin_;

    //! Relative position.
    /*!
     * Relative position given as a CartesianPositionElements object.
     */
    CartesianPositionElements relativePosition_;

private:
};

#endif // GRAVITYFIELDMODEL_H

// End of file.
