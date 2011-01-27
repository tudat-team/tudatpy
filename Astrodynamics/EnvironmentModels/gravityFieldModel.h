/*! \file gravityFieldModel.h
 *    Header file that defines the gravity field model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModel/
 *    Version           : 6
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
 *    Last modified     : 13 January, 2011
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
 *      101110    K. Kumar            File created.
 *      101116    K. Kumar            Changed filename and class name.
 *      101215    K. Kumar            Added virtual functions and missing
 *                                    Doxygen comments.
 *      101216    K. Kumar            Added set/get functions for origin.
 *      110107    K. Kumar            Removed reference radius get/set virtual
 *                                    functions.
 *      110113    K. Kumar            Updated arguments of get-functions with
 *                                    const.
 */

#ifndef GRAVITYFIELDMODEL_H
#define GRAVITYFIELDMODEL_H

// Include statements.
#include "environmentModel.h"

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
    GravityFieldModel( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~GravityFieldModel( );

    //! Set the gravitational parameter.
    /*!
     * Define the gravitational parameter in meter^3 per second^2.
     * \param gravitationalParameter
     */
    void setGravitationalParameter( const double& gravitationalParameter );

    //! Set origin of gravity field.
    /*!
     * Set origin of gravity field.
     * \param positionVectorOfOrigin Position vector of origin.
     */
    void setOrigin( Vector3d& positionVectorOfOrigin );

    //! Get the gravitational parameter.
    /*!
     * Return the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    double getGravitationalParameter( );

    //! Get origin of gravity field.
    /*!
     * Get origin of gravity field.
     * \return Position vector of origin.
     */
    Vector3d getOrigin( );

    //! Get the potential.
    /*!
     * Returns the potential for the gravity field selected.
     * \return Potential.
     */
    virtual double getPotential( const Vector3d& positionVector ) = 0;

    //! Get the gradient of the potential.
    /*!
     * Returns the gradient of the potential for the gravity field selected.
     * \return Gradient of potential.
     */
    virtual Vector3d getGradientOfPotential(
            const Vector3d& positionVector ) = 0;

    //! Get the Laplacian of the potential.
    /*!
     * Returns the Laplacian of the potential for the gravity field selected.
     * \return Laplacian of potential.
     */
    virtual Matrix3d getLaplacianOfPotential(
            const Vector3d& positionVector ) = 0;

protected:

    //! Gravitational parameter.
    /*!
     * The gravitational parameter in meter^3 per second^2.
     */
    double gravitationalParameter_;

    //! Origin of gravity field.
    /*!
     * Origin of gravity field.
     */
    Vector3d positionVectorOfOrigin_;

private:
};

#endif // GRAVITYFIELDMODEL_H

// End of file.
