/*! \file sphericalHarmonicsGravityField.h
 *    Header file that defines the spherical harmonics gravity field model
 *    included in Tudat.
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
 *    Date created      : 16 November, 2010
 *    Last modified     : 4 February, 2011
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
 *      YYMMDD    Author            Comment
 *      101116    K. Kumar          File created.
 *      101117    K. Kumar          Added getPotential() and
 *                                  getGradientOfPotential functions.
 *      101214    K. Kumar          Updated getGradientOfPotential() and
 *                                  getLaplacianOfPotential().
 *      110106    K. Kumar          Added order and degree of expansion
 *                                  variables and set/get functions.
 *      110202    K. Kumar          Updated code to make use of the
 *                                  CartesianPositionElements class.
 *      110204    K. Kumar          Removed "vector" from naming.
 */

#ifndef SPHERICALHARMONICSGRAVITYFIELD_H
#define SPHERICALHARMONICSGRAVITYFIELD_H

// Include statements.
#include "gravityFieldModel.h"

//! SphericalHarmonicsGravityField class.
/*!
 * Spherical harmonics gravity field model class included in Tudat.
 */
class SphericalHarmonicsGravityField : public GravityFieldModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SphericalHarmonicsGravityField( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SphericalHarmonicsGravityField( );

    //! Set the reference radius.
    /*!
     * Define the reference radius used for the spherical harmonics expansion
     * in meters.
     *  \param referenceRadius Reference radius.
     */
    void setReferenceRadius( const double& referenceRadius );

    //! Set degree of spherical harmonics gravity field expansion.
    /*!
     * This function sets the degree of the spherical harmonics
     * gravity field expansion.
     * \param degreeOfExpansion Degree of spherical harmonics expansion.
     */
    void setDegreeOfExpansion( const int& degreeOfExpansion );

    //! Set order of spherical harmonics gravity field expansion.
    /*!
     * This function sets the order of the spherical harmonics
     * gravity field expansion.
     * \param orderOfExpansion Order of spherical harmonics expansion.
     */
    void setOrderOfExpansion( const int& orderOfExpansion );

    //! Get the reference radius.
    /*!
     * Return the reference radius used for the spherical harmonics expansion
     * in meters.
     * \return Reference radius.
     */
    double getReferenceRadius( );

    //! Get degree of spherical harmonics gravity field expansion.
    /*!
     * This function gets the degree of the spherical harmonics
     * gravity field expansion.
     * \return Degree of spherical harmonics expansion.
     */
    double getDegreeOfExpansion( );

    //! Get order of spherical harmonics gravity field expansion.
    /*!
     * This function gets the order of the spherical harmonics
     * gravity field expansion.
     * \return Order of spherical harmonics expansion.
     */
    double getOrderOfExpansion( );

    //! Get the gravitational potential.
    /*!
     * Get the value of the gravitational potential, expressed in spherical
     * harmonics for the given position.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Gravitational potential.
     */
    double getPotential( CartesianPositionElements* pointerToPosition );

    //! Get the gradient of the gravitational potential.
    /*!
     * Get the value of the gradient of the gravitational potential, expressed
     * in spherical harmonics for the given position.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Gradient of gravitational potential.
     */
    Vector3d getGradientOfPotential( CartesianPositionElements*
                                     pointerToPosition );

    //! Get the Laplacian of the gravitational potential.
    /*!
     * Get the value of the Laplacian of the gravitational potential
     * expressed in spherical harmonics for the given position.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Laplacian of gravitational potential.
     */
    Matrix3d getLaplacianOfPotential( CartesianPositionElements*
                                      pointerToPosition );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param pointerToSphericalHarmonicsGravityField Pointer to
     *          spherical harmonics gravity field.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     SphericalHarmonicsGravityField*
                                     pointerToSphericalHarmonicsGravityField );

protected:

private:

    //! Degree of spherical harmonics expansion.
    /*!
     * Degree of spherical harmonics expansion.
     */
    int degreeOfExpansion_;

    //! Order of spherical harmonics expansion.
    /*!
     * Order of spherical harmonics expansion.
     */
    int orderOfExpansion_;

    //! Reference radius.
    /*!
     * The reference radius used for the spherical harmonics expansion in
     * meters.
     */
    double referenceRadius_;
};

#endif // SPHERICALHARMONICSGRAVITYFIELD_H

// End of file.
