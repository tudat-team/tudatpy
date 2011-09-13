/*! \file sphericalHarmonicsGravityField.h
 *    Header file that defines the spherical harmonics gravity field model
 *    included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 8
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
 *    Last modified     : 5 August, 2011
 *
 *    References
 *      Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. Revisiting
 *          Spacetrack Report #3: Rev 1, Proceedings of the AIAA/AAS Astro-
 *          dynamics Specialist Conference. Keystone, CO, 2006.
 *
 *    Notes
 *      The coefficients J2, J3, and J4 have been hardcoded for now, but in
 *      future these variables should be removed and the data should be stored
 *      in data files.
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
 *      110310    K. Kumar          Changed naming from Laplacian to gradient
 *                                  tensor.
 *      110805    K. Kumar          Added predefined functionality with WGS-72
 *                                  and WGS-84 predefined Earth gravity fields.
 */

#ifndef SPHERICALHARMONICSGRAVITYFIELD_H
#define SPHERICALHARMONICSGRAVITYFIELD_H

// Include statements.
#include "Astrodynamics/EnvironmentModels/GravityFieldModel/gravityFieldModel.h"

//! SphericalHarmonicsGravityField class.
/*!
 * Spherical harmonics gravity field model class included in Tudat.
 */
class SphericalHarmonicsGravityField : public GravityFieldModel
{
public:

    //! Bodies with predefined spherical harmonics gravity fields.
    /*!
     * Bodies with predefined spherical harmonics gravity fields.
     */
    enum BodiesWithPredefinedSphericalHarmonicsGravityFields
    {
        earthWorldGeodeticSystem72,
        earthWorldGeodeticSystem84
    };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SphericalHarmonicsGravityField( ):degreeOfExpansion_( -0 ), orderOfExpansion_( -0 ),
        referenceRadius_( -0.0 ), j2SphericalHarmonicsGravityFieldCoefficient_( -0.0 ),
        j3SphericalHarmonicsGravityFieldCoefficient_( -0.0 ),
        j4SphericalHarmonicsGravityFieldCoefficient_( -0.0 ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~SphericalHarmonicsGravityField( ){ }

    //! Set predefined spherical harmonics gravity field settings.
    /*!
     * Sets predefined spherical harmonics gravity field settings.
     * \param bodyWithPredefinedSphericalHarmonicsGravityField Body with
     *          predefined spherical harmonics gravity field.
     */
    void setPredefinedSphericalHarmonicsGravityFieldSettings(
        BodiesWithPredefinedSphericalHarmonicsGravityFields
        bodyWithPredefinedSphericalHarmonicsGravityField );

    //! Set the reference radius.
    /*!
     * Define the reference radius used for the spherical harmonics expansion
     * in meters.
     *  \param referenceRadius Reference radius.
     */
    void setReferenceRadius( const double& referenceRadius ){ referenceRadius_ = referenceRadius; }

    //! Set degree of spherical harmonics gravity field expansion.
    /*!
     * This function sets the degree of the spherical harmonics
     * gravity field expansion.
     * \param degreeOfExpansion Degree of spherical harmonics expansion.
     */
    void setDegreeOfExpansion( const unsigned int& degreeOfExpansion )
        { degreeOfExpansion_ = degreeOfExpansion; }

    //! Set order of spherical harmonics gravity field expansion.
    /*!
     * This function sets the order of the spherical harmonics
     * gravity field expansion.
     * \param orderOfExpansion Order of spherical harmonics expansion.
     */
    void setOrderOfExpansion( const unsigned int& orderOfExpansion )
        { orderOfExpansion_ = orderOfExpansion; }

    //! Get the reference radius.
    /*!
     * Return the reference radius used for the spherical harmonics expansion
     * in meters.
     * \return Reference radius.
     */
    double getReferenceRadius( ){ return referenceRadius_; }

    //! Get degree of spherical harmonics gravity field expansion.
    /*!
     * This function gets the degree of the spherical harmonics
     * gravity field expansion.
     * \return Degree of spherical harmonics expansion.
     */
    double getDegreeOfExpansion( ){ return degreeOfExpansion_; }

    //! Get order of spherical harmonics gravity field expansion.
    /*!
     * This function gets the order of the spherical harmonics
     * gravity field expansion.
     * \return Order of spherical harmonics expansion.
     */
    double getOrderOfExpansion( ){ return orderOfExpansion_; }

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

    //! Get gradient tensor of the gravitational potential.
    /*!
     * Get the value of the gradient tensor of the gravitational potential
     * expressed in spherical harmonics for the given position.
     * \param pointerToPosition Position given as a pointer to a
     *          CartesianPositionElements object.
     * \return Gradient tensor of gravitational potential.
     */
    Matrix3d getGradientTensorOfPotential( CartesianPositionElements*
                                           pointerToPosition );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param sphericalHarmonicsGravityField Spherical harmonics gravity field.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     SphericalHarmonicsGravityField&
                                     sphericalHarmonicsGravityField );

protected:

    //! Degree of spherical harmonics expansion.
    /*!
     * Degree of spherical harmonics expansion.
     */
    unsigned int degreeOfExpansion_;

    //! Order of spherical harmonics expansion.
    /*!
     * Order of spherical harmonics expansion.
     */
    unsigned int orderOfExpansion_;

    //! Reference radius.
    /*!
     * The reference radius used for the spherical harmonics expansion in
     * meters.
     */
    double referenceRadius_;

private:

    //! J2 spherical harmonics gravity field coefficient.
    /*!
     * J2 spherical harmonics gravity field coefficient.
     */
    double j2SphericalHarmonicsGravityFieldCoefficient_;

    //! J3 spherical harmonics gravity field coefficient.
    /*!
     * J3 spherical harmonics gravity field coefficient.
     */
    double j3SphericalHarmonicsGravityFieldCoefficient_;

    //! J4 spherical harmonics gravity field coefficient.
    /*!
     * J4 spherical harmonics gravity field coefficient.
     */
    double j4SphericalHarmonicsGravityFieldCoefficient_;
};

#endif // SPHERICALHARMONICSGRAVITYFIELD_H

// End of file.
