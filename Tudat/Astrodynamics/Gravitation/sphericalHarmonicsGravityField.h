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
 *      101116    K. Kumar          File created.
 *      101117    K. Kumar          Added getPotential( ) and getGradientOfPotential functions.
 *      101214    K. Kumar          Updated getGradientOfPotential( ) and getLaplacianOfPotential( ).
 *      110106    K. Kumar          Added order and degree of expansion variables and set/get
 *                                  functions.
 *      110202    K. Kumar          Updated code to make use of the CartesianPositionElements
 *                                  class.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *      110805    K. Kumar          Added predefined functionality with WGS-72 and WGS-84
 *                                  predefined Earth gravity fields.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *      Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. Revisiting Spacetrack Report #3:
 *          Rev 1, Proceedings of the AIAA/AAS Astrodynamics Specialist Conference. Keystone, CO,
 *          2006.
 *
 *    The coefficients J2, J3, and J4 have been hardcoded for now, but in future these variables
 *    should be removed and the data should be stored in data files.
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H

#include <iostream>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{

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
    SphericalHarmonicsGravityField( const unsigned int degreeOfExpansion = 0,
                                    const unsigned int orderOfExpansion = 0,
                                    const double referenceRadius = 0.0 )
        : degreeOfExpansion_( degreeOfExpansion ),
          orderOfExpansion_( orderOfExpansion ),
          referenceRadius_( referenceRadius ),
          j2SphericalHarmonicsGravityFieldCoefficient_( TUDAT_NAN ),
          j3SphericalHarmonicsGravityFieldCoefficient_( TUDAT_NAN ),
          j4SphericalHarmonicsGravityFieldCoefficient_( TUDAT_NAN )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~SphericalHarmonicsGravityField( ) { }

    //! Set predefined spherical harmonics gravity field settings.
    /*!
     * Sets predefined spherical harmonics gravity field settings.
     * \param bodyWithPredefinedSphericalHarmonicsGravityField Body with
     *          predefined spherical harmonics gravity field.
     */
    void setPredefinedSphericalHarmonicsGravityFieldSettings(
        BodiesWithPredefinedSphericalHarmonicsGravityFields
        bodyWithPredefinedSphericalHarmonicsGravityField );

    //! Get the reference radius.
    /*!
     * Returns the reference radius used for the spherical harmonics expansion in meters.
     * \return Reference radius.
     */
    double getReferenceRadius( ) { return referenceRadius_; }

    //! Get degree of spherical harmonics gravity field expansion.
    /*!
     * Returns the degree of the spherical harmonics gravity field expansion.
     * \return Degree of spherical harmonics expansion.
     */
    double getDegreeOfExpansion( ) { return degreeOfExpansion_; }

    //! Get order of spherical harmonics gravity field expansion.
    /*!
     * Returns the order of the spherical harmonics gravity field expansion.
     * \return Order of spherical harmonics expansion.
     */
    double getOrderOfExpansion( ) { return orderOfExpansion_; }

    //! Get the gravitational potential.
    /*!
     * Returns the value of the gravitational potential, expressed in spherical harmonics for the
     * given position.
     * \param position Position at which potential is to be determined
     * \return Gravitational potential.
     */
    double getPotential( const Eigen::Vector3d& position )
    {
        relativePosition_ = position - positionOfOrigin_;
        return gravitationalParameter_ / relativePosition_.norm( );
    }

    //! Get the gradient of the gravitational potential.
    /*!
     * Returns the value of the gradient of the gravitational potential, expressed in spherical
     * harmonics for the given position.
     * \param position Position at which gradient of potential is to be determined
     * \return Gradient of gravitational potential.
     */
    Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& position )
    {
        relativePosition_ = position - positionOfOrigin_;
        return -gravitationalParameter_ * relativePosition_
                / pow( relativePosition_.norm( ), 3.0 );
    }

    //! Get gradient tensor of the gravitational potential.
    /*!
     * Returns the value of the gradient tensor of the gravitational potential expressed in
     * spherical harmonics for the given position.
     * \param position Position at which gradient tensor of potential is to be determined
     * \return Gradient tensor of gravitational potential.
     */
    Eigen::Matrix3d getGradientTensorOfPotential( const Eigen::Vector3d& position );

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

} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H
