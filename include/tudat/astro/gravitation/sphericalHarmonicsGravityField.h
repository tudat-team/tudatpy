/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. Revisiting Spacetrack Report #3:
 *          Rev 1, Proceedings of the AIAA/AAS astro Specialist Conference. Keystone, CO,
 *          2006.
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

namespace tudat
{

namespace gravitation
{

//! Function to calculate the gravitational potential from a spherical harmonic field expansion.
/*!
 *  Function to calculate the gravitational potential from a spherical harmonic field expansion.
 *  \param bodyFixedPosition Position of point at which potential is to be calculated wrt the
 *  massive body, in the frame in which the expansion is defined (typically body-fixed).
 *  \param gravitationalParameter Gravitational parameter of massive body.
 *  \param referenceRadius Reference radius of spherical harmonic field expansion.
 *  \param cosineCoefficients Cosine spherical harmonic coefficients (geodesy normalized).
 *  \param sineCoefficients Sine spherical harmonic coefficients (geodesy normalized).
 *  \param sphericalHarmonicsCache Cache object containing current values of trigonometric funtions of latitude anf longitude,
 *  as well as legendre polynomials at current state.
 *  \param minimumumDegree Maximum degree of spherical harmonic expansion.
 *  \param minimumumOrder Maximum order of spherical harmonic expansion.

 *  \return Gravitational potential at position defined by bodyFixedPosition
 */
double calculateSphericalHarmonicGravitationalPotential(
        const Eigen::Vector3d& bodyFixedPosition, const double gravitationalParameter,
        const double referenceRadius,
        const Eigen::MatrixXd& cosineCoefficients, const Eigen::MatrixXd& sineCoefficients,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const int minimumumDegree = 0, const int minimumumOrder = 0 );

//! Class to represent a spherical harmonic gravity field expansion.
/*!
 *  Class to represent a spherical harmonic gravity field expansion of a massive body with
 *  time-independent spherical harmonic gravity field coefficients.
 */
class SphericalHarmonicsGravityField: public GravityFieldModel
{
public:

    //! Class constructor.
    /*!
     *  Class constructor.
     *  \param gravitationalParameter Gravitational parameter of massive body
     *  \param referenceRadius Reference radius of spherical harmonic field expansion
     *  \param cosineCoefficients Cosine spherical harmonic coefficients (geodesy normalized)
     *  \param sineCoefficients Sine spherical harmonic coefficients (geodesy normalized)
     *  \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional).
     *  \param updateInertiaTensor Function that is to be called to update the inertia tensor (typicaly in Body class; default
     *  empty)
     */
    SphericalHarmonicsGravityField(
            const double gravitationalParameter,
            const double referenceRadius,
            const Eigen::MatrixXd& cosineCoefficients = Eigen::MatrixXd::Identity( 1, 1 ),
            const Eigen::MatrixXd& sineCoefficients = Eigen::MatrixXd::Zero( 1, 1 ),
            const std::string& fixedReferenceFrame = "",
            const std::function< void( ) > updateInertiaTensor = std::function< void( ) > ( ) )
        : GravityFieldModel( gravitationalParameter, updateInertiaTensor ), referenceRadius_( referenceRadius ),
          cosineCoefficients_( cosineCoefficients ), sineCoefficients_( sineCoefficients ),
          fixedReferenceFrame_( fixedReferenceFrame )
    {
        sphericalHarmonicsCache_ = std::make_shared< basic_mathematics::SphericalHarmonicsCache >( );
        sphericalHarmonicsCache_->resetMaximumDegreeAndOrder( cosineCoefficients_.rows( ) + 1,
                                                              cosineCoefficients_.cols( ) + 1 );
    }

    //! Virtual destructor.
    /*!
     *  Virtual destructor.
     */
    virtual ~SphericalHarmonicsGravityField( ) { }

    //! Function to get the reference radius.
    /*!
     *  Returns the reference radius used for the spherical harmonics expansion in meters.
     *  \return Reference radius of spherical harmonic field expansion
     */
    double getReferenceRadius( )
    {
        return referenceRadius_;
    }

    //! Function to get the cosine spherical harmonic coefficients (geodesy normalized)
    /*!
     *  Function to get the cosine spherical harmonic coefficients (geodesy normalized)
     *  \return Cosine spherical harmonic coefficients (geodesy normalized)
     */
    Eigen::MatrixXd getCosineCoefficients( )
    {
        return cosineCoefficients_;
    }

    //! Function to get the sine spherical harmonic coefficients (geodesy normalized)
    /*!
     *  Function to get the sine spherical harmonic coefficients (geodesy normalized)
     *  \return Sine spherical harmonic coefficients (geodesy normalized)
     */
    Eigen::MatrixXd getSineCoefficients( )
    {
        return sineCoefficients_;
    }

    //! Function to reset the cosine spherical harmonic coefficients (geodesy normalized)
    /*!
     *  Function to reset the cosine spherical harmonic coefficients (geodesy normalized)
     *  \param cosineCoefficients New cosine spherical harmonic coefficients (geodesy normalized)
     */
    void setCosineCoefficients( const Eigen::MatrixXd& cosineCoefficients )
    {
        cosineCoefficients_ = cosineCoefficients;

        if( !( updateInertiaTensor_ == nullptr ) )
        {
            updateInertiaTensor_( );
        }
    }

    //! Function to reset the cosine spherical harmonic coefficients (geodesy normalized)
    /*!
     *  Function to reset the cosine spherical harmonic coefficients (geodesy normalized)
     *  \param sineCoefficients New sine spherical harmonic coefficients (geodesy normalized)
     */
    void setSineCoefficients( const Eigen::MatrixXd& sineCoefficients )
    {
        sineCoefficients_ = sineCoefficients;
        if( !( updateInertiaTensor_ == nullptr ) )
        {
            updateInertiaTensor_( );
        }
    }

    //! Function to get a cosine spherical harmonic coefficient block (geodesy normalized)
    /*!
     *  Function to get a cosine spherical harmonic coefficient block (geodesy normalized)
    *   up to a given degree and order
     *  \param maximumDegree Maximum degree of coefficient block
     *  \param maximumOrder Maximum order of coefficient block
     *  \return Cosine spherical harmonic coefficients (geodesy normalized) up to given
     *  degree and order
     */
    Eigen::MatrixXd getCosineCoefficientsBlock( const int maximumDegree, const int maximumOrder )
    {
        return cosineCoefficients_.block( 0, 0, maximumDegree + 1, maximumOrder + 1 );
    }

    //! Function to get a sine spherical harmonic coefficient block (geodesy normalized)
    /*!
     *  Function to get a sine spherical harmonic coefficient block (geodesy normalized)
     *  up to a given degree and order
     *  \param maximumDegree Maximum degree of coefficient block
     *  \param maximumOrder Maximum order of coefficient block
     *  \return Sine spherical harmonic coefficients (geodesy normalized) up to given
     *  degree and order
     */
    Eigen::MatrixXd getSineCoefficientsBlock( const int maximumDegree, const int maximumOrder )
    {
        return sineCoefficients_.block( 0, 0, maximumDegree + 1, maximumOrder + 1 );
    }

    //! Get maximum degree of spherical harmonics gravity field expansion.
    /*!
     *  Returns the maximum degree of the spherical harmonics gravity field expansion.
     *  \return Degree of spherical harmonics expansion.
     */
    double getDegreeOfExpansion( )
    {
        return cosineCoefficients_.rows( ) + 1;
    }

    //! Get maximum order of spherical harmonics gravity field expansion.
    /*!
     *  Returns the maximum order of the spherical harmonics gravity field expansion.
     *  \return Order of spherical harmonics expansion.
     */
    double getOrderOfExpansion( )
    {
        return cosineCoefficients_.cols( ) + 1;
    }

    //! Function to calculate the gravitational potential at a given point
    /*!
     *  Function to calculate the gravitational potential due to this body at a given point.
     *  Note that this function, which has the same interface as in the base class and
     *  expands the gravity field to its maximum degree and order.
     *  \param bodyFixedPosition of point at which potential is to be calculated, in body-fixed
     *  frame.
     *  \return Gravitational potential at requested point.
     */
    double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        return getGravitationalPotential( bodyFixedPosition, cosineCoefficients_.rows( ) - 1,
                                          sineCoefficients_.cols( ) - 1 );
    }

    //! Function to calculate the gravitational potential due to terms up to given degree and
    //! order at a given point
    /*!
     *  Function to calculate the gravitational potential due to terms up to given degree and
     *  order due to this body at a given point.
     *  \param bodyFixedPosition Position of point at which potential is to be calculate,
     *  in body-fixed frame.
     *  \param maximumDegree Maximum degree of spherical harmonic coefficients to include.
     *  \param maximumOrder Maximum order of spherical harmonic coefficients to include.
     *  \param minimumDegree Minimum degree of spherical harmonic coefficients to include, default 0
     *  \param minimumOrder Maximum order of spherical harmonic coefficients to include, default 0
     *  \return Gravitational potential due to terms up to given degree and order at
     *  requested point.
     */
    double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition,
                                      const double maximumDegree,
                                      const double maximumOrder,
                                      const double minimumDegree = 0,
                                      const double minimumOrder = 0 )
    {
        return calculateSphericalHarmonicGravitationalPotential(
                    bodyFixedPosition, gravitationalParameter_, referenceRadius_,
                    cosineCoefficients_.block( 0, 0, maximumDegree + 1, maximumOrder + 1 ),
                    sineCoefficients_.block( 0, 0, maximumDegree + 1, maximumOrder + 1 ),
                    sphericalHarmonicsCache_,
                    minimumDegree, minimumOrder );
    }

    //! Get the gradient of the potential.
    /*!
     * Returns the gradient of the potential for the gravity field selected.
     *  Note that this function, which has the same interface as in the base class and
     *  expands the gravity field to its maximum degree and order.
     * \param bodyFixedPosition Position at which gradient of potential is to be determined
     * \return Gradient of potential.
     */
    Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        return getGradientOfPotential( bodyFixedPosition, cosineCoefficients_.rows( ),
                                       sineCoefficients_.cols( ) );
    }

    //! Get the gradient of the potential.
    /*!
     *  Returns the gradient of the potential for the gravity field selected.
     *  \param bodyFixedPosition Position at which gradient of potential is to be determined
     *  \param maximumDegree Maximum degree of spherical harmonic coefficients to include.
     *  \param maximumOrder Maximum order of spherical harmonic coefficients to include.
     *  \return Gradient of potential.
     */
    Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition,
                                            const double maximumDegree,
                                            const double maximumOrder )
    {
        std::map< std::pair< int, int >, Eigen::Vector3d > dummyMap;

        return computeGeodesyNormalizedGravitationalAccelerationSum(
                    bodyFixedPosition, gravitationalParameter_, referenceRadius_,
                    cosineCoefficients_.block( 0, 0, maximumDegree, maximumOrder ),
                    sineCoefficients_.block( 0, 0, maximumDegree, maximumOrder ), sphericalHarmonicsCache_, dummyMap );
    }

    //! Function to retrieve the tdentifier for body-fixed reference frame
    /*!
     *  Function to retrieve the tdentifier for body-fixed reference frame
     *  \return Function to retrieve the tdentifier for body-fixed reference frame.
     */
    std::string getFixedReferenceFrame( )
    {
        return fixedReferenceFrame_;
    }

    //! Function to retrieve if spherical harmonic coefficients are normalized
    /*!
     * Function to retrieve if spherical harmonic coefficients are normalized
     * \return Boolean stating whether spherical harmonic coefficients are normalized
     */
    bool areCoefficientsGeodesyNormalized( )
    {
        return true;
    }

    //! Function to get inertia tensor normalization factor
    /*!
     * Function to get inertia tensor normalization factor (M*R^2)
     * \return Inertia tensor normalization factor (M*R^2)
     */
    double getInertiaTensorNormalizationFactor( )
    {
        return gravitationalParameter_ * referenceRadius_ * referenceRadius_ / physical_constants::GRAVITATIONAL_CONSTANT;
    }

protected:

    //! Reference radius of spherical harmonic field expansion
    /*!
     *  Reference radius of spherical harmonic field expansion
     */
    double referenceRadius_;

    //! Cosine spherical harmonic coefficients (geodesy normalized)
    /*!
     *  Cosine spherical harmonic coefficients (geodesy normalized)
     */
    Eigen::MatrixXd cosineCoefficients_;

    //! Sine spherical harmonic coefficients (geodesy normalized)
    /*!
     *  Sine spherical harmonic coefficients (geodesy normalized)
     */
    Eigen::MatrixXd sineCoefficients_;

    //! Identifier for body-fixed reference frame
    /*!
     *  Identifier for body-fixed reference frame
     */
    std::string fixedReferenceFrame_;

    //! Cache object for potential calculations.
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;
};

//! Function to determine a body's inertia tensor from its degree two unnormalized gravity field coefficients
/*!
 * Function to determine a body's inertia tensor from its degree two unnormalized gravity field coefficients, and the mean
 * moment of inertia.
 * \param c20Coefficient Degree 2, order 0, unnormalized cosine spherical harmonic gravity field coefficient
 * \param c21Coefficient Degree 2, order 1, unnormalized cosine spherical harmonic gravity field coefficient
 * \param c22Coefficient Degree 2, order 2, unnormalized cosine spherical harmonic gravity field coefficient
 * \param s21Coefficient Degree 2, order 1, unnormalized sine spherical harmonic gravity field coefficient
 * \param s22Coefficient Degree 2, order 2, unnormalized sine spherical harmonic gravity field coefficient
 * \param scaledMeanMomentOfInertia  Mean moment of inertial, divided by (M*R^2)
 * \param bodyMass Mass M of body
 * \param referenceRadius Reference radius R of the spherical harmonic coefficients
 * \return Inertia tensor of body
 */
Eigen::Matrix3d getInertiaTensor(
        const double c20Coefficient,
        const double c21Coefficient,
        const double c22Coefficient,
        const double s21Coefficient,
        const double s22Coefficient,
        const double scaledMeanMomentOfInertia,
        const double bodyMass,
        const double referenceRadius );

//! Function to determine a body's inertia tensor from its unnormalized gravity field coefficients
/*!
 * Function to determine a body's inertia tensor from unnormalized gravity field coefficients, and the mean
 * moment of inertia.
 * \param unnormalizedCosineCoefficients Block of cosine spherical harmonic coefficients, degree/order stored in row/column
 * \param unnormalizedSineCoefficients Block of sine spherical harmonic coefficients, degree/order stored in row/column
 * \param scaledMeanMomentOfInertia  Mean moment of inertial, divided by (M*R^2)
 * \param bodyMass Mass M of body
 * \param referenceRadius Reference radius R of the spherical harmonic coefficients
 * \return Inertia tensor of body
 */
Eigen::Matrix3d getInertiaTensor(
        const Eigen::MatrixXd& unnormalizedCosineCoefficients,
        const Eigen::MatrixXd& unnormalizedSineCoefficients,
        const double scaledMeanMomentOfInertia,
        const double bodyMass,
        const double referenceRadius );

//! Function to determine a body's inertia tensor from its gravity field model
/*!
 * Function to determine a body's inertia tensor from its gravity field model, and the mean
 * moment of inertia.
 * \param sphericalHarmonicGravityField Spherical harmonic gravity field from which the inertia tensor is to be computed
 * \param scaledMeanMomentOfInertia  Mean moment of inertial, divided by (M*R^2), with M the mass of the body and R the
 * reference radius of the gravity field
 * \return Inertia tensor of body
 */
Eigen::Matrix3d getInertiaTensor(
        const std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicGravityField,
        const double scaledMeanMomentOfInertia );

//! Retrieve degree 2 spherical harmonic coefficients from inertia tensor and assiciated parameters
/*!
 * Retrieve degree 2 spherical harmonic coefficients from inertia tensor and assiciated parameters
 * \param inertiaTensor Inertia tensor
 * \param bodyGravitationalParameter Gravitational paramater of gravity field
 * \param referenceRadius Reference radius of gravity field
 * \param useNormalizedCoefficients Boolean stating whether spherical harmonic coefficients are normalized
 * \param cosineCoefficients Cosine coefficients of  of gravity field (returned by reference)
 * \param sineCoefficients Sine coefficients of  of gravity field (returned by reference)
 * \param scaledMeanMomentOfInertia Scaled mean moment of inertia (returned by reference)
 */
void getDegreeTwoSphericalHarmonicCoefficients(
        const Eigen::Matrix3d inertiaTensor, const double bodyGravitationalParameter, const double referenceRadius,
        const bool useNormalizedCoefficients,
        Eigen::MatrixXd& cosineCoefficients, Eigen::MatrixXd& sineCoefficients, double& scaledMeanMomentOfInertia );

std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > getDegreeTwoSphericalHarmonicCoefficients(
        const Eigen::Matrix3d inertiaTensor, const double bodyGravitationalParameter, const double referenceRadius,
        const int maximumCoefficientDegree = 2, const bool useNormalizedCoefficients = true );

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H
