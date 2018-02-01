/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Heiskanen, W.A., Moritz, H. Physical geodesy. Freeman, 1967.
 *
 *    Notes
 *      The class implementation currently only wraps the geodesy-normalized free function to
 *      compute the gravitational acceleration. Maybe in future, using an enum, the user can be
 *      given the choice of the free function to wrap.
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITY_MODEL_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITY_MODEL_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
/*!
 * This function computes the acceleration caused by gravitational spherical harmonics, with the
 * coefficients expressed using a geodesy-normalization. This acceleration is the summation of all
 * harmonic terms from degree and order zero, up to a user-specified highest degree and order. The
 * harmonic coefficients for the function must be provided in geodesy-normalized format. This
 * geodesy-normalization is defined as:
 * \f{eqnarray*}{
 *     \bar{ C }_{ n, m } = \Pi_{ n, m } C_{ n, m } \\
 *     \bar{ S }_{ n, m } = \Pi_{ n, m } S_{ n, m }
 * \f}
 * in which \f$ \bar{ C }_{ n, m } \f$ and \f$ \bar{ S }_{ n, m } \f$ are a geodesy-normalized
 * cosine and sine harmonic coefficient respectively (of degree \f$ n \f$ and order \f$ m \f$). The
 * unnormalized harmonic coefficients are represented by \f$ C_{ n, m } \f$ and \f$ S_{ n, m } \f$.
 * The normalization factor \f$ \Pi_{ n, m } \f$ is given by Heiskanen & Moritz [1967] as:
 * \f[
 *     \Pi_{ n, m } = \sqrt{ \frac{ ( n + m )! }{ ( 2 - \delta_{ 0, m } ) ( 2 n + 1 ) ( n - m )! } }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order and \f$ \delta_{ 0, m } \f$ is the
 * Kronecker delta.
 * \param positionOfBodySubjectToAcceleration Cartesian position vector with respect to the
 *          reference frame that is associated with the harmonic coefficients.
 *          The order is important!
 *          position( 0 ) = x coordinate [m],
 *          position( 1 ) = y coordinate [m],
 *          position( 2 ) = z coordinate [m].
 * \param gravitationalParameter Gravitational parameter associated with the spherical harmonics
 *          [m^3 s^-2].
 * \param equatorialRadius Reference radius of the spherical harmonics [m].
 * \param cosineHarmonicCoefficients Matrix with <B>geodesy-normalized</B> cosine harmonic
 *          coefficients. The row index indicates the degree and the column index indicates the order
 *          of coefficients.
 * \param sineHarmonicCoefficients Matrix with <B>geodesy-normalized</B> sine harmonic coefficients.
 *          The row index indicates the degree and the column index indicates the order of
 *          coefficients. The matrix must be equal in size to cosineHarmonicCoefficients.
 * \param sphericalHarmonicsCache Cache object for computing/retrieving repeated terms in spherical harmonics potential
 *          gradient calculation.
 * \return Cartesian acceleration vector resulting from the summation of all harmonic terms.
 *           The order is important!
 *           acceleration( 0 ) = x acceleration [m s^-2],
 *           acceleration( 1 ) = y acceleration [m s^-2],
 *           acceleration( 2 ) = z acceleration [m s^-2].
 */
Eigen::Vector3d computeGeodesyNormalizedGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameter,
        const double equatorialRadius,
        const Eigen::MatrixXd& cosineHarmonicCoefficients,
        const Eigen::MatrixXd& sineHarmonicCoefficients,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Compute gravitational acceleration due to single spherical harmonics term.
/*!
 * This function computes the acceleration caused by a single gravitational spherical harmonics
 * term, with the coefficients expressed using a geodesy-normalization. The harmonic coefficients
 * for the function must be provided in geodesy-normalized format. This geodesy-normalization is
 * defined as:
 * \f{eqnarray*}{
 *     \bar{ C }_{ n, m } = \Pi_{ n, m } C_{ n, m } \\
 *     \bar{ S }_{ n, m } = \Pi_{ n, m } S_{ n, m }
 * \f}
 * in which \f$ \bar{ C }_{ n, m } \f$ and \f$ \bar{ S }_{ n, m } \f$ are a geodesy-normalized
 * cosine and sine harmonic coefficient respectively (of degree \f$ n \f$ and order \f$ m \f$). The
 * unnormalized harmonic coefficients are represented by \f$ C_{ n, m } \f$ and \f$ S_{ n, m } \f$.
 * The normalization factor \f$ \Pi_{ n, m } \f$ is given by Heiskanen & Moritz [1967] as:
 * \f[
 *     \Pi_{ n, m } = \sqrt{ \frac{ ( n + m )! }{ ( 2 - \delta_{ 0, m } ) ( 2 n + 1 ) ( n - m )! } }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order and \f$ \delta_{ 0, m } \f$ is the
 * Kronecker delta.
 * \param positionOfBodySubjectToAcceleration Cartesian position vector with respect to the
 *          reference frame that is associated with the harmonic coefficients.
 *          The order is important!
 *          position( 0 ) = x coordinate [m],
 *          position( 1 ) = y coordinate [m],
 *          position( 2 ) = z coordinate [m].
 * \param degree Degree of the harmonic term.
 * \param order Order of the harmonic term.
 *  * \param cosineHarmonicCoefficient <B>Geodesy-normalized</B> cosine harmonic
 *          coefficient.
 * \param sineHarmonicCoefficient <B>Geodesy-normalized</B> sine harmonic coefficient.
 * \param gravitationalParameter Gravitational parameter associated with the spherical harmonic
 *          [m^3 s^-2].
 * \param equatorialRadius Reference radius of the spherical harmonic [m].
 * \param sphericalHarmonicsCache Cache object for computing/retrieving repeated terms in spherical harmonics potential
 *          gradient calculation.
 * \return Cartesian acceleration vector resulting from the spherical harmonic term.
 *           The order is important!
 *           acceleration( 0 ) = x acceleration [m s^-2],
 *           acceleration( 1 ) = y acceleration [m s^-2],
 *           acceleration( 2 ) = z acceleration [m s^-2].
 */
Eigen::Vector3d computeSingleGeodesyNormalizedGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameter,
        const double equatorialRadius,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Template class for general spherical harmonics gravitational acceleration model.
/*!
 * This templated class implements a general spherical harmonics gravitational acceleration model.
 * The acceleration computed with this class is based on the geodesy-normalization described by
 * (Heiskanen & Moritz, 1967), implemented in the
 * computeGeodesyNormalizedGravitationalAccelerationSum() function. The acceleration computed is a
 * sum, based on the matrix of coefficients of the model provided.
 */
class SphericalHarmonicsGravitationalAccelerationModel
        : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >,
        public SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d >
{
private:

    //! Typedef for base class.
    typedef SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > Base;

    //! Typedef for coefficient-matrix-returning function.
    typedef boost::function< Eigen::MatrixXd( ) > CoefficientMatrixReturningFunction;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of spherical
    //! harmonics expansion.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter and equatorial radius of the
     * body exerting the acceleration, constant coefficient matrices for the spherical harmonics
     * expansion, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * Boost::lambda library to create a function on-the-fly that returns the constant
     * gravitational parameter, equatorial radius and coefficient matrices provided. The
     * constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param anEquatorialRadius A (constant) equatorial radius [m].
     * \param aCosineHarmonicCoefficientMatrix A (constant) cosine harmonic coefficient matrix.
     * \param aSineHarmonicCoefficientMatrix A (constant) sine harmonic coefficient matrix.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param sphericalHarmonicsCache Cache object for computing/retrieving repeated terms in spherical harmonics potential
     *          gradient calculation.
     */
    SphericalHarmonicsGravitationalAccelerationModel(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const double anEquatorialRadius,
            const Eigen::MatrixXd aCosineHarmonicCoefficientMatrix,
            const Eigen::MatrixXd aSineHarmonicCoefficientMatrix,
            const StateFunction positionOfBodyExertingAccelerationFunction
            = boost::lambda::constant( Eigen::Vector3d::Zero( ) ),
            const boost::function< Eigen::Quaterniond( ) >
            rotationFromBodyFixedToIntegrationFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
            const bool isMutualAttractionUsed = 0,
            boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
            boost::make_shared< basic_mathematics::SphericalHarmonicsCache >( ) )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameter,
                positionOfBodyExertingAccelerationFunction,
                isMutualAttractionUsed ),
          equatorialRadius( anEquatorialRadius ),
          getCosineHarmonicsCoefficients(
              boost::lambda::constant(aCosineHarmonicCoefficientMatrix ) ),
          getSineHarmonicsCoefficients( boost::lambda::constant(aSineHarmonicCoefficientMatrix ) ),
          rotationFromBodyFixedToIntegrationFrameFunction_(
              rotationFromBodyFixedToIntegrationFrameFunction ),
          sphericalHarmonicsCache_( sphericalHarmonicsCache ),
          currentAcceleration_( Eigen::Vector3d::Zero( ) )

    {
        sphericalHarmonicsCache_->resetMaximumDegreeAndOrder(
                    std::max< int >( static_cast< int >( getCosineHarmonicsCoefficients( ).rows( ) ), sphericalHarmonicsCache_->getMaximumDegree( ) ),
                    std::max< int >( static_cast< int >( getCosineHarmonicsCoefficients( ).cols( ) ), sphericalHarmonicsCache_->getMaximumOrder( ) ) + 1 );
        this->updateMembers( );
    }

    //! Constructor taking functions for position of bodies, and parameters of spherical harmonics
    //! expansion.
    /*!
     * Constructor taking pointer to functions returning the position of the body subject to
     * gravitational acceleration, the gravitational parameter of the body exerting the
     * acceleration (central body), the equatorial radius of the central body, the coefficient
     * matrices of the spherical harmonics expansion, and the position of the central body. The
     * constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameterFunction Pointer to function returning gravitational parameter.
     * \param anEquatorialRadius Pointer to function returning equatorial radius.
     * \param cosineHarmonicCoefficientsFunction Pointer to function returning matrix of
                cosine-coefficients of spherical harmonics expansion.
     * \param sineHarmonicCoefficientsFunction Pointer to function returning matrix of
                sine-coefficients of spherical harmonics expansion.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param sphericalHarmonicsCache Cache object for computing/retrieving repeated terms in spherical harmonics potential
     */
    SphericalHarmonicsGravitationalAccelerationModel(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const boost::function< double( ) > aGravitationalParameterFunction,
            const double anEquatorialRadius,
            const CoefficientMatrixReturningFunction cosineHarmonicCoefficientsFunction,
            const CoefficientMatrixReturningFunction sineHarmonicCoefficientsFunction,
            const StateFunction positionOfBodyExertingAccelerationFunction
            = boost::lambda::constant( Eigen::Vector3d::Zero( ) ),
            const boost::function< Eigen::Quaterniond( ) >
            rotationFromBodyFixedToIntegrationFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
            const bool isMutualAttractionUsed = 0,
            boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache
            = boost::make_shared< basic_mathematics::SphericalHarmonicsCache >( ) )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameterFunction,
                positionOfBodyExertingAccelerationFunction,
                isMutualAttractionUsed ),
          equatorialRadius( anEquatorialRadius ),
          getCosineHarmonicsCoefficients( cosineHarmonicCoefficientsFunction ),
          getSineHarmonicsCoefficients( sineHarmonicCoefficientsFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
          sphericalHarmonicsCache_( sphericalHarmonicsCache ),
          currentAcceleration_( Eigen::Vector3d::Zero( ) )
    {
        sphericalHarmonicsCache_->resetMaximumDegreeAndOrder(
                    std::max< int >( static_cast< int >( getCosineHarmonicsCoefficients( ).rows( ) ), sphericalHarmonicsCache_->getMaximumDegree( ) ),
                    std::max< int >( static_cast< int >( getCosineHarmonicsCoefficients( ).cols( ) ), sphericalHarmonicsCache_->getMaximumOrder( ) ) + 1 );


        this->updateMembers( );
    }

    //! Get gravitational acceleration.
    /*!
     * Returns the gravitational acceleration computed using the input parameters provided to the
     * class. This function serves as a wrapper for the
     * computeGeodesyNormalizedGravitationalAccelerationSum() function.
     * \return Computed gravitational acceleration vector.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class
     * members of this class.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            cosineHarmonicCoefficients = getCosineHarmonicsCoefficients( );
            sineHarmonicCoefficients = getSineHarmonicsCoefficients( );
            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
            this->updateBaseMembers( );
            currentAcceleration_ = rotationToIntegrationFrame_ *
                    computeGeodesyNormalizedGravitationalAccelerationSum(
                        rotationToIntegrationFrame_.inverse( ) * (
                            this->positionOfBodySubjectToAcceleration - this->positionOfBodyExertingAcceleration ),
                        gravitationalParameter,
                        equatorialRadius,
                        cosineHarmonicCoefficients,
                        sineHarmonicCoefficients, sphericalHarmonicsCache_ );
        }
    }

    //! Function to retrieve the spherical harmonics cache for this acceleration.
    /*!
     *  Function to retrieve the spherical harmonics cache for this acceleration.
     *  \return Spherical harmonics cache for this acceleration
     */
    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > getSphericalHarmonicsCache( )
    {
        return sphericalHarmonicsCache_;
    }

    //! Function to retrieve the spherical harmonics reference radius.
    /*!
     *  Function to retrieve the spherical harmonics reference radius.
     *  \return Spherical harmonics reference radius.
     */
    double getReferenceRadius( )
    {
        return equatorialRadius;
    }

    //! Matrix of cosine coefficients.
    /*!
     * Matrix containing coefficients of cosine terms for spherical harmonics expansion.
     */
    CoefficientMatrixReturningFunction getCosineHarmonicCoefficientsFunction( )
    {
        return getCosineHarmonicsCoefficients;
    }

    //! Matrix of sine coefficients.
    /*!
     * Matrix containing coefficients of sine terms for spherical harmonics expansion.
     */
    CoefficientMatrixReturningFunction getSineHarmonicCoefficientsFunction( )
    {
        return getSineHarmonicsCoefficients;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
    /*!
     *  Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
     *  \return current rotation from body-fixed frame to integration frame, in the form of a quaternion.
     */
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
    /*!
     *  Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
     *  \return current rotation from body-fixed frame to integration frame, as a rotation matrix.
     */
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

protected:

private:

    //! Equatorial radius [m].
    /*!
     * Current value of equatorial (planetary) radius used for spherical harmonics expansion [m].
    */
    const double equatorialRadius;

    //! Matrix of cosine coefficients.
    /*!
     * Matrix containing coefficients of cosine terms for spherical harmonics expansion.
     */
    Eigen::MatrixXd cosineHarmonicCoefficients;

    //! Matrix of sine coefficients.
    /*!
     * Matrix containing coefficients of sine terms for spherical harmonics expansion.
     */
    Eigen::MatrixXd sineHarmonicCoefficients;

    //! Pointer to function returning cosine harmonics coefficients matrix.
    /*!
     * Pointer to function that returns the current coefficients of the cosine terms of the
     * spherical harmonics expansion.
     */
    const CoefficientMatrixReturningFunction getCosineHarmonicsCoefficients;

    //! Pointer to function returning sine harmonics coefficients matrix.
    /*!
     * Pointer to function that returns the current coefficients of the sine terms of the
     * spherical harmonics expansion.
     */
    const CoefficientMatrixReturningFunction getSineHarmonicsCoefficients;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    boost::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //!  Spherical harmonics cache for this acceleration
    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    //! Current acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAcceleration_;

};


//! Typedef for shared-pointer to SphericalHarmonicsGravitationalAccelerationModel.
typedef boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel >
SphericalHarmonicsGravitationalAccelerationModelPointer;


} // namespace gravitation

} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITY_MODEL_H
