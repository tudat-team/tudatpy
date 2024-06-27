/*
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

#include <iostream>

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <memory>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModelBase.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
{
namespace gravitation
{

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
    typedef std::function< Eigen::MatrixXd( ) > CoefficientMatrixReturningFunction;

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
            const StateFunction positionOfBodyExertingAccelerationFunction =
            [ ]( Eigen::Vector3d& input ){ input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) >
            rotationFromBodyFixedToIntegrationFrameFunction =
            [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
            std::make_shared< basic_mathematics::SphericalHarmonicsCache >( ) )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameter,
                positionOfBodyExertingAccelerationFunction,
                isMutualAttractionUsed ),
          equatorialRadius( anEquatorialRadius ),
          getCosineHarmonicsCoefficients( [ = ]( ){ return aCosineHarmonicCoefficientMatrix; } ),
          getSineHarmonicsCoefficients( [ = ]( ){ return aSineHarmonicCoefficientMatrix; } ),
          rotationFromBodyFixedToIntegrationFrameFunction_(
              rotationFromBodyFixedToIntegrationFrameFunction ),
          sphericalHarmonicsCache_( sphericalHarmonicsCache ),
          saveSphericalHarmonicTermsSeparately_( false )
    {
        maximumDegree_ = static_cast< int >( getCosineHarmonicsCoefficients( ).rows( ) ) - 1 ;
        maximumOrder_ = static_cast< int >( getCosineHarmonicsCoefficients( ).cols( ) )- 1 ;
        sphericalHarmonicsCache_->resetMaximumDegreeAndOrder(
                    std::max< int >( maximumDegree_,
                                     sphericalHarmonicsCache_->getMaximumDegree( ) ) + 1,
                    std::max< int >( maximumOrder_,
                                     sphericalHarmonicsCache_->getMaximumOrder( ) ) + 1 );
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
            const std::function< double( ) > aGravitationalParameterFunction,
            const double anEquatorialRadius,
            const CoefficientMatrixReturningFunction cosineHarmonicCoefficientsFunction,
            const CoefficientMatrixReturningFunction sineHarmonicCoefficientsFunction,
            const StateFunction positionOfBodyExertingAccelerationFunction =
            [ ]( Eigen::Vector3d& input ){ input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) >
            rotationFromBodyFixedToIntegrationFrameFunction =
            [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache
            = std::make_shared< basic_mathematics::SphericalHarmonicsCache >( ) )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameterFunction,
                positionOfBodyExertingAccelerationFunction,
                isMutualAttractionUsed ),
          equatorialRadius( anEquatorialRadius ),
          getCosineHarmonicsCoefficients( cosineHarmonicCoefficientsFunction ),
          getSineHarmonicsCoefficients( sineHarmonicCoefficientsFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
          sphericalHarmonicsCache_( sphericalHarmonicsCache ),
          saveSphericalHarmonicTermsSeparately_( false )
    {
        maximumDegree_ = static_cast< int >( getCosineHarmonicsCoefficients( ).rows( ) - 1 );
        maximumOrder_ = static_cast< int >( getCosineHarmonicsCoefficients( ).cols( ) - 1 );
        sphericalHarmonicsCache_->resetMaximumDegreeAndOrder(
                    std::max< int >( maximumDegree_ + 1,
                                     sphericalHarmonicsCache_->getMaximumDegree( ) ),
                    std::max< int >( maximumOrder_ + 1,
                                     sphericalHarmonicsCache_->getMaximumOrder( ) ) + 1 );


    }

    //! Get gravitational acceleration in body-fixed frame of body undergoing acceleration.
    /*!
     * Returns the gravitational acceleration in body-fixed frame of body undergoing acceleration computed
     * computed by the updateMembers function.
     * \return Computed gravitational acceleration vector.
     */
    Eigen::Vector3d getAccelerationInBodyFixedFrame( )
    {
        return currentAccelerationInBodyFixedFrame_;
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

            currentInertialRelativePosition_ =
                    this->positionOfBodySubjectToAcceleration - this->positionOfBodyExertingAcceleration ;

            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * (
                        currentInertialRelativePosition_ );

            currentAcceleration_ =
                    computeGeodesyNormalizedGravitationalAccelerationSum(
                        currentRelativePosition_,
                        gravitationalParameter,
                        equatorialRadius,
                        cosineHarmonicCoefficients,
                        sineHarmonicCoefficients, sphericalHarmonicsCache_,
                        accelerationPerTerm_,
                        saveSphericalHarmonicTermsSeparately_,
                        rotationToIntegrationFrame_.toRotationMatrix( ) );
            currentAccelerationInBodyFixedFrame_ = rotationToIntegrationFrame_.inverse( ) * currentAcceleration_;

            if ( this->updatePotential_ )
            {
                this->currentPotential_ = gravitation::calculateSphericalHarmonicGravitationalPotential(
                        currentRelativePosition_,
                        gravitationalParameter,
                        equatorialRadius,
                        cosineHarmonicCoefficients,
                        sineHarmonicCoefficients,
                        sphericalHarmonicsCache_ );
            }
        }
    }

    //! Function to retrieve total spherical harmonic acceleration in inertial frame, with alternative coefficients
    /*!
     * Function to retrieve total spherical harmonic acceleration in inertial frame, with alternative coefficients, e.g.
     * different from those defined in this class
     * \param cosineCoefficients Cosine coefficients to use
     * \param sineCoefficients Sine coefficients to use
     * \return Total spherical harmonic acceleration in inertial frame, with alternative coefficients
     */
    Eigen::VectorXd getAccelerationWithAlternativeCoefficients(
            const Eigen::MatrixXd& cosineCoefficients, const Eigen::MatrixXd& sineCoefficients)
    {
        std::map< std::pair< int, int >, Eigen::Vector3d > dummy;
        return computeGeodesyNormalizedGravitationalAccelerationSum(
                    currentRelativePosition_,
                    gravitationalParameter,
                    equatorialRadius,
                    cosineCoefficients,
                    sineCoefficients, sphericalHarmonicsCache_,
                    dummy,
                    false,
                    rotationToIntegrationFrame_.toRotationMatrix( ) );
    }

    //! Function to retrieve spherical harmonic acceleration in inertial frame, with alternative coefficients, per term
    /*!
     * Function to retrieve spherical harmonic acceleration in inertial frame, with alternative coefficients, per term.
     * The acceleration is calculated separately for the contribution of each degree/order of the spherical harmonic field
     * \param cosineCoefficients Cosine coefficients to use
     * \param sineCoefficients Sine coefficients to use
     * \param coefficientIndices Degrees and orders for which acceleration contributions are to be determined.
     * \return Total spherical harmonic acceleration in inertial frame, with alternative coefficients
     */
    Eigen::VectorXd getAccelerationComponentsWithAlternativeCoefficients(
            const Eigen::MatrixXd& cosineCoefficients, const Eigen::MatrixXd& sineCoefficients,
            const std::vector< std::pair< int, int > >& coefficientIndices )
    {
        std::map< std::pair< int, int >, Eigen::Vector3d > accelerationPerTerm;

        computeGeodesyNormalizedGravitationalAccelerationSum(
                    currentRelativePosition_,
                    gravitationalParameter,
                    equatorialRadius,
                    cosineCoefficients,
                    sineCoefficients, sphericalHarmonicsCache_,
                    accelerationPerTerm,
                    true,
                    rotationToIntegrationFrame_.toRotationMatrix( ) );


        Eigen::VectorXd returnVector = Eigen::VectorXd( 3 * coefficientIndices.size( ) );
        for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
        {
            returnVector.segment( i * 3, 3 ) = accelerationPerTerm.at( coefficientIndices.at( i ) );
        }
        return returnVector;
   }

    //! Function to retrieve the spherical harmonics cache for this acceleration.
    /*!
     *  Function to retrieve the spherical harmonics cache for this acceleration.
     *  \return Spherical harmonics cache for this acceleration
     */
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > getSphericalHarmonicsCache( )
    {
        return sphericalHarmonicsCache_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
    //! fixed to body undergoing acceleration
    /*!
     * Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
     * fixed to bodyundergoing acceleration
     * \return Current position vector from body exerting acceleration to body undergoing acceleration, in frame
     * fixed to bodyundergoing acceleration
     */
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
    //! frame
    /*!
     * Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
     * frame
     * \return Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
     */
    Eigen::Vector3d getCurrentInertialRelativePosition( )
    {
        return currentInertialRelativePosition_;
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

    //! Function to set whether each of the separate spherical harmonic terms should be saved
    /*!
     * Function to set whether each of the separate spherical harmonic terms should be saved (in accelerationPerTerm_ member
     * variable of this class)
     * \param saveSphericalHarmonicTermsSeparately Boolean denoting whether each of the separate spherical harmonic terms should
     * be saved (in accelerationPerTerm_ member variable of this class
     */
    void setSaveSphericalHarmonicTermsSeparately( const bool saveSphericalHarmonicTermsSeparately )
    {
        saveSphericalHarmonicTermsSeparately_ = saveSphericalHarmonicTermsSeparately;
    }

    //! Function to retrieve the contributions of separate degrees/ordesr to the acceleration, concatenated in a single vector
    /*!
     * Function to retrieve the contributions of specific separate degree/order to the acceleration, concatenated in a single
     * vector
     * \param coefficientIndices List of degree/order at which the contributions to the full acceleration are to be retrieved
     * \return Contributions of separate degrees/ordesr to the acceleration, concatenated in a single vector
     */
    Eigen::VectorXd getConcatenatedAccelerationComponents( const std::vector< std::pair< int, int > >& coefficientIndices )
    {
        if( !saveSphericalHarmonicTermsSeparately_ )
        {
            throw std::runtime_error( "Error when retrieving component accelerations from spherial harmonic acceleration, components not saved" );
        }

        Eigen::VectorXd returnVector = Eigen::VectorXd( 3 * coefficientIndices.size( ) );
        for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
        {
            if( accelerationPerTerm_.count( coefficientIndices.at( i ) ) != 0 )
            {
                returnVector.segment( i * 3, 3 ) = accelerationPerTerm_.at( coefficientIndices.at( i ) );
            }
            else
            {
                throw std::runtime_error( "Error when retrieving spherical harmonic acceleration at degree/order: " +
                                          std::to_string( coefficientIndices.at( i ).first ) + "/" +
                                          std::to_string( coefficientIndices.at( i ).second ) +
                                          ". This degree/order combination is not within the selected range of the current acceleration model." );
            }

        }
        return returnVector;
    }

    Eigen::VectorXd getConcatenatedAccelerationComponentNorms( const std::vector< std::pair< int, int > >& coefficientIndices )
    {
        if( !saveSphericalHarmonicTermsSeparately_ )
        {
            throw std::runtime_error( "Error when retrieving component accelerations from spherial harmonic acceleration, components not saved" );
        }

        Eigen::VectorXd returnVector = Eigen::VectorXd( coefficientIndices.size( ) );
        for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
        {
            if( accelerationPerTerm_.count( coefficientIndices.at( i ) ) != 0 )
            {
                returnVector( i ) = accelerationPerTerm_.at( coefficientIndices.at( i ) ).norm( );
            }
            else
            {
                throw std::runtime_error( "Error when retrieving spherical harmonic acceleration at degree/order: " +
                                          std::to_string( coefficientIndices.at( i ).first ) + "/" +
                                          std::to_string( coefficientIndices.at( i ).second ) +
                                          ". This degree/order combination is not within the selected range of the current acceleration model." );
            }
        }
        return returnVector;
    }

    //! Function to retrieve maximum degree of gravity field expansion
    /*!
     * Function to retrieve maximum degree of gravity field expansion
     * \return Maximum degree of gravity field expansion
     */
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    //! Function to retrieve maximum order of gravity field expansion
    /*!
     * Function to retrieve maximum order of gravity field expansion
     * \return Maximum order of gravity field expansion
     */
    int getMaximumOrder( )
    {
        return maximumOrder_;
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
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    Eigen::Vector3d currentRelativePosition_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
    Eigen::Vector3d currentInertialRelativePosition_;

    //!  Spherical harmonics cache for this acceleration
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    //! Current acceleration in frame fixed to body undergoing acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

    //! List of contributions to accelerations at given degrees/orders, represented by first/second entry of map key pair.
    std::map< std::pair< int, int >, Eigen::Vector3d > accelerationPerTerm_;

    //! Boolean that denotes whether each of the separate spherical harmonic terms should be saved (in accelerationPerTerm_)
    bool saveSphericalHarmonicTermsSeparately_;

    //! Maximum degree of gravity field expansion
    int maximumDegree_;

    //! Maximum order of gravity field expansion
    int maximumOrder_;

};


//! Typedef for shared-pointer to SphericalHarmonicsGravitationalAccelerationModel.
typedef std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel >
SphericalHarmonicsGravitationalAccelerationModelPointer;


} // namespace gravitation

} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITY_MODEL_H
