#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
#define MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H

#include <tuple>
#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/legendrePolynomials.h"

#include "tudat/astro/gravitation/mutualForcePotential.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/sphericalHarmonicTransformations.h"

namespace tudat
{

namespace gravitation
{

//! Full two-body mutual spherical-harmonic acceleration model.
/*!
 * Computes the translational acceleration from the mutual potential using the effective one-body formulation
 * of Dirkx et al. (2019), in particular the effective coefficient mapping (Eqs. (47)-(48)) and potential
 * summation (Eq. (49)), evaluated as a Cartesian gradient for the translational equations (Eq. (55)).
 */
class FullTwoBodySphericalHarmonicAcceleration : public basic_astrodynamics::AccelerationModel3d
{
public:
    //! Constructor for the full two-body mutual spherical-harmonic acceleration model.
    FullTwoBodySphericalHarmonicAcceleration(
            const std::function< Eigen::Vector3d( ) > positionOfBody1Function,
            const std::function< Eigen::Vector3d( ) > positionOfBody2Function,
            const std::function< double( ) > gravitationalParameterFunction,
            const double equatorialRadiusOfBody1,
            const double equatorialRadiusOfBody2,
            const std::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody1Function,
            const std::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody1Function,
            const std::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody2Function,
            const std::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody2Function,
            const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
            const std::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation,
            const std::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation,
            const bool isMutualAttractionUsed,
            const bool areCoefficientsNormalized = 1 );

    ~FullTwoBodySphericalHarmonicAcceleration( ) {}

    //! Update current acceleration and intermediate states for the provided epoch.
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime == currentTime_ ) )
        {
            currentRotationFromInertialToBody1_ = toLocalFrameOfBody1Transformation_( );
            currentRotationFromBody2ToBody1_ = currentRotationFromInertialToBody1_ * toLocalFrameOfBody2Transformation_( ).inverse( );

            const Eigen::Vector3d positionBody1 = positionOfBody1Function_( );
            const Eigen::Vector3d positionBody2 = positionOfBody2Function_( );
            currentRelativePosition_ = positionBody1 - positionBody2;
            currentBodyFixedRelativePosition_ = currentRotationFromInertialToBody1_ * currentRelativePosition_;
            const double currentGravitationalParameter = gravitationalParameterFunction_( );

            effectiveMutualPotentialField_->computeCurrentEffectiveCoefficients( currentRotationFromBody2ToBody1_ );

            const double currentDistance = currentRelativePosition_.norm( );
            for( int i = 0; i <= effectiveMutualPotentialField_->getMaximumDegree1( ); i++ )
            {
                radius1Powers_[ i ] = basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody1_ / currentDistance, i );
            }

            for( int i = 0; i <= effectiveMutualPotentialField_->getMaximumDegree2( ); i++ )
            {
                radius2Powers_[ i ] = basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody2_ / currentDistance, i );
            }

            if( areCoefficientsNormalized_ )
            {
                mutualPotentialGradient_ =
                        computeGeodesyNormalizedMutualGravitationalAccelerationSum( currentBodyFixedRelativePosition_,
                                                                                    currentGravitationalParameter,
                                                                                    equatorialRadiusOfBody1_,
                                                                                    equatorialRadiusOfBody2_,
                                                                                    effectiveCosineCoefficientFunction_,
                                                                                    effectiveSineCoefficientFunction_,
                                                                                    coefficientCombinationsToUse_,
                                                                                    effectiveMutualPotentialField_->getMaximumDegree1( ),
                                                                                    effectiveMutualPotentialField_->getMaximumDegree2( ),
                                                                                    maximumDegree_,
                                                                                    radius1Powers_,
                                                                                    radius2Powers_,
                                                                                    sphericalHarmonicsCache_ );
            }
            else
            {
                mutualPotentialGradient_ = computeUnnormalizedMutualGravitationalAccelerationSum( currentBodyFixedRelativePosition_,
                                                                                                  currentGravitationalParameter,
                                                                                                  equatorialRadiusOfBody1_,
                                                                                                  equatorialRadiusOfBody2_,
                                                                                                  effectiveCosineCoefficientFunction_,
                                                                                                  effectiveSineCoefficientFunction_,
                                                                                                  coefficientCombinationsToUse_,
                                                                                                  sphericalHarmonicsCache_ );
            }

            currentAcceleration_ = currentRotationFromInertialToBody1_.inverse( ) * mutualPotentialGradient_;
            basic_astrodynamics::AccelerationModel3d::currentAcceleration_ = currentAcceleration_;
            currentTime_ = currentTime;
        }
    }

    //! Return the most recently computed acceleration.
    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    //! Return whether mutual attraction is included in the gravitational parameter.
    /*!
     *  If true, the effective gravitational parameter is the sum of both bodies' gravitational parameters.
     *  If false, only the accelerating body's parameter is used.
     */
    bool getIsMutualAttractionUsed( )
    {
        return isMutualAttractionUsed_;
    }

    //! Return the current relative position (body 1 minus body 2) in inertial coordinates.
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Return the current relative position expressed in the body-1-fixed frame.
    Eigen::Vector3d getCurrentBodyFixedRelativePosition( )
    {
        return currentBodyFixedRelativePosition_;
    }

    //! Return the current mutual potential gradient in body-1-fixed Cartesian coordinates.
    Eigen::Vector3d getMutualPotentialGradient( )
    {
        return mutualPotentialGradient_;
    }

    //! Return the current inertial-to-body-1 quaternion.
    Eigen::Quaterniond getCurrentRotationFromInertialToBody1( )
    {
        return currentRotationFromInertialToBody1_;
    }

    //! Return the current body-2-to-body-1 quaternion.
    Eigen::Quaterniond getCurrentRotationFromBody2ToBody1( )
    {
        return currentRotationFromBody2ToBody1_;
    }

    //! Return the effective mutual potential field object used by this model.
    std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > getEffectiveMutualPotentialField( )
    {
        return effectiveMutualPotentialField_;
    }

    //! Return the spherical harmonics cache used during acceleration evaluation.
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > getSphericalHarmonicsCache( )
    {
        return sphericalHarmonicsCache_;
    }

    //! Return the reference radius of body 1.
    double getEquatorialRadiusOfBody1( )
    {
        return equatorialRadiusOfBody1_;
    }

    //! Return the reference radius of body 2.
    double getEquatorialRadiusOfBody2( )
    {
        return equatorialRadiusOfBody2_;
    }

    //! Return whether input/output coefficients are geodesy normalized.
    bool getAreCoefficientsNormalized( )
    {
        return areCoefficientsNormalized_;
    }

    //! Return the current effective gravitational parameter.
    double getCurrentGravitationalParameter( )
    {
        return gravitationalParameterFunction_( );
    }

    //! Return precomputed power (R1/r)^index for the current state.
    double getRadius1Power( const int index )
    {
        return radius1Powers_.at( index );
    }

    //! Return precomputed power (R2/r)^index for the current state.
    double getRadius2Power( const int index )
    {
        return radius2Powers_.at( index );
    }

    //! Return full list of precomputed powers (R1/r)^l.
    const std::vector< double >& getRadius1Powers( ) const
    {
        return radius1Powers_;
    }

    //! Return full list of precomputed powers (R2/r)^l.
    const std::vector< double >& getRadius2Powers( ) const
    {
        return radius2Powers_;
    }

private:
    //! Function returning the position of body 1
    std::function< Eigen::Vector3d( ) > positionOfBody1Function_;

    //! Function returning the position of body 2
    std::function< Eigen::Vector3d( ) > positionOfBody2Function_;

    //! Function returning the effective gravitational parameter
    std::function< double( ) > gravitationalParameterFunction_;

    //! Equatorial radius of body 1.
    double equatorialRadiusOfBody1_;

    //! Equatorial radius of body 2.
    double equatorialRadiusOfBody2_;

    //! List of degrees/orders that are to be used for the series expansion.
    /*!
     * List of degrees/orders that are to be used for the series expansion. Each tuple contains: (degree of body 1, order of
     * body 1, degree of body 2, order of body 2)
     */
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    //! Function that returns the effective one-body spherical harmonic cosine coefficient
    /*!
     *  Function that returns the effective one-body spherical harmonic cosine coefficient as a function of (degree of body 1,
     *  order of body 1, degree of body 2, order of body 2)
     */
    std::function< double( int, int, int, int ) > effectiveCosineCoefficientFunction_;

    //! Function that returns the effective one-body spherical harmonic sine coefficient
    /*!
     *  Function that returns the effective one-body spherical harmonic sine coefficient as a function of (degree of body 1,
     *  order of body 1, degree of body 2, order of body 2)
     */
    std::function< double( int, int, int, int ) > effectiveSineCoefficientFunction_;

    //! Function that returns the rotation from the inertial to body-fixed frame of body 1
    std::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation_;

    //! Function that returns the rotation from the inertial to body-fixed frame of body 2
    std::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation_;

    //! Maximum degree of effective one-body spherical harmonic expansion
    unsigned int maximumDegree_;

    //! Maximum order of effective one-body spherical harmonic expansion
    unsigned int maximumOrder_;

    //! Object that computes the effect one-body coefficients from the two one-body gravoity fields and their relative orientation
    std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;

    //! Function that computes the coefficients of body, in the body-fixed frame of body 1
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    //! Acceleration, as computed by last call to updateMembers
    Eigen::Vector3d currentAcceleration_;

    //! Position of body 1 w.r.t. body 2, expressed in inertial frame, as computed by last call to updateMembers
    Eigen::Vector3d currentRelativePosition_;

    //! Position of body 1 w.r.t. body 2, expressed in body-fixed frame of body 1, as computed by last call to updateMembers
    Eigen::Vector3d currentBodyFixedRelativePosition_;

    //! Current gradient of mutual potential, in body-fixed frame of body 1 in Cartesian coordinates
    Eigen::Vector3d mutualPotentialGradient_;

    //! Current rotation from body-fixed frame of body 2 to that of body 1, as computed by last call to updateMembers
    Eigen::Quaterniond currentRotationFromBody2ToBody1_;

    //! Current rotation from inertial frame to body-fixed frame of body 1, as computed by last call to updateMembers
    Eigen::Quaterniond currentRotationFromInertialToBody1_;

    //! Boolean denoting whether mutual attraction is used in the effective gravitational parameter.
    bool isMutualAttractionUsed_;

    //! Boolean denoting whether spherical harmonic coefficients are normalized
    bool areCoefficientsNormalized_;

    //! List of integer powers of (distance / equatorial radius of body 1)
    std::vector< double > radius1Powers_;

    //! List of integer powers of (distance / equatorial radius of body 2)
    std::vector< double > radius2Powers_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
