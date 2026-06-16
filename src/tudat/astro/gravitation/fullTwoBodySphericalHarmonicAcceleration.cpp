#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"

#include <algorithm>

namespace tudat
{

namespace gravitation
{

FullTwoBodySphericalHarmonicAcceleration::FullTwoBodySphericalHarmonicAcceleration(
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
        const bool areCoefficientsNormalized ):
    positionOfBody1Function_( positionOfBody1Function ), positionOfBody2Function_( positionOfBody2Function ),
    gravitationalParameterFunction_( gravitationalParameterFunction ), equatorialRadiusOfBody1_( equatorialRadiusOfBody1 ),
    equatorialRadiusOfBody2_( equatorialRadiusOfBody2 ), coefficientCombinationsToUse_( coefficientCombinationsToUse ),
    toLocalFrameOfBody1Transformation_( toLocalFrameOfBody1Transformation ),
    toLocalFrameOfBody2Transformation_( toLocalFrameOfBody2Transformation ), isMutualAttractionUsed_( isMutualAttractionUsed ),
    areCoefficientsNormalized_( areCoefficientsNormalized )
{
    // Determine the maximum effective degree/order l,m used in the Eq. (49) summation of Dirkx et al. (2019),
    // with l=l1+l2 and m=|m1+m2| for each selected (l1,m1,l2,m2) interaction tuple.
    maximumDegree_ = 0;
    maximumOrder_ = 0;

    unsigned int degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse.at( i ) );
        orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse.at( i ) );
        degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse.at( i ) );
        orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse.at( i ) );

        if( degreeOfBody1 + degreeOfBody2 > maximumDegree_ )
        {
            maximumDegree_ = degreeOfBody1 + degreeOfBody2;
        }

        if( orderOfBody1 + orderOfBody2 > maximumOrder_ )
        {
            maximumOrder_ = orderOfBody1 + orderOfBody2;
        }
    }

    // The geodesy-normalized acceleration evaluation precomputes P_lm and dP_lm/dphi terms used in Eq. (55),
    // and requires access up to m=l+1 for each l due to derivative evaluation.
    maximumOrder_ = std::max( maximumOrder_, maximumDegree_ );

    // Initialize cache and effective-field object that converts two-body coefficients to effective one-body
    // coefficients (Dirkx et al. (2019), Eqs. (47)-(48)).
    sphericalHarmonicsCache_ = std::make_shared< basic_mathematics::SphericalHarmonicsCache >( maximumDegree_ + 1, maximumOrder_ + 1 );
    effectiveMutualPotentialField_ = std::make_shared< EffectiveMutualSphericalHarmonicsField >( coefficientCombinationsToUse_,
                                                                                                 cosineHarmonicCoefficientsOfBody1Function,
                                                                                                 sineHarmonicCoefficientsOfBody1Function,
                                                                                                 cosineHarmonicCoefficientsOfBody2Function,
                                                                                                 sineHarmonicCoefficientsOfBody2Function,
                                                                                                 gravitationalParameterFunction,
                                                                                                 equatorialRadiusOfBody1_,
                                                                                                 equatorialRadiusOfBody2_,
                                                                                                 areCoefficientsNormalized );
    effectiveCosineCoefficientFunction_ = std::bind( &EffectiveMutualSphericalHarmonicsField::getEffectiveCosineCoefficient,
                                                     effectiveMutualPotentialField_,
                                                     std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4 );
    effectiveSineCoefficientFunction_ = std::bind( &EffectiveMutualSphericalHarmonicsField::getEffectiveSineCoefficient,
                                                   effectiveMutualPotentialField_,
                                                   std::placeholders::_1,
                                                   std::placeholders::_2,
                                                   std::placeholders::_3,
                                                   std::placeholders::_4 );
    radius1Powers_.resize( effectiveMutualPotentialField_->getMaximumDegree1( ) + 1 );
    radius2Powers_.resize( effectiveMutualPotentialField_->getMaximumDegree2( ) + 1 );
}

void FullTwoBodySphericalHarmonicAcceleration::updateMembers( const double currentTime )
{
    if( !( currentTime == currentTime_ ) )
    {
        // Step 1: update current frame rotations.
        currentRotationFromInertialToBody1_ = toLocalFrameOfBody1Transformation_( );
        currentRotationFromBody2ToBody1_ = currentRotationFromInertialToBody1_ * toLocalFrameOfBody2Transformation_( ).inverse( );

        // Step 2: compute relative state and express it in body-1 frame, in which Eq. (49) is evaluated.
        currentRelativePosition_ = positionOfBody1Function_( ) - positionOfBody2Function_( );
        currentBodyFixedRelativePosition_ = currentRotationFromInertialToBody1_ * ( currentRelativePosition_ );

        // Step 3: transform body-2 coefficients to frame F1 and build effective coefficients used in Eq. (49)
        // through the Eq. (47)-(48) mapping.
        effectiveMutualPotentialField_->computeCurrentEffectiveCoefficients( currentRotationFromBody2ToBody1_ );

        // Step 4: precompute (R1/r)^l1 and (R2/r)^l2 factors from the radius substitutions in Dirkx Eq. (44).
        double currentDistance = currentRelativePosition_.norm( );
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
            // Step 5a: evaluate the Cartesian gradient of the effective potential in F1 for normalized coefficients.
            // This applies Eq. (55) with potential terms from Eq. (49).
            mutualPotentialGradient_ =
                    computeGeodesyNormalizedMutualGravitationalAccelerationSum( currentBodyFixedRelativePosition_,
                                                                                gravitationalParameterFunction_( ),
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
            // Step 5b: same as above, but for unnormalized coefficients.
            mutualPotentialGradient_ = computeUnnormalizedMutualGravitationalAccelerationSum( currentBodyFixedRelativePosition_,
                                                                                              gravitationalParameterFunction_( ),
                                                                                              equatorialRadiusOfBody1_,
                                                                                              equatorialRadiusOfBody2_,
                                                                                              effectiveCosineCoefficientFunction_,
                                                                                              effectiveSineCoefficientFunction_,
                                                                                              coefficientCombinationsToUse_,
                                                                                              sphericalHarmonicsCache_ );
        }

        // Step 6: output acceleration in inertial orientation.
        // The mutual potential gradient is assembled in body-1-fixed coordinates, so we always rotate it back
        // to inertial orientation for translational propagation consistency.
        currentAcceleration_ = currentRotationFromInertialToBody1_.inverse( ) * mutualPotentialGradient_;
        currentTime_ = currentTime;
    }
}

}  // namespace gravitation

}  // namespace tudat
