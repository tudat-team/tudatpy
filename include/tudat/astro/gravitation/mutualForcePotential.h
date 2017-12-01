#ifndef MUTUALFORCEPOTENTIAL_H
#define MUTUALFORCEPOTENTIAL_H

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"

namespace tudat
{

namespace gravitation
{

std::pair< int, int > getMaximumDegrees(
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse );

double getGammaCoefficientForMutualForcePotential(
        const int l, const int m, const int j, const int k );

double computeSingleMutualForcePotentialTerm(
        const double effectiveCosineCoefficient,
        const double effectiveSineCoefficient,
        const double polynomialParameter,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const int degreeOfBody1,
        const int orderOfBody1,
        const int degreeOfBody2,
        const int orderOfBody2 );

double getMutualPotentialEffectiveCoefficientMultiplier(
        const int degree1, const int order1, const int degree2, const int order2, const bool areCoefficientsNormalized );

double computeMutualForcePotential(
        const Eigen::Vector3d& bodyFixedPosition,
        const double effectiveGravitationalParameterOfBody1,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const int maximumDegreeOfBody1,
        const int maximumDegreeOfBody2,
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
Eigen::Vector3d computeGeodesyNormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const int maximumDegree1,
        const int maximumDegree2,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
Eigen::Vector3d computeUnnormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

void computePartialDerivativesOfPotentialComponentsWrtFullCoefficients(
        std::vector< Eigen::Matrix< double, 1, 2 > >& potentialComponentsWrtFullCoefficients,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const double distance,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const boost::function< int( const int, const int, const int, const int )> effectiveIndexFunction );

inline double getSigmaSignFunction( const int order )
{
    return ( ( order >= 0 ) ? ( 1.0 ) : ( ( std::abs( order ) % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) ) );
}

class EffectiveMutualSphericalHarmonicsField
{
public:
    EffectiveMutualSphericalHarmonicsField(
            const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
            const boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody1,
            const boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody1,
            const boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody2,
            const boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody2,
            const boost::function< double( ) > gravitationalParameterFunction,
            const double equatorialRadiusOfBody1,
            const double equatorialRadiusOfBody2,
            const bool areCoefficientsNormalized = 1 ):
        coefficientCombinationsToUse_( coefficientCombinationsToUse ),
        cosineCoefficientFunctionOfBody1_( cosineCoefficientFunctionOfBody1 ),
        sineCoefficientFunctionOfBody1_( sineCoefficientFunctionOfBody1 ),
        cosineCoefficientFunctionOfBody2_( cosineCoefficientFunctionOfBody2 ),
        sineCoefficientFunctionOfBody2_( sineCoefficientFunctionOfBody2 ),
        gravitationalParameterFunction_( gravitationalParameterFunction ),
        equatorialRadiusOfBody1_( equatorialRadiusOfBody1 ),
        equatorialRadiusOfBody2_( equatorialRadiusOfBody2 ),
        areCoefficientsNormalized_( areCoefficientsNormalized ),
        calculatePartials_( false )
    {
        std::pair< int, int > maximumDegrees = getMaximumDegrees( coefficientCombinationsToUse_ );
        maximumDegree1_ = maximumDegrees.first;
        maximumDegree2_ = maximumDegrees.second;

        effectiveCosineCoefficients_.resize(
                    getTotalVectorSize( ) );
        effectiveSineCoefficients_.resize(
                    getTotalVectorSize( ) );
        multipliers_.resize(
                    getTotalVectorSize( ) );
        transformationCache_ = boost::make_shared< basic_mathematics::SphericalHarmonicTransformationCache >(
                    cosineCoefficientFunctionOfBody1( ).rows( ) + cosineCoefficientFunctionOfBody2( ).rows( ),
                    cosineCoefficientFunctionOfBody1( ).cols( ) + cosineCoefficientFunctionOfBody2( ).cols( ) );
        initializeMultipliers( );
    }

    int getMaximumDegree1( )
    {
        return maximumDegree1_;
    }

    int getMaximumDegree2( )
    {
        return maximumDegree2_;
    }

    int getTotalVectorSize( )
    {
        return ( 2 * maximumDegree1_ + 1 ) * ( maximumDegree1_ + 1 ) * ( 2 * maximumDegree2_ + 1 ) * ( maximumDegree2_ + 1 );
    }

    void getCurrentEffectiveCoefficients(
            const int degree1, const int order1, const int degree2, const int order2, const int effectiveIndex,
            double& cosineCoefficient, double& sineCoefficient );

    void computeCurrentEffectiveCoefficients( const Eigen::Quaterniond coefficientRotationQuaterion );

    void computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
            const Eigen::MatrixXd& transformedCosineCoefficients,
            const Eigen::MatrixXd& transformedSineCoefficients );

    void updateEffectiveMutualPotential( );

    int getEffectiveIndex(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return degreeOfBody1 + ( maximumDegree1_ + 1 ) * ( orderOfBody1 + maximumDegree1_ +( 2 * maximumDegree1_ + 1 ) *
                                                ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * ( orderOfBody2 + maximumDegree2_ ) ) );
    }

    double getGravitationalPotential(
            const Eigen::Vector3d& bodyFixedPosition,
            boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
    {

        return computeMutualForcePotential(
                    bodyFixedPosition, gravitationalParameterFunction_( ), equatorialRadiusOfBody1_, equatorialRadiusOfBody2_,
                    maximumDegree1_, maximumDegree2_,
                    boost::bind(
                        &EffectiveMutualSphericalHarmonicsField::getEffectiveCosineCoefficient,
                        this, _1, _2, _3, _4 ),
                    boost::bind(
                        &EffectiveMutualSphericalHarmonicsField::getEffectiveSineCoefficient,
                        this, _1, _2, _3, _4 ), coefficientCombinationsToUse_,
                    sphericalHarmonicsCache );
    }


    double getEffectiveCosineCoefficient(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveCosineCoefficients_[ degreeOfBody1 + ( maximumDegree1_ + 1 ) * ( orderOfBody1 + maximumDegree1_+ ( 2 * maximumDegree1_ + 1 ) * ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * ( maximumDegree2_ + orderOfBody2 ) ) ) ];
    }

    double getEffectiveSineCoefficient(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        std::cout<<"Getting sine: "<<degreeOfBody1<<" "<<orderOfBody1<<" "<<degreeOfBody2<<" "<<orderOfBody2<<" "<<
                   degreeOfBody1 + ( maximumDegree1_ + 1 ) * ( orderOfBody1 + ( 2 * maximumDegree1_ + 1 ) * ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * orderOfBody2 ) )<<" "<<
                   effectiveSineCoefficients_[ degreeOfBody1 + ( maximumDegree1_ + 1 ) * ( orderOfBody1 + ( 2 * maximumDegree1_ + 1 ) * ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * orderOfBody2 ) ) ]<<std::endl;
        return effectiveSineCoefficients_[ degreeOfBody1 + ( maximumDegree1_ + 1 ) * ( orderOfBody1 + ( 2 * maximumDegree1_ + 1 ) * ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * orderOfBody2 ) ) ];
    }

    boost::shared_ptr< basic_mathematics::SphericalHarmonicTransformationCache > getTransformationCache( )
    {
        return transformationCache_;
    }

    void setCalculatePartials( )
    {
        transformedCosineCoefficientsOfBody2Partials_.resize( 3 );
        transformedSineCoefficientsOfBody2Partials_.resize( 3 );

        transformationCache_->setUpdatePartials( );

        calculatePartials_ = true;
    }

    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getCoefficientCombinationsToUse( )
    {
        return coefficientCombinationsToUse_;
    }

    Eigen::MatrixXd getCosineCoefficientsOfBody1( )
    {
        return cosineCoefficientsOfBody1_;
    }

    Eigen::MatrixXd getSineCoefficientsOfBody1( )
    {
        return sineCoefficientsOfBody1_;
    }

    double getMultiplier( const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return multipliers_[ getEffectiveIndex(
                    degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ] ;
    }

//    void getTransformedCoefficientPartialsWrtEulerAngles(
//            std::vector< Eigen::MatrixXd >& transformedCosineCoefficientsOfBody2Partials,
//            std::vector< Eigen::MatrixXd >& transformedSineCoefficientsOfBody2Partials )
//    {
//        transformationCache_->getPartialDerivativesOfTransformedCoefficientsWrtEulerAngles(
//                    cosineCoefficientsOfBody2_,
//                    sineCoefficientsOfBody2_,
//                    transformedCosineCoefficientsOfBody2Partials_,
//                    transformedSineCoefficientsOfBody2Partials_,
//                    areCoefficientsNormalized_ );

//        transformedCosineCoefficientsOfBody2Partials = transformedCosineCoefficientsOfBody2Partials_;
//        transformedSineCoefficientsOfBody2Partials = transformedSineCoefficientsOfBody2Partials_;

//    }

    void computePartialsOfFullCoefficientsWrtTransformedCoefficients(
            std::vector< Eigen::Matrix2d >& fullCoefficientsWrtBody2CoefficientsList );

    Eigen::MatrixXd getTransformedCosineCoefficientsOfBody2( )
    {
        return transformedCosineCoefficientsOfBody2_;
    }

    Eigen::MatrixXd getTransformedSineCoefficientsOfBody2( )
    {
        return transformedSineCoefficientsOfBody2_;
    }

private:

    double getEffectiveCoefficientMultiplier(
            const int degree1, const int order1, const int degree2, const int order2 );

    void initializeMultipliers( );

    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody1_;
    boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody1_;
    boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody2_;
    boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody2_;


    Eigen::MatrixXd cosineCoefficientsOfBody1_;
    Eigen::MatrixXd sineCoefficientsOfBody1_;
    Eigen::MatrixXd cosineCoefficientsOfBody2_;
    Eigen::MatrixXd sineCoefficientsOfBody2_;

    Eigen::MatrixXd transformedCosineCoefficientsOfBody2_;
    Eigen::MatrixXd transformedSineCoefficientsOfBody2_;

    std::vector< Eigen::MatrixXd > transformedCosineCoefficientsOfBody2Partials_;
    std::vector< Eigen::MatrixXd > transformedSineCoefficientsOfBody2Partials_;

    std::vector< double > effectiveCosineCoefficients_;
    std::vector< double > effectiveSineCoefficients_;

    std::vector< double > multipliers_;

    boost::shared_ptr< basic_mathematics::SphericalHarmonicTransformationCache > transformationCache_;

    boost::function< double( ) > gravitationalParameterFunction_;
    double equatorialRadiusOfBody1_;
    double equatorialRadiusOfBody2_;

    bool areCoefficientsNormalized_;

    bool calculatePartials_;

    int maximumDegree1_;

    int maximumDegree2_;

};


}

}

#endif // MUTUALFORCEPOTENTIAL_H
