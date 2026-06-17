/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MUTUALFORCEPOTENTIAL_H
#define TUDAT_MUTUALFORCEPOTENTIAL_H

#include <boost/math/special_functions/factorials.hpp>
#include <tuple>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <functional>
#include <memory>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/basic/sphericalHarmonics.h"

#include "tudat/math/basic/sphericalHarmonicTransformations.h"

namespace tudat
{

namespace gravitation
{

//! Function to get maximum degrees of used for the spherical harmonic expansions of the two bodies
/*!
 *  Function to get maximum degrees of used for the spherical harmonic expansions of the two bodies
 *  \param coefficientCombinationsToUse st of degrees/orders that are to be used for the series expansion.
 *  Each tuple contains: (degree of body 1, order of body 1, degree of body 2, order of body 2)
 *  \return Maximum degree of body 1 and body 2 (as a pair)
 */
std::pair< int, int > getMaximumDegrees(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse );

//! Function to compute cross-body normalization terms for mutual two-body potential
/*!
 * Function to compute cross-body normalization term gamma for mutual two-body potential, according to formulation of Compere &
 * Lemaitre (2014)
 * \param l Parameter l used by Compere & Lemaitre (2014); degree of body 1
 * \param m Parameter m used by Compere & Lemaitre (2014); order of body 1
 * \param j Parameter j used by Compere & Lemaitre (2014); degree of body 2
 * \param k Parameter k used by Compere & Lemaitre (2014); order of body 2
 * \return Term gamma for mutual two-body potential, according to formulation of Compere & Lemaitre (2014)
 */
double getGammaCoefficientForMutualForcePotential( const int l, const int m, const int j, const int k );

//! Function to compute cross-body normalization terms for mutual two-body potential, for unnormalized or fully normalized
//! coefficients
/*!
 * Function to compute cross-body normalization terms for mutual two-body potential, for unnormalized or fully normalized
 * coefficients
 * \param degree1 Degree of body 1
 * \param order1 Order of body 1
 * \param degree2 Degree of body 2
 * \param order2 Order of body 2
 * \param areCoefficientsNormalized Boolean denoting whether the coefficients are fully normalized or unnormalized
 * \return Cross-body normalization terms for mutual two-body potential
 */
double getMutualPotentialEffectiveCoefficientMultiplier( const int degree1,
                                                         const int order1,
                                                         const int degree2,
                                                         const int order2,
                                                         const bool areCoefficientsNormalized );

//! Function to compute single-term in two-body potential, from effective one-body formulation of Dirkx et al. (2018)
/*!
 * Function to compute single-term in two-body potential, from effective one-body formulation of Dirkx et al. (2018),
 * omitting the radius power term, and the common multiplier for all terms
 * \param effectiveCosineCoefficient Effective one-body cosine coefficient
 * \param effectiveSineCoefficient Effective one-body sine coefficient
 * \param sphericalHarmonicsCache Cache object used to pre-compute Legendre polynomials and other spherical harmonic terms
 * \param degree1 Degree of body 1
 * \param order1 Order of body 1
 * \param degree2 Degree of body 2
 * \param order2 Order of body 2
 * \return Single-term in effective one-body formulation, omitting the radius power term, and the common multiplier for all terms
 */
double computeSingleMutualForcePotentialTerm( const double effectiveCosineCoefficient,
                                              const double effectiveSineCoefficient,
                                              std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
                                              const int degreeOfBody1,
                                              const int orderOfBody1,
                                              const int degreeOfBody2,
                                              const int orderOfBody2 );

double computeMutualForcePotential(
        const Eigen::Vector3d& bodyFixedPosition,
        const double effectiveGravitationalParameterOfBody1,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const int maximumDegreeOfBody1,
        const int maximumDegreeOfBody2,
        const std::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const std::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Compute full two-body mutual acceleration for normalized coefficients.
/*!
 * Evaluates the mutual-potential gradient using the effective one-body mapping from Dirkx et al. (2019),
 * combining effective coefficients (Eqs. (47)-(48)) in the potential form of Eq. (49), as used in the
 * translational equation Eq. (55).
 */
Eigen::Vector3d computeGeodesyNormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const std::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const std::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const int maximumDegree1,
        const int maximumDegree2,
        const int maximumEvaluationDegree,
        const std::vector< double >& radius1Powers,
        const std::vector< double >& radius2Powers,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Compute full two-body mutual acceleration for unnormalized coefficients.
/*!
 * Unnormalized counterpart of computeGeodesyNormalizedMutualGravitationalAccelerationSum.
 */
Eigen::Vector3d computeUnnormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const std::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const std::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

void computePartialDerivativesOfPotentialComponentsWrtFullCoefficients(
        std::vector< Eigen::Matrix< double, 1, 2 > >& potentialComponentsWrtFullCoefficients,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const double distance,
        const std::vector< double >& radius1Powers,
        const std::vector< double >& radius2Powers,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::function< int( const int, const int, const int, const int ) > effectiveIndexFunction );

inline double getSigmaSignFunction( const int order )
{
    // Implements sigma_m from Dirkx et al. (2019), Eq. (22), for signed-order handling.
    return ( ( order >= 0 ) ? ( 1.0 ) : ( ( std::abs( order ) % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) ) );
}

//! Map a compact sign-combination case index to the signed (m1,m2) pair.
/*!
 * Generates the four signed-order combinations needed when the coefficient combination list stores non-negative
 * orders only: (+m1,+m2), (-m1,+m2), (+m1,-m2), (-m1,-m2). Terms that are not distinct for m=0 are skipped.
 * This centralizes the signed-order branching used by the full two-body acceleration/potential evaluation in the
 * Eq. (49) summation and the Eq. (47)-(48) effective-coefficient construction.
 */
inline bool getSignedOrdersForCombinationCase( const int combinationCase,
                                               const int orderOfBody1,
                                               const int orderOfBody2,
                                               int& signedOrderOfBody1,
                                               int& signedOrderOfBody2 )
{
    // This switch enumerates the signed-order variants needed by the Eq. (49) summation
    // and the effective-coefficient combinations in Eqs. (47)-(48).
    switch( combinationCase )
    {
        case 0:
            signedOrderOfBody1 = orderOfBody1;
            signedOrderOfBody2 = orderOfBody2;
            return true;
        case 1:
            signedOrderOfBody1 = -orderOfBody1;
            signedOrderOfBody2 = orderOfBody2;
            return ( orderOfBody1 != 0 );
        case 2:
            signedOrderOfBody1 = orderOfBody1;
            signedOrderOfBody2 = -orderOfBody2;
            return ( orderOfBody2 != 0 );
        case 3:
            signedOrderOfBody1 = -orderOfBody1;
            signedOrderOfBody2 = -orderOfBody2;
            return ( orderOfBody1 != 0 && orderOfBody2 != 0 );
        default:
            return false;
    }
}

class EffectiveMutualSphericalHarmonicsField
{
public:
    EffectiveMutualSphericalHarmonicsField(
            const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
            const std::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody1,
            const std::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody1,
            const std::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody2,
            const std::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody2,
            const std::function< double( ) > gravitationalParameterFunction,
            const double equatorialRadiusOfBody1,
            const double equatorialRadiusOfBody2,
            const bool areCoefficientsNormalized = 1 ):
        coefficientCombinationsToUse_( coefficientCombinationsToUse ),
        cosineCoefficientFunctionOfBody1_( cosineCoefficientFunctionOfBody1 ),
        sineCoefficientFunctionOfBody1_( sineCoefficientFunctionOfBody1 ),
        cosineCoefficientFunctionOfBody2_( cosineCoefficientFunctionOfBody2 ),
        sineCoefficientFunctionOfBody2_( sineCoefficientFunctionOfBody2 ),
        gravitationalParameterFunction_( gravitationalParameterFunction ), equatorialRadiusOfBody1_( equatorialRadiusOfBody1 ),
        equatorialRadiusOfBody2_( equatorialRadiusOfBody2 ), areCoefficientsNormalized_( areCoefficientsNormalized ),
        calculatePartials_( false )
    {
        std::pair< int, int > maximumDegrees = getMaximumDegrees( coefficientCombinationsToUse_ );
        maximumDegree1_ = maximumDegrees.first;
        maximumDegree2_ = maximumDegrees.second;

        effectiveCosineCoefficients_.resize( getTotalVectorSize( ) );
        effectiveSineCoefficients_.resize( getTotalVectorSize( ) );
        multipliers_.resize( getTotalVectorSize( ) );
        transformationCache_ = std::make_shared< basic_mathematics::SphericalHarmonicTransformationCache >(
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

    void getCurrentEffectiveCoefficients( const int degree1,
                                          const int order1,
                                          const int degree2,
                                          const int order2,
                                          const int effectiveIndex,
                                          double& cosineCoefficient,
                                          double& sineCoefficient );

    void computeCurrentEffectiveCoefficients( const Eigen::Quaterniond coefficientRotationQuaterion );

    void computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients( const Eigen::MatrixXd& transformedCosineCoefficients,
                                                                               const Eigen::MatrixXd& transformedSineCoefficients );

    void updateEffectiveMutualPotential( );

    int getEffectiveIndex( const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return degreeOfBody1 +
                ( maximumDegree1_ + 1 ) *
                ( orderOfBody1 + maximumDegree1_ +
                  ( 2 * maximumDegree1_ + 1 ) * ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * ( orderOfBody2 + maximumDegree2_ ) ) );
    }

    double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition,
                                      std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
    {
        return computeMutualForcePotential( bodyFixedPosition,
                                            gravitationalParameterFunction_( ),
                                            equatorialRadiusOfBody1_,
                                            equatorialRadiusOfBody2_,
                                            maximumDegree1_,
                                            maximumDegree2_,
                                            std::bind( &EffectiveMutualSphericalHarmonicsField::getEffectiveCosineCoefficient,
                                                       this,
                                                       std::placeholders::_1,
                                                       std::placeholders::_2,
                                                       std::placeholders::_3,
                                                       std::placeholders::_4 ),
                                            std::bind( &EffectiveMutualSphericalHarmonicsField::getEffectiveSineCoefficient,
                                                       this,
                                                       std::placeholders::_1,
                                                       std::placeholders::_2,
                                                       std::placeholders::_3,
                                                       std::placeholders::_4 ),
                                            coefficientCombinationsToUse_,
                                            sphericalHarmonicsCache );
    }

    double getEffectiveCosineCoefficient( const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveCosineCoefficients_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    double getEffectiveSineCoefficient( const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveSineCoefficients_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    std::shared_ptr< basic_mathematics::SphericalHarmonicTransformationCache > getTransformationCache( )
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

    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& getCoefficientCombinationsToUse( ) const
    {
        return coefficientCombinationsToUse_;
    }

    const Eigen::MatrixXd& getCosineCoefficientsOfBody1( ) const
    {
        return cosineCoefficientsOfBody1_;
    }

    const Eigen::MatrixXd& getSineCoefficientsOfBody1( ) const
    {
        return sineCoefficientsOfBody1_;
    }

    const Eigen::MatrixXd& getCosineCoefficientsOfBody2( ) const
    {
        return cosineCoefficientsOfBody2_;
    }

    const Eigen::MatrixXd& getSineCoefficientsOfBody2( ) const
    {
        return sineCoefficientsOfBody2_;
    }

    double getMultiplier( const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return multipliers_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    void computePartialsOfFullCoefficientsWrtTransformedCoefficients(
            std::vector< Eigen::Matrix2d >& fullCoefficientsWrtBody2CoefficientsList );

    const Eigen::MatrixXd& getTransformedCosineCoefficientsOfBody2( ) const
    {
        return transformedCosineCoefficientsOfBody2_;
    }

    const Eigen::MatrixXd& getTransformedSineCoefficientsOfBody2( ) const
    {
        return transformedSineCoefficientsOfBody2_;
    }

private:
    double getEffectiveCoefficientMultiplier( const int degree1, const int order1, const int degree2, const int order2 );

    void initializeMultipliers( );

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    std::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody1_;
    std::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody1_;
    std::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody2_;
    std::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody2_;

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

    std::shared_ptr< basic_mathematics::SphericalHarmonicTransformationCache > transformationCache_;

    std::function< double( ) > gravitationalParameterFunction_;
    double equatorialRadiusOfBody1_;
    double equatorialRadiusOfBody2_;

    bool areCoefficientsNormalized_;

    bool calculatePartials_;

    int maximumDegree1_;

    int maximumDegree2_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_MUTUALFORCEPOTENTIAL_H
