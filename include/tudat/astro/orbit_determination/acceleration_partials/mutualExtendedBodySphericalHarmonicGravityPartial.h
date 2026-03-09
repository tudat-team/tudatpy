/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MUTUALEXTENDEDBODYSPHERICALHARMONICGRAVITYPARTIAL_H
#define TUDAT_MUTUALEXTENDEDBODYSPHERICALHARMONICGRAVITYPARTIAL_H

#include <tuple>

#include "tudat/astro/gravitation/mutualExtendedBodySphericalHarmonicAcceleration.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class for calculating partial derivatives of a mutual extended-body spherical harmonic acceleration.
class MutualExtendedBodySphericalHarmonicsGravityPartial : public AccelerationPartial
{
public:
    MutualExtendedBodySphericalHarmonicsGravityPartial(
            const std::string& acceleratedBody,
            const std::string& acceleratingBody,
            const std::shared_ptr< gravitation::MutualExtendedBodySphericalHarmonicAcceleration > accelerationModel );

    ~MutualExtendedBodySphericalHarmonicsGravityPartial( ) { }

    void update( const double currentTime = TUDAT_NAN ) override;

    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 ) override
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }

    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 3 ) override
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
    }

    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 0 ) override
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }

    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 3 ) override
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

private:
    void updateCurrentPositionPartial( );

    void updateCurrentPartialsWrtEffectiveCoefficients( );

    void updateCurrentTransformedBody2CoefficientPartials(
            const int degree,
            const int order,
            const bool wrtCosineCoefficient,
            Eigen::MatrixXd& transformedCosinePartials,
            Eigen::MatrixXd& transformedSinePartials );

    void wrtCosineCoefficientBlockOfBody1(
            const std::vector< std::pair< int, int > >& blockIndices,
            Eigen::MatrixXd& partialMatrix );

    void wrtSineCoefficientBlockOfBody1(
            const std::vector< std::pair< int, int > >& blockIndices,
            Eigen::MatrixXd& partialMatrix );

    void wrtCosineCoefficientBlockOfBody2(
            const std::vector< std::pair< int, int > >& blockIndices,
            Eigen::MatrixXd& partialMatrix );

    void wrtSineCoefficientBlockOfBody2(
            const std::vector< std::pair< int, int > >& blockIndices,
            Eigen::MatrixXd& partialMatrix );

    std::shared_ptr< gravitation::MutualExtendedBodySphericalHarmonicAcceleration > accelerationModel_;

    std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;

    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    std::vector< int > effectiveIndicesForCoefficientCombinations_;

    std::vector< Eigen::Matrix< double, 3, 2 > > currentPartialsWrtEffectiveCoefficients_;

    std::vector< Eigen::Matrix2d > currentEffectiveCoefficientsWrtTransformedBody2Coefficients_;

    Eigen::Matrix3d currentPartialWrtPosition_;

    Eigen::Matrix3d currentBodyFixedPartialWrtPosition_;

    Eigen::Matrix3d currentPartialWrtVelocity_;

    Eigen::Matrix3d currentRotationToBodyFixedFrame_;

    Eigen::Matrix3d currentRotationToInertialFrame_;

    Eigen::Vector3d currentBodyFixedRelativePosition_;

    Eigen::Vector3d currentSphericalPosition_;

    Eigen::Vector3d currentBodyFixedSphericalGradient_;

    double currentDistance_;

    double currentGravitationalParameter_;

    std::vector< double > currentRadius1Powers_;

    std::vector< double > currentRadius2Powers_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_MUTUALEXTENDEDBODYSPHERICALHARMONICGRAVITYPARTIAL_H
