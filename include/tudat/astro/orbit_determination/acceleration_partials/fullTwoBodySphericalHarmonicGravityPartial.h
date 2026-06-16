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

#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class for analytical partial derivatives of the full two-body spherical-harmonic acceleration.
/*!
 * Evaluates derivatives of the acceleration model based on Dirkx et al. (2019) effective-coefficient
 * formulation (Eqs. (47)-(49)) and translational dynamics expression (Eq. (55)).
 */
class FullTwoBodySphericalHarmonicsGravityPartial : public AccelerationPartial
{
public:
    //! Constructor.
    FullTwoBodySphericalHarmonicsGravityPartial(
            const std::string& acceleratedBody,
            const std::string& acceleratingBody,
            const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > accelerationModel );

    ~FullTwoBodySphericalHarmonicsGravityPartial( ) {}

    //! Update all cached partial terms to the current model state.
    void update( const double currentTime = TUDAT_NAN ) override;

    //! Insert partial w.r.t. position of body undergoing acceleration.
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

    //! Insert partial w.r.t. velocity of body undergoing acceleration (zero for this conservative model).
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

    //! Insert partial w.r.t. position of body exerting acceleration.
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

    //! Insert partial w.r.t. velocity of body exerting acceleration (zero for this conservative model).
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

    //! Retrieve partial function for scalar parameters (none implemented for this model).
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Retrieve partial function for vector parameters (spherical-harmonic blocks).
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    const Eigen::Matrix3d& getCurrentBodyFixedPartialWrtPosition( ) const
    {
        return currentBodyFixedPartialWrtPosition_;
    }

    const std::vector< Eigen::Matrix< double, 3, 2 > >& getCurrentBodyFixedPartialsWrtEffectiveCoefficients( ) const
    {
        return currentBodyFixedPartialsWrtEffectiveCoefficients_;
    }

    const std::vector< Eigen::Matrix< double, 3, 2 > >& getCurrentPartialsWrtEffectiveCoefficients( ) const
    {
        return currentPartialsWrtEffectiveCoefficients_;
    }

    const std::vector< Eigen::Matrix2d >& getCurrentEffectiveCoefficientsWrtTransformedBody2Coefficients( ) const
    {
        return currentEffectiveCoefficientsWrtTransformedBody2Coefficients_;
    }

    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& getCoefficientCombinationsToUse( ) const
    {
        return coefficientCombinationsToUse_;
    }

    //! Convenience wrapper used by torque partials to get transformed body-2 coefficient derivatives.
    void calculateCurrentTransformedBody2CoefficientPartials( const int degree,
                                                              const int order,
                                                              const bool wrtCosineCoefficient,
                                                              Eigen::MatrixXd& transformedCosinePartials,
                                                              Eigen::MatrixXd& transformedSinePartials )
    {
        updateCurrentTransformedBody2CoefficientPartials(
                degree, order, wrtCosineCoefficient, transformedCosinePartials, transformedSinePartials );
    }

private:
    //! Update partial of acceleration w.r.t. Cartesian relative position.
    void updateCurrentPositionPartial( );

    //! Update partials of acceleration w.r.t. effective coefficients.
    void updateCurrentPartialsWrtEffectiveCoefficients( );

    //! Update derivatives of transformed body-2 coefficients w.r.t. a single original body-2 coefficient.
    void updateCurrentTransformedBody2CoefficientPartials( const int degree,
                                                           const int order,
                                                           const bool wrtCosineCoefficient,
                                                           Eigen::MatrixXd& transformedCosinePartials,
                                                           Eigen::MatrixXd& transformedSinePartials );

    //! Partial w.r.t. cosine coefficient block of body 1.
    void wrtCosineCoefficientBlockOfBody1( const std::vector< std::pair< int, int > >& blockIndices, Eigen::MatrixXd& partialMatrix );

    //! Partial w.r.t. sine coefficient block of body 1.
    void wrtSineCoefficientBlockOfBody1( const std::vector< std::pair< int, int > >& blockIndices, Eigen::MatrixXd& partialMatrix );

    //! Partial w.r.t. cosine coefficient block of body 2.
    void wrtCosineCoefficientBlockOfBody2( const std::vector< std::pair< int, int > >& blockIndices, Eigen::MatrixXd& partialMatrix );

    //! Partial w.r.t. sine coefficient block of body 2.
    void wrtSineCoefficientBlockOfBody2( const std::vector< std::pair< int, int > >& blockIndices, Eigen::MatrixXd& partialMatrix );

    std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > accelerationModel_;

    std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;

    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    std::vector< Eigen::Matrix< double, 3, 2 > > currentPartialsWrtEffectiveCoefficients_;

    std::vector< Eigen::Matrix< double, 3, 2 > > currentBodyFixedPartialsWrtEffectiveCoefficients_;

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

    Eigen::MatrixXd body2BasisCosineCoefficients_;
    Eigen::MatrixXd body2BasisSineCoefficients_;
    Eigen::MatrixXd transformedCosineBody2CoefficientPartialsScratch_;
    Eigen::MatrixXd transformedSineBody2CoefficientPartialsScratch_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_MUTUALEXTENDEDBODYSPHERICALHARMONICGRAVITYPARTIAL_H
