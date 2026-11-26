/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Marsden, B. G., Sekanina, Z., Yeomans, D. K. (1973), Comets and nongravitational forces. V,
 *          The Astronomical Journal, 78, 211.
 */

#ifndef TUDAT_MARSDENACCELERATIONPARTIAL_H
#define TUDAT_MARSDENACCELERATIONPARTIAL_H

#include "tudat/astro/basic_astro/marsdenAcceleration.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Partial of the Marsden acceleration (inertial frame) w.r.t. the position of the accelerated body.
/*!
 *  a = R_rsw( r, v ) * [A1, A2, A3]^T * g( |r - v * tau| / AU ).
 *  \param relativeState Cartesian state of accelerated body w.r.t. the Sun (inertial frame) [m, m/s].
 *  \param A1 Radial coefficient [m/s^2].
 *  \param A2 Transverse coefficient [m/s^2].
 *  \param A3 Normal coefficient [m/s^2].
 *  \param alpha, r0, m, n, k Parameters of g(r).
 *  \param tau Time offset between perihelion and peak activity [s].
 */

Eigen::Matrix3d calculatePartialOfMarsdenAccelerationWrtPositionOfAcceleratedBody( const Eigen::Vector6d& relativeState,
                                                                                   const double A1,
                                                                                   const double A2,
                                                                                   const double A3,
                                                                                   const double alpha,
                                                                                   const double r0,
                                                                                   const double m,
                                                                                   const double n,
                                                                                   const double k,
                                                                                   const double tau );

//! Partial of the Marsden acceleration (inertial frame) w.r.t. the velocity of the accelerated body.
Eigen::Matrix3d calculatePartialOfMarsdenAccelerationWrtVelocityOfAcceleratedBody( const Eigen::Vector6d& relativeState,
                                                                                   const double A1,
                                                                                   const double A2,
                                                                                   const double A3,
                                                                                   const double alpha,
                                                                                   const double r0,
                                                                                   const double m,
                                                                                   const double n,
                                                                                   const double k,
                                                                                   const double tau );

//! Class to calculate the partials of the Marsden acceleration w.r.t. parameters and states.
class MarsdenAccelerationPartial : public AccelerationPartial
{
public:
    MarsdenAccelerationPartial( const std::shared_ptr< basic_astrodynamics::MarsdenAcceleration > marsdenAcceleration,
                                const std::string acceleratedBody,
                                const std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, marsdenAcceleration, basic_astrodynamics::marsden_acceleration ),
        marsdenAcceleration_( marsdenAcceleration )
    {
        currentPartialWrtPosition_.setZero( );
        currentPartialWrtVelocity_.setZero( );
    }

    //! Partial of the acceleration w.r.t. the position of body undergoing acceleration (added to partialMatrix block).
    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 )
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

    //! Partial of the acceleration w.r.t. the position of body exerting acceleration (= minus the accelerated-body partial).
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 0 )
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

    //! Partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 )
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

    //! Partial of the acceleration w.r.t. the velocity of body exerting acceleration (= minus the accelerated-body partial).
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 0 )
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

    //! No dependency on non-translational integrated states is implemented.
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& stateReferencePoint,
                                                                     const propagators::IntegratedStateType integratedStateType )
    {
        return 0;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter (A1, A2, A3, tau).
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! No vector parameters influence the Marsden acceleration.
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for updating the state partials to the current time.
    void update( const double currentTime = TUDAT_NAN )
    {
        marsdenAcceleration_->updateMembers( currentTime );

        if( !( currentTime_ == currentTime ) )
        {
            currentPartialWrtPosition_ =
                    calculatePartialOfMarsdenAccelerationWrtPositionOfAcceleratedBody( marsdenAcceleration_->getCurrentRelativeState( ),
                                                                                       marsdenAcceleration_->getA1( ),
                                                                                       marsdenAcceleration_->getA2( ),
                                                                                       marsdenAcceleration_->getA3( ),
                                                                                       marsdenAcceleration_->getAlpha( ),
                                                                                       marsdenAcceleration_->getR0( ),
                                                                                       marsdenAcceleration_->getM( ),
                                                                                       marsdenAcceleration_->getN( ),
                                                                                       marsdenAcceleration_->getK( ),
                                                                                       marsdenAcceleration_->getTau( ) );

            currentPartialWrtVelocity_ =
                    calculatePartialOfMarsdenAccelerationWrtVelocityOfAcceleratedBody( marsdenAcceleration_->getCurrentRelativeState( ),
                                                                                       marsdenAcceleration_->getA1( ),
                                                                                       marsdenAcceleration_->getA2( ),
                                                                                       marsdenAcceleration_->getA3( ),
                                                                                       marsdenAcceleration_->getAlpha( ),
                                                                                       marsdenAcceleration_->getR0( ),
                                                                                       marsdenAcceleration_->getM( ),
                                                                                       marsdenAcceleration_->getN( ),
                                                                                       marsdenAcceleration_->getK( ),
                                                                                       marsdenAcceleration_->getTau( ) );

            currentTime_ = currentTime;
        }
    }

protected:
    //! Partial w.r.t. radial coefficient A1 (= R_hat * g(r)), taken from the acceleration model.
    void wrtA1( Eigen::MatrixXd& partial )
    {
        partial = marsdenAcceleration_->getCurrentPartialWrtA1( );
    }

    //! Partial w.r.t. transverse coefficient A2 (= S_hat * g(r)).
    void wrtA2( Eigen::MatrixXd& partial )
    {
        partial = marsdenAcceleration_->getCurrentPartialWrtA2( );
    }

    //! Partial w.r.t. normal coefficient A3 (= W_hat * g(r)).
    void wrtA3( Eigen::MatrixXd& partial )
    {
        partial = marsdenAcceleration_->getCurrentPartialWrtA3( );
    }

    //! Partial w.r.t. the time offset tau (computed in the .cpp, as it is not linear in tau).
    void wrtTau( Eigen::MatrixXd& tauPartial );

    const std::shared_ptr< basic_astrodynamics::MarsdenAcceleration > marsdenAcceleration_;

    //! Current partial w.r.t. position of body undergoing acceleration (= -partial w.r.t. position of body exerting acceleration).
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! Current partial w.r.t. velocity of body undergoing acceleration.
    Eigen::Matrix3d currentPartialWrtVelocity_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_MARSDENACCELERATIONPARTIAL_H