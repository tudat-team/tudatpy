/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITATIONALPARAMETERRATEACCELERATIONPARTIAL_H
#define TUDAT_GRAVITATIONALPARAMETERRATEACCELERATIONPARTIAL_H

#include "tudat/astro/gravitation/gravitationalParameterRateAcceleration.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to calculate the partials of the gravitational parameter rate acceleration w.r.t. parameters and states.
class GravitationalParameterRateAccelerationPartial : public AccelerationPartial
{
public:
    GravitationalParameterRateAccelerationPartial(
            const std::shared_ptr< gravitation::GravitationalParameterRateAcceleration > accelerationModel,
            const std::string acceleratedBody,
            const std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody,
                             acceleratingBody,
                             accelerationModel,
                             basic_astrodynamics::gravitational_parameter_rate_acceleration ),
        accelerationModel_( accelerationModel )
    {}

    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 )
    {
        partialMatrix.block( startRow, startColumn, 3, 3 ) += ( addContribution ? 1.0 : -1.0 ) * currentPartialWrtPosition_;
    }

    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 0 )
    {
        partialMatrix.block( startRow, startColumn, 3, 3 ) -= ( addContribution ? 1.0 : -1.0 ) * currentPartialWrtPosition_;
    }

    //! Partials w.r.t. the relative gravitational parameter rate and the gravitational parameter of the exerting body.
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        // The acceleration is the product of the two parameters and -(t - t0) r / |r|^3
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        if( parameter->getParameterName( ).second.first == acceleratingBody_ )
        {
            if( parameter->getParameterName( ).first == estimatable_parameters::relative_gravitational_parameter_rate )
            {
                partialFunction = [ this ]( Eigen::MatrixXd& partial ) {
                    partial = accelerationModel_->getGravitationalParameterFunction( )( ) * accelerationPerUnitParameterProduct( );
                };
            }
            else if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
            {
                partialFunction = [ this ]( Eigen::MatrixXd& partial ) {
                    partial = accelerationModel_->getRelativeGravitationalParameterRate( ) * accelerationPerUnitParameterProduct( );
                };
            }
        }
        return std::make_pair( partialFunction, partialFunction ? 1 : 0 );
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    void update( const double currentTime = TUDAT_NAN )
    {
        accelerationModel_->updateMembers( currentTime );
        if( !( currentTime_ == currentTime ) )
        {
            currentPartialWrtPosition_ = calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                    accelerationModel_->getCurrentRelativePosition( ),
                    accelerationModel_->getGravitationalParameterFunction( )( ) *
                            accelerationModel_->getRelativeGravitationalParameterRate( ) *
                            accelerationModel_->getCurrentTimeSinceReferenceEpoch( ) );
            currentTime_ = currentTime;
        }
    }

protected:
    Eigen::Vector3d accelerationPerUnitParameterProduct( )
    {
        const Eigen::Vector3d relativePosition = accelerationModel_->getCurrentRelativePosition( );
        return -accelerationModel_->getCurrentTimeSinceReferenceEpoch( ) * relativePosition / std::pow( relativePosition.norm( ), 3 );
    }

    std::shared_ptr< gravitation::GravitationalParameterRateAcceleration > accelerationModel_;

    Eigen::Matrix3d currentPartialWrtPosition_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_GRAVITATIONALPARAMETERRATEACCELERATIONPARTIAL_H
