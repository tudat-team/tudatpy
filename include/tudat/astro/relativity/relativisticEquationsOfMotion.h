/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTICEQUATIONSOFMOTION_H
#define TUDAT_RELATIVISTICEQUATIONSOFMOTION_H

#include <vector>

#include <memory> 

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/relativity/metric.h"


namespace tudat
{

namespace relativity
{

Eigen::Vector4d evaluateFourAcceleration(
        const std::vector< Eigen::Matrix4d >& christoffelSymbols, const Eigen::Vector4d& fourVelocity );

Eigen::Vector4d evaluateFourAcceleration(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric );

Eigen::Vector4d evaluateFourAccelerationWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 8, 1 >& currentState,
        const double currentTime );

Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTime(
        const std::vector< Eigen::Matrix4d >& christoffelSymbols, const Eigen::Vector3d& coordinateVelocity );

Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTime(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState );

Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTimeWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime );

double evaluateProperTimeEquation(
        const Eigen::Matrix4d& metricPerturbation, const Eigen::Vector3d& coordinateVelocity, const int squareRootOrderExpansion = 1 );

double evaluateProperTimeEquation(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const int squareRootOrderExpansion = 1 );

double evaluateProperTimeEquationWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime, const int squareRootOrderExpansion = 1 );

class DirectRelativisticAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    DirectRelativisticAcceleration( const std::shared_ptr< relativity::Metric > spaceTimeMetric,
                                    const std::function< Eigen::Vector6d( ) > acceleratedBodyStateFunction ):
        spaceTimeMetric_( spaceTimeMetric ), acceleratedBodyStateFunction_( acceleratedBodyStateFunction ){ }

    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            this->currentTime_ = currentTime;

            currentAcceleratedBodyState_ =  acceleratedBodyStateFunction_( );
            spaceTimeMetric_->update( currentAcceleratedBodyState_, currentTime, 0, 1 );
            currentAcceleration_ = evaluateRelativisticEquationsOfMotionInCoordinateTime( spaceTimeMetric_, currentAcceleratedBodyState_ );
        }
    }

private:
    std::shared_ptr< relativity::Metric > spaceTimeMetric_;

    std::function< Eigen::Vector6d( ) > acceleratedBodyStateFunction_;

    Eigen::Vector6d currentAcceleratedBodyState_;

    Eigen::Vector3d currentAcceleration_;

};

}

}

#endif // TUDAT_RELATIVISTICEQUATIONSOFMOTION_H
