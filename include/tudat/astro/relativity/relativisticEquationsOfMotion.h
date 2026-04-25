/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

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

//! Evaluate four-acceleration from Christoffel symbols and four-velocity.
/*!
 *  \param christoffelSymbols Christoffel symbols \f$\Gamma^\mu_{\alpha\beta}\f$ at evaluation point.
 *  \param fourVelocity Four-velocity \f$u^\mu\f$.
 *  \return Four-acceleration \f$du^\mu/d\tau\f$.
 */
Eigen::Vector4d evaluateFourAcceleration(
        const std::vector< Eigen::Matrix4d >& christoffelSymbols, const Eigen::Vector4d& fourVelocity );

//! Evaluate four-acceleration using a metric object and current 8D state.
/*!
 *  \param spaceTimeMetric Metric object used in four-acceleration evaluation.
 *  \param currentState Current 8D state used in metric/four-velocity evaluation.
 *  \return Four-acceleration \f$du^\mu/d\tau\f$.
 */
Eigen::Vector4d evaluateFourAcceleration(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric,
        const Eigen::Matrix< double, 8, 1 >& currentState );

//! Update metric and evaluate four-acceleration at a given state and epoch.
/*!
 *  \param spaceTimeMetric Metric object to update.
 *  \param currentState Current 8D state used in metric/four-velocity evaluation.
 *  \param currentTime Evaluation epoch.
 *  \return Four-acceleration \f$du^\mu/d\tau\f$.
 */
Eigen::Vector4d evaluateFourAccelerationWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 8, 1 >& currentState,
        const double currentTime );

//! Evaluate coordinate-time equations of motion from Christoffel symbols.
/*!
 *  \param christoffelSymbols Christoffel symbols \f$\Gamma^\mu_{\alpha\beta}\f$ at evaluation point.
 *  \param coordinateVelocity Coordinate velocity \f$d\mathbf{x}/dt\f$.
 *  \return Coordinate acceleration \f$d^2\mathbf{x}/dt^2\f$.
 */
Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTime(
        const std::vector< Eigen::Matrix4d >& christoffelSymbols, const Eigen::Vector3d& coordinateVelocity );

//! Evaluate coordinate-time equations of motion using current values in a metric object.
/*!
 *  \param spaceTimeMetric Metric object with already-updated state and metric quantities.
 *  \param currentState Current translational state.
 *  \return Coordinate acceleration \f$d^2\mathbf{x}/dt^2\f$.
 */
Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTime(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState );

//! Update metric and evaluate coordinate-time equations of motion.
/*!
 *  \param spaceTimeMetric Metric object to update.
 *  \param currentState Current translational state.
 *  \param currentTime Evaluation epoch.
 *  \return Coordinate acceleration \f$d^2\mathbf{x}/dt^2\f$.
 */
Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTimeWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime );

//! Evaluate proper-time rate equation from metric perturbation and coordinate velocity.
/*!
 *  \param metricPerturbation Metric perturbation tensor at evaluation point.
 *  \param coordinateVelocity Coordinate velocity \f$d\mathbf{x}/dt\f$.
 *  \param squareRootOrderExpansion Expansion order used in the square-root approximation.
 *  \return Proper-time rate \f$d\tau/dt - 1\f$ contribution used in propagation.
 */
double evaluateProperTimeEquation(
        const Eigen::Matrix4d& metricPerturbation, const Eigen::Vector3d& coordinateVelocity, const int squareRootOrderExpansion = 1 );

//! Evaluate proper-time rate equation using current values in a metric object.
/*!
 *  \param spaceTimeMetric Metric object with already-updated state and metric quantities.
 *  \param currentState Current translational state.
 *  \param squareRootOrderExpansion Expansion order used in the square-root approximation.
 *  \return Proper-time rate \f$d\tau/dt - 1\f$ contribution used in propagation.
 */
double evaluateProperTimeEquation(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const int squareRootOrderExpansion = 1 );

//! Update metric and evaluate proper-time rate equation.
/*!
 *  \param spaceTimeMetric Metric object to update.
 *  \param currentState Current translational state.
 *  \param currentTime Evaluation epoch.
 *  \param squareRootOrderExpansion Expansion order used in the square-root approximation.
 *  \return Proper-time rate \f$d\tau/dt - 1\f$ contribution used in propagation.
 */
double evaluateProperTimeEquationWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime, const int squareRootOrderExpansion = 1 );

//! Direct relativistic acceleration model evaluated from a metric.
class DirectRelativisticAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    //! Constructor.
    /*!
     *  \param spaceTimeMetric Metric used to evaluate relativistic acceleration.
     *  \param acceleratedBodyStateFunction Function returning state of accelerated body.
     */
    DirectRelativisticAcceleration( const std::shared_ptr< relativity::Metric > spaceTimeMetric,
                                    const std::function< Eigen::Vector6d( ) > acceleratedBodyStateFunction ):
        spaceTimeMetric_( spaceTimeMetric ), acceleratedBodyStateFunction_( acceleratedBodyStateFunction ){ }

    std::shared_ptr< relativity::Metric > getSpaceTimeMetric( ) const
    {
        return spaceTimeMetric_;
    }

    //! Update internal acceleration value for a new time.
    /*!
     *  \param currentTime Current epoch.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            this->currentTime_ = currentTime;

            currentAcceleratedBodyState_ =  acceleratedBodyStateFunction_( );
            spaceTimeMetric_->update( currentAcceleratedBodyState_, currentTime, 0, 1 );
            this->currentAcceleration_ =
                    evaluateRelativisticEquationsOfMotionInCoordinateTime( spaceTimeMetric_, currentAcceleratedBodyState_ );

        }
    }

private:
    std::shared_ptr< relativity::Metric > spaceTimeMetric_;

    std::function< Eigen::Vector6d( ) > acceleratedBodyStateFunction_;

    Eigen::Vector6d currentAcceleratedBodyState_;

};

}

}

#endif // TUDAT_RELATIVISTICEQUATIONSOFMOTION_H
