/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CAMERA_PIXELS_PARTIAL_H
#define TUDAT_CAMERA_PIXELS_PARTIAL_H

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include <functional>
#include <Eigen/Core>

namespace tudat
{
namespace observation_partials
{

Eigen::Matrix< double, 1, 3 > calculatePartialOfPixelUWrtObserverPositionCameraFrame(
        Eigen::Vector3d relativeRangeVectorCameraFrame,
        Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix );

Eigen::Matrix< double, 1, 3 > calculatePartialOfPixelVWrtObserverPositionCameraFrame(
        Eigen::Vector3d relativeRangeVectorCameraFrame,
        Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix );

Eigen::Matrix< double, 2, 3 > calculatePartialOfPixelsWrtLinkEndPositionCameraFrame( Eigen::Vector3d relativeRangeVectorCameraFrame,
                                                                                     Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix,
                                                                                     const bool isLinkEndObserver );

class CameraPixelsScaling : public DirectPositionPartialScaling< 2 >
{
public:
    CameraPixelsScaling( std::function< Eigen::Quaterniond( const double epoch ) > rotationFromInertialToCameraFrameFunction,
                         Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix ):
        DirectPositionPartialScaling< 2 >( observation_models::camera_pixels ),
        rotationFromInertialToCameraFrameFunction_( rotationFromInertialToCameraFrameFunction ), focalLengthsMatrix_( focalLengthsMatrix )
    {
        observerIndex_ =
                observation_models::getSingleLinkStateEntryIndices( observation_models::camera_pixels ).at( observation_models::observer );
        observedBodyIndex_ = observation_models::getSingleLinkStateEntryIndices( observation_models::camera_pixels )
                                     .at( observation_models::observed_body );
    }

    ~CameraPixelsScaling( ) {}

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states.
     *  Considers observations along the line of sight coinciding with the camera frame z-axis (camera boresight).
     *  Considers gnomonic projection of observed object onto camera focal plane (source: Owen Jr, 2024).
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    //! Get the position scaling factor (accounting for light-time correction) for a given link end type
    /*!
     *  Get the position scaling factor (accounting for light-time correction) for a given link end type.
     *  Considers observations along the line of sight coinciding with the camera frame z-axis (camera boresight).
     *  Considers gnomonic projection of observed object onto camera focal plane (source: Owen Jr, 2024).
     *
     *  \param linkEndType Link end type for which scaling factor is to be retrieved.
     *  \return Position scaling factor for given link end type.
     */
    Eigen::Matrix< double, 2, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return realPositionScalingFactor_ * ( ( linkEndType == observation_models::observer ) ? ( 1.0 ) : ( -1.0 ) );
    }

    //! Get the position scaling factor (not accounting for light-time correction) for a given link end type
    /*!  Get the position scaling factor (not accounting for light-time correction) for a given link end type.
     *  Considers observations along the line of sight coinciding with the camera frame z-axis (camera boresight).
     *  Considers gnomonic projection of observed object onto camera focal plane (source: Owen Jr, 2024).
     *  \param linkEndType Link end type for which scaling factor is to be retrieved.
     *  \return Position scaling factor for given link end type.
     */
    Eigen::Matrix< double, 2, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return positionScalingFactor_ * ( ( linkEndType == observation_models::observer ) ? ( 1.0 ) : ( -1.0 ) );
    }

    //! Get the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!  Get the factor by which the light-time partials should be scaled in one-way observation partial.
     *  \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    Eigen::Vector2d getLightTimePartialScalingFactor( )
    {
        return lightTimeCorrectionScalingFactor_;
    }

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:
    Eigen::Matrix< double, 2, 3 > positionScalingFactor_;
    Eigen::Matrix< double, 2, 3 > realPositionScalingFactor_;
    Eigen::Vector2d lightTimeCorrectionScalingFactor_;
    std::function< Eigen::Quaterniond( const double epoch ) > rotationFromInertialToCameraFrameFunction_;
    Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix_;
    observation_models::LinkEndType currentLinkEndType_;
    int observerIndex_;
    int observedBodyIndex_;
};
}  // namespace observation_partials
}  // namespace tudat
#endif  // TUDAT_CAMERA_PIXELS_PARTIAL_H
