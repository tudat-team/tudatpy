/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PIXEL_COORDINATES_PARTIAL_H
#define TUDAT_PIXEL_COORDINATES_PARTIAL_H

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/system_models/camera.h"
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

class PixelCoordinatesScaling : public DirectPositionPartialScaling< 2 >
{
public:
    PixelCoordinatesScaling( std::function< Eigen::Quaterniond( const double epoch ) > rotationFromInertialToCameraFrameFunction,
                             std::function< Eigen::Matrix< double, 2, 3 >( const Eigen::Vector3d& cameraFramePosition, const double time ) >
                                     pixelLinePartialWrtCameraFramePositionFunction,
                             const bool correctForStellarAberration = false ):
        DirectPositionPartialScaling< 2 >( observation_models::pixel_coordinates ),
        rotationFromInertialToCameraFrameFunction_( rotationFromInertialToCameraFrameFunction ),
        pixelLinePartialWrtCameraFramePositionFunction_( pixelLinePartialWrtCameraFramePositionFunction ),
        correctForStellarAberration_( correctForStellarAberration )
    {
        observerIndex_ = observation_models::getSingleLinkStateEntryIndices( observation_models::pixel_coordinates )
                                 .at( observation_models::receiver );
        observedBodyIndex_ = observation_models::getSingleLinkStateEntryIndices( observation_models::pixel_coordinates )
                                     .at( observation_models::transmitter );
    }

    ~PixelCoordinatesScaling( ) {}

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
                 const Eigen::VectorXd currentObservation,
                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr ) override;

    //! Get the position scaling factor (accounting for light-time correction) for a given link end type
    /*!
     *  Get the position scaling factor (accounting for light-time correction) for a given link end type.
     *  Considers observations along the line of sight coinciding with the camera frame z-axis (camera boresight).
     *  Considers gnomonic projection of observed object onto camera focal plane (source: Owen Jr, 2024).
     *
     *  \param linkEndType Link end type for which scaling factor is to be retrieved.
     *  \return Position scaling factor for given link end type.
     */
    Eigen::Matrix< double, 2, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType ) override
    {
        return realPositionScalingFactor_ * ( ( linkEndType == observation_models::receiver ) ? ( 1.0 ) : ( -1.0 ) );
    }

    //! Get the position scaling factor (not accounting for light-time correction) for a given link end type
    /*!  Get the position scaling factor (not accounting for light-time correction) for a given link end type.
     *  Considers observations along the line of sight coinciding with the camera frame z-axis (camera boresight).
     *  Considers gnomonic projection of observed object onto camera focal plane (source: Owen Jr, 2024).
     *  \param linkEndType Link end type for which scaling factor is to be retrieved.
     *  \return Position scaling factor for given link end type.
     */
    Eigen::Matrix< double, 2, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType ) override
    {
        return positionScalingFactor_ * ( ( linkEndType == observation_models::receiver ) ? ( 1.0 ) : ( -1.0 ) );
    }

    //! Only the receiver velocity contributes, through stellar aberration.
    bool isVelocityScalingNonZero( const observation_models::LinkEndType linkEndType ) override
    {
        return correctForStellarAberration_ && linkEndType == observation_models::receiver;
    }

    Eigen::Matrix< double, 2, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType ) override
    {
        if( linkEndType == observation_models::receiver )
        {
            return velocityScalingFactor_;
        }
        return Eigen::Matrix< double, 2, 3 >::Zero( );
    }

    Eigen::Matrix< double, 2, 3 > getFixedTimeVelocityScalingFactor( const observation_models::LinkEndType linkEndType ) override
    {
        return getVelocityScalingFactor( linkEndType );
    }

    //! Supply the reference-body rotation partials contributing directly to the camera orientation.
    void setRotationMatrixPartials( const std::string& referenceBodyName,
                                    const std::function< Eigen::Quaterniond( const double ) >& rotationFromReferenceFrameToInertialFunction,
                                    const RotationMatrixPartialNamedList& rotationMatrixPartials )
    {
        pointingReferenceBodyName_ = referenceBodyName;
        rotationFromReferenceFrameToInertialFunction_ = rotationFromReferenceFrameToInertialFunction;
        rotationMatrixPartials_ = rotationMatrixPartials;
    }

    //! Whether the camera orientation depends on the given parameter (pointing correction, or reference-body rotation).
    bool useLinkIndependentPartialsForParameter(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier ) override;

    //! Identify the receiver camera, whose pointing corrections are used for the pointing correction partials.
    void setPointingCamera( const observation_models::LinkEndId& cameraLinkEnd, const std::shared_ptr< system_models::Camera >& camera )
    {
        pointingCorrectionParameterIdentifier_ =
                std::make_pair( estimatable_parameters::camera_pointing_correction,
                                std::make_pair( cameraLinkEnd.bodyName_, cameraLinkEnd.getReferencePointName( ) ) );
        pointingCamera_ = camera;
    }

    std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > getLinkIndependentPartials(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ) override;

    //! Get the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!  Get the factor by which the light-time partials should be scaled in one-way observation partial.
     *  \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    Eigen::Vector2d getLightTimePartialScalingFactor( ) override
    {
        return lightTimeCorrectionScalingFactor_;
    }

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( ) override
    {
        return currentLinkEndType_;
    }

private:
    Eigen::Matrix< double, 2, 3 > positionScalingFactor_;
    Eigen::Matrix< double, 2, 3 > realPositionScalingFactor_;
    Eigen::Matrix< double, 2, 3 > velocityScalingFactor_ = Eigen::Matrix< double, 2, 3 >::Zero( );
    Eigen::Vector2d lightTimeCorrectionScalingFactor_;
    //! Apparent unit direction when aberration is enabled; otherwise the geometric range vector.
    Eigen::Vector3d currentRelativeRangeVectorCameraFrame_;
    Eigen::Matrix< double, 2, 3 > currentPixelLinePartialWrtCameraFramePosition_;
    Eigen::Vector3d currentInertialDirectionToProject_;
    Eigen::Matrix3d currentRotationFromInertialToCameraFrame_;
    double currentReceiverTime_ = TUDAT_NAN;
    std::function< Eigen::Quaterniond( const double epoch ) > rotationFromInertialToCameraFrameFunction_;
    std::function< Eigen::Matrix< double, 2, 3 >( const Eigen::Vector3d& cameraFramePosition, const double time ) >
            pixelLinePartialWrtCameraFramePositionFunction_;
    std::string pointingReferenceBodyName_;
    std::function< Eigen::Quaterniond( const double ) > rotationFromReferenceFrameToInertialFunction_;
    RotationMatrixPartialNamedList rotationMatrixPartials_;
    estimatable_parameters::EstimatebleParameterIdentifier pointingCorrectionParameterIdentifier_;
    std::shared_ptr< system_models::Camera > pointingCamera_;
    observation_models::LinkEndType currentLinkEndType_;
    int observerIndex_;
    int observedBodyIndex_;
    bool hasCurrentCameraFrameQuantities_ = false;
    bool correctForStellarAberration_;
};
}  // namespace observation_partials
}  // namespace tudat
#endif  // TUDAT_PIXEL_COORDINATES_PARTIAL_H
