/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CAMERA_PIXELS_OBSERVATION_MODEL_H
#define TUDAT_CAMERA_PIXELS_OBSERVATION_MODEL_H

#include <map>
#include <functional>
#include <Eigen/Core>

#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/system_models/camera.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating camera pixel observables.
/*!
 *  Class for simulating camera pixel observables, using light-time (with light-time corrections)
 *  to determine the states of the link ends (source and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 *  The mapping from relative inertial position to pixel coordinates can be obtained from a tudat::system_models::Camera object, or be
 * defined as a custom function.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class CameraPixelsObservationModel : public ObservationModel< 2, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between source and receiver
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    CameraPixelsObservationModel(
            const LinkEnds linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const std::shared_ptr< system_models::Camera > camera,
            const std::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = nullptr ):
        ObservationModel< 2, ObservationScalarType, TimeType >( camera_pixels, linkEnds, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator ), camera_( camera )
    {}

    //! Destructor
    ~CameraPixelsObservationModel( ) {}

    //! Function to compute ideal camera pixel observation at given time.
    /*!
     * Function to compute ideal camera pixel observation at given time. The times and states of the link ends are also returned in full
     * precision (determined by class template arguments). These states and times are returned by reference.
     * \param time Time at which observable is to be evaluated.
     * \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time is kept constant (to
     * input value) \param linkEndTimes List of times at each link end during observation (returned by reference). \param linkEndStates List
     * of states at each link end during observation (returned by reference). \return Ideal camera pixel observable at given time.
     */
    Eigen::Matrix< ObservationScalarType, 2, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetingsInput = nullptr ) override
    {
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error in CameraPixelsObservationModel, time must be associated with receiver link end." );
        }

        std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetings;
        this->setFrequencyProperties( time, linkEndAssociatedWithTime, lightTimeCalculator_, ancillarySetingsInput, ancillarySetings );
        const bool isTimeAtObserver = true;

        if( ancillarySetings != nullptr )
        {
            throw std::runtime_error( "Error in CameraPixelsObservationModel, ancillary settings are not supported." );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > sourceState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > observerState;

        // Compute light-time and observer and observed body states.
        // States are modified by light-time calculator.
        ObservationScalarType lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                observerState, sourceState, time, isTimeAtObserver, ancillarySetings );

        Eigen::Matrix< ObservationScalarType, 3, 1 > inertialRelativePosition = sourceState.segment( 0, 3 ) - observerState.segment( 0, 3 );

        linkEndTimes.clear( );
        linkEndStates.clear( );

        // Link-end order for camera_pixels is [transmitter, receiver].
        linkEndStates.push_back( sourceState.template cast< double >( ) );
        linkEndStates.push_back( observerState.template cast< double >( ) );

        linkEndTimes.push_back( static_cast< double >( time - lightTime ) );
        linkEndTimes.push_back( static_cast< double >( time ) );

        Eigen::Vector2d pixelCoordinates = camera_->calculateObservableFromInertial( inertialRelativePosition.template cast< double >( ),
                                                                                     static_cast< double >( time ) );
        return pixelCoordinates.template cast< ObservationScalarType >( );
    }

    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

protected:
    //! Object to compute the light-time (including any corrections w.r.t. Euclidean case) between source and receiver
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator_;

    //! Camera object to map relative position of source and receiver to camera pixel coordinates.
    std::shared_ptr< system_models::Camera > camera_;
};
}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CAMERA_PIXELS_OBSERVATION_MODEL_H
