/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H
#define TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H

#include <map>
#include <stdexcept>
#include <Eigen/Core>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/celestialReferenceFrames.h"
#include "tudat/astro/observation_models/stellarAberrationCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Direction from receiver to transmitter, optionally aberrated, rotated from the global frame to the frame selected by ancillary settings
inline Eigen::Vector3d getAngularPositionDirectionInSelectedFrame(
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& receiverVelocity,
        const double receptionTime,
        const Eigen::Matrix3d& j2000ToGlobalFrameTransformation,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    const Eigen::Vector3d direction = getPositionAngleUseStellarAberration( ancillarySettings )
            ? calculateApparentDirectionWithStellarAberration( relativePosition, receiverVelocity )
            : relativePosition;
    return getPositionAngleReferenceFrameRotationFromJ2000( receptionTime, ancillarySettings ) *
            j2000ToGlobalFrameTransformation.transpose( ) * direction;
}

//! Class for simulating angular position (right ascension/declination) observables.
/*!
 *  Class for simulating angular position (right ascension/declination), using light-time (with light-time corrections)
 *  to determine the states of the link ends (source and receiver).
 *  Without position-angle ancillary settings (position_angle_reference_frame and/or position_angle_use_stellar_aberration), the
 *  angles are computed w.r.t. the axes of the global frame orientation (so, for ECLIPJ2000, they are ecliptic longitude and
 *  latitude). With these settings, the angles are right ascension and declination w.r.t. the selected celestial frame (J2000 by
 *  default), optionally including stellar aberration (which uses the receiver velocity w.r.t. the global frame origin).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class AngularPositionObservationModel : public ObservationModel< 2, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between source and receiver
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     *  \param normalizeRightAscension Whether the right ascension is multiplied by the cosine of the declination
     *  \param j2000ToGlobalFrameTransformation Rotation from ICRF/J2000 to the frame of the link-end states (NaN if not defined)
     *  \param isGlobalFrameOriginBarycentric Whether the global frame origin is the solar system barycentre (a warning is printed
     *  when stellar aberration is used otherwise)
     */
    AngularPositionObservationModel(
            const LinkEnds linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const std::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = nullptr,
            const bool normalizeRightAscension = false,
            const Eigen::Matrix3d& j2000ToGlobalFrameTransformation = Eigen::Matrix3d::Identity( ),
            const bool isGlobalFrameOriginBarycentric = true ):
        ObservationModel< 2, ObservationScalarType, TimeType >(
                angular_position,
                linkEnds,
                observationBiasCalculator,
                std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{ std::make_shared<
                        FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                        std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                                lightTimeCalculator },
                        std::make_shared< LightTimeConvergenceCriteria >( ),
                        false ) } ),
        normalizeRightAscension_( normalizeRightAscension ), j2000ToGlobalFrameTransformation_( j2000ToGlobalFrameTransformation ),
        isGlobalFrameOriginBarycentric_( isGlobalFrameOriginBarycentric )
    {}

    //! Destructor
    ~AngularPositionObservationModel( ) {}

    //! Function to compute ideal angular position observation at given time.
    /*!
     *  This function compute ideal angular position observation at a given time. The time argument can be either the
     *  reception or transmission time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated angular position observable values.
     */
    Eigen::Matrix< ObservationScalarType, 2, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetingsInput = nullptr ) override
    {
        std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetings;
        std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator =
                getLightTimeCalculator( );
        this->setFrequencyProperties( time, linkEndAssociatedWithTime, lightTimeCalculator, ancillarySetingsInput, ancillarySetings );

        this->getFullLinkLightTimeCalculatorFromBase( )->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySetings );

        Eigen::Matrix< ObservationScalarType, 3, 1 > relativePosition =
                linkEndStates.at( 0 ).template cast< ObservationScalarType >( ).segment( 0, 3 ) -
                linkEndStates.at( 1 ).template cast< ObservationScalarType >( ).segment( 0, 3 );
        // Rotate to the celestial frame (and apply stellar aberration) selected by the position-angle ancillary settings, if any
        if( ancillarySetingsInput != nullptr &&
            ( ancillarySetingsInput->hasAncillaryIntData( position_angle_reference_frame ) ||
              ancillarySetingsInput->hasAncillaryIntData( position_angle_use_stellar_aberration ) ) )
        {
            if( !j2000ToGlobalFrameTransformation_.allFinite( ) )
            {
                throw std::runtime_error(
                        "Error, angular position observations with position-angle ancillary settings require a J2000, ICRF or "
                        "ECLIPJ2000 global frame orientation." );
            }
            if( getPositionAngleUseStellarAberration( ancillarySetingsInput ) )
            {
                warnIfStellarAberrationOriginIsNotBarycentric( isGlobalFrameOriginBarycentric_ );
            }
            relativePosition = getAngularPositionDirectionInSelectedFrame( relativePosition.template cast< double >( ),
                                                                           linkEndStates.at( 1 ).segment( 3, 3 ),
                                                                           linkEndTimes.at( 1 ),
                                                                           j2000ToGlobalFrameTransformation_,
                                                                           ancillarySetingsInput )
                                       .template cast< ObservationScalarType >( );
        }

        // Return observable
        double rightAscension = 2.0 *
                std::atan( relativePosition[ 1 ] /
                           ( std::sqrt( relativePosition[ 0 ] * relativePosition[ 0 ] + relativePosition[ 1 ] * relativePosition[ 1 ] ) +
                             relativePosition[ 0 ] ) );
        double declination = mathematical_constants::PI / 2.0 - std::acos( relativePosition[ 2 ] / relativePosition.norm( ) );
        //        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << sphericalRelativeCoordinates.z( ),
        //                 mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ) ).finished( );
        if( !normalizeRightAscension_ )
        {
            return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << rightAscension, declination ).finished( );
        }
        else
        {
            return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << rightAscension * std::cos( declination ), declination ).finished( );
        }
    }

    //! Function to get the object to calculate light time.
    /*!
     * Function to get the object to calculate light time.
     * \return Object to calculate light time.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return this->getSingleLegLightTimeCalculator( );
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( ) } } };
    }

    bool getNormalizeRightAscension( ) const
    {
        return normalizeRightAscension_;
    }

    const Eigen::Matrix3d& getJ2000ToGlobalFrameTransformation( ) const
    {
        return j2000ToGlobalFrameTransformation_;
    }

    ResidualWrappingSettings getResidualWrappingSettings( ) const override
    {
        ResidualWrappingSettings residualWrappingSettings;
        residualWrappingSettings.normalizeRightAscension = normalizeRightAscension_;
        return residualWrappingSettings;
    }

private:
    bool normalizeRightAscension_;

    //! Rotation from ICRF/J2000 to the frame of the link-end states, used when a celestial frame is selected
    Eigen::Matrix3d j2000ToGlobalFrameTransformation_;

    //! Whether the global frame origin is the solar system barycentre
    bool isGlobalFrameOriginBarycentric_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H
