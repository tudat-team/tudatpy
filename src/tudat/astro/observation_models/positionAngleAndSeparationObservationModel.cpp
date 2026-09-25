/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

namespace tudat
{
namespace observation_models
{
namespace
{

double getPositionAngleReferenceEpoch( const double observationTime,
                                       const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    double referenceEpoch = observationTime;
    if( ancillarySettings != nullptr )
    {
        const double ancillaryReferenceEpoch = ancillarySettings->getAncillaryDoubleData( position_angle_reference_epoch, false );
        if( !std::isnan( ancillaryReferenceEpoch ) )
        {
            referenceEpoch = ancillaryReferenceEpoch;
        }
    }

    if( !std::isfinite( referenceEpoch ) )
    {
        throw std::runtime_error( "Position-angle reference epoch must be finite." );
    }
    return referenceEpoch;
}

Eigen::Vector3d getMeanOfDatePoleInJ2000( const double terrestrialTime, const bool useIau2006Precession )
{
    const sofa_interface::PrecessionNutationModel model = useIau2006Precession ? sofa_interface::PrecessionNutationModel::iau_2006_2000a
                                                                               : sofa_interface::PrecessionNutationModel::iau_1976_1980;
    return sofa_interface::getPrecessionMatrix( terrestrialTime, model ).transpose( ) * Eigen::Vector3d::UnitZ( );
}

Eigen::Vector3d getTrueOfDatePoleInJ2000( const double terrestrialTime, const bool useIau2006PrecessionNutation )
{
    const sofa_interface::PrecessionNutationModel model = useIau2006PrecessionNutation
            ? sofa_interface::PrecessionNutationModel::iau_2006_2000a
            : sofa_interface::PrecessionNutationModel::iau_1976_1980;
    return sofa_interface::getPrecessionNutationMatrix( terrestrialTime, model ).transpose( ) * Eigen::Vector3d::UnitZ( );
}

}  // namespace

Eigen::Vector3d getPositionAngleReferencePoleInJ2000( const double observationTime,
                                                      const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    PositionAngleReferenceFrame referenceFrame = j2000_position_angle_reference_frame;
    if( ancillarySettings != nullptr && ancillarySettings->hasAncillaryIntData( position_angle_reference_frame ) )
    {
        const int referenceFrameValue = ancillarySettings->getAncillaryIntData( position_angle_reference_frame );
        if( referenceFrameValue < j2000_position_angle_reference_frame || referenceFrameValue > custom_position_angle_reference_pole )
        {
            throw std::runtime_error( "Position-angle reference-frame ancillary setting contains an unsupported enum value: " +
                                      std::to_string( referenceFrameValue ) + "." );
        }
        referenceFrame = static_cast< PositionAngleReferenceFrame >( referenceFrameValue );
    }

    switch( referenceFrame )
    {
        case j2000_position_angle_reference_frame:
            return Eigen::Vector3d::UnitZ( );
        case b1950_position_angle_reference_frame:
            // NAIF's built-in B1950-to-J2000 frame rotation, applied to the B1950 positive z-axis.
            return Eigen::Vector3d( -4.8590038153592703e-3, -2.7162594714247041e-5, 9.9998819460237420e-1 );
        case mean_of_date_iau_1976_position_angle_reference_frame:
        case true_of_date_iau_1976_1980_position_angle_reference_frame:
        case mean_of_date_iau_2006_position_angle_reference_frame:
        case true_of_date_iau_2006_2000a_position_angle_reference_frame: {
            const double referenceEpoch = getPositionAngleReferenceEpoch( observationTime, ancillarySettings );
            const double terrestrialTime = sofa_interface::convertTDBtoTT( referenceEpoch, Eigen::Vector3d::Zero( ) );
            if( referenceFrame == mean_of_date_iau_1976_position_angle_reference_frame )
            {
                return getMeanOfDatePoleInJ2000( terrestrialTime, false );
            }
            if( referenceFrame == true_of_date_iau_1976_1980_position_angle_reference_frame )
            {
                return getTrueOfDatePoleInJ2000( terrestrialTime, false );
            }
            if( referenceFrame == mean_of_date_iau_2006_position_angle_reference_frame )
            {
                return getMeanOfDatePoleInJ2000( terrestrialTime, true );
            }
            return getTrueOfDatePoleInJ2000( terrestrialTime, true );
        }
        case custom_position_angle_reference_pole: {
            if( ancillarySettings == nullptr )
            {
                throw std::runtime_error( "Custom position-angle reference pole selected without ancillary settings." );
            }
            const std::vector< double > referencePoleValues =
                    ancillarySettings->getAncillaryDoubleVectorData( position_angle_custom_reference_pole );
            if( referencePoleValues.size( ) != 3 )
            {
                throw std::runtime_error( "Custom position-angle reference pole must contain exactly three values." );
            }
            const Eigen::Vector3d referencePole( referencePoleValues.at( 0 ), referencePoleValues.at( 1 ), referencePoleValues.at( 2 ) );
            if( !referencePole.allFinite( ) || referencePole.norm( ) <= 100.0 * std::numeric_limits< double >::epsilon( ) )
            {
                throw std::runtime_error( "Custom position-angle reference pole must be finite and non-zero." );
            }
            return referencePole.normalized( );
        }
        default:
            throw std::runtime_error( "Position-angle reference-frame ancillary setting contains an unsupported enum value: " +
                                      std::to_string( static_cast< int >( referenceFrame ) ) + "." );
    }
}

PositionAngleDirectionType getPositionAngleDirectionType(
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    PositionAngleDirectionType directionType = astrometric_position_angle_direction;
    if( ancillarySettings != nullptr && ancillarySettings->hasAncillaryIntData( position_angle_direction_type ) )
    {
        const int directionTypeValue = ancillarySettings->getAncillaryIntData( position_angle_direction_type );
        if( directionTypeValue < astrometric_position_angle_direction || directionTypeValue > aberrated_position_angle_direction )
        {
            throw std::runtime_error( "Position-angle direction-type ancillary setting contains an unsupported enum value: " +
                                      std::to_string( directionTypeValue ) + "." );
        }
        directionType = static_cast< PositionAngleDirectionType >( directionTypeValue );
    }

    switch( directionType )
    {
        case astrometric_position_angle_direction:
        case aberrated_position_angle_direction:
            return directionType;
        default:
            throw std::runtime_error( "Position-angle direction-type ancillary setting contains an unsupported enum value: " +
                                      std::to_string( static_cast< int >( directionType ) ) + "." );
    }
}

}  // namespace observation_models
}  // namespace tudat
