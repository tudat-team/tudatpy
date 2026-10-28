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

PositionAngleReferenceFrame getPositionAngleReferenceFrame( const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    if( ancillarySettings == nullptr || !ancillarySettings->hasAncillaryIntData( position_angle_reference_frame ) )
    {
        return j2000_position_angle_reference_frame;
    }
    const int referenceFrameValue = ancillarySettings->getAncillaryIntData( position_angle_reference_frame );
    if( referenceFrameValue < j2000_position_angle_reference_frame || referenceFrameValue > custom_position_angle_reference_pole )
    {
        throw std::runtime_error( "Position-angle reference-frame ancillary setting contains an unsupported enum value: " +
                                  std::to_string( referenceFrameValue ) + "." );
    }
    return static_cast< PositionAngleReferenceFrame >( referenceFrameValue );
}

}  // namespace

Eigen::Matrix3d getPositionAngleReferenceFrameRotationFromJ2000(
        const double observationTime,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    const PositionAngleReferenceFrame referenceFrame = getPositionAngleReferenceFrame( ancillarySettings );
    switch( referenceFrame )
    {
        case j2000_position_angle_reference_frame:
            return Eigen::Matrix3d::Identity( );
        case b1950_position_angle_reference_frame:
            // Transpose of NAIF's built-in B1950-to-J2000 frame rotation.
            return ( Eigen::Matrix3d( ) << 9.9992567949568767e-1,
                     -1.1181483220466210e-2,
                     -4.8590038153592703e-3,
                     1.1181483239171700e-2,
                     9.9993748489331347e-1,
                     -2.7162594714247041e-5,
                     4.8590037723143067e-3,
                     -2.7170293744021070e-5,
                     9.9998819460237420e-1 )
                    .finished( )
                    .transpose( );
        case mean_of_date_iau_1976_position_angle_reference_frame:
        case true_of_date_iau_1976_1980_position_angle_reference_frame:
        case mean_of_date_iau_2006_position_angle_reference_frame:
        case true_of_date_iau_2006_2000a_position_angle_reference_frame: {
            const double referenceEpoch = getPositionAngleReferenceEpoch( observationTime, ancillarySettings );
            const double terrestrialTime = sofa_interface::convertTDBtoTT( referenceEpoch, Eigen::Vector3d::Zero( ) );
            const bool useIau2006 = ( referenceFrame == mean_of_date_iau_2006_position_angle_reference_frame ||
                                      referenceFrame == true_of_date_iau_2006_2000a_position_angle_reference_frame );
            if( referenceFrame == mean_of_date_iau_1976_position_angle_reference_frame ||
                referenceFrame == mean_of_date_iau_2006_position_angle_reference_frame )
            {
                return sofa_interface::getPrecessionMatrix( terrestrialTime, basic_astrodynamics::JULIAN_DAY_ON_J2000, useIau2006 );
            }
            return sofa_interface::getPrecessionNutationMatrix( terrestrialTime, basic_astrodynamics::JULIAN_DAY_ON_J2000, useIau2006 );
        }
        default:
            throw std::runtime_error( "Error, a full rotation is not defined for position-angle reference frame setting " +
                                      std::to_string( static_cast< int >( referenceFrame ) ) + "." );
    }
}

Eigen::Vector3d getPositionAngleReferencePoleInJ2000( const double observationTime,
                                                      const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    if( getPositionAngleReferenceFrame( ancillarySettings ) != custom_position_angle_reference_pole )
    {
        return getPositionAngleReferenceFrameRotationFromJ2000( observationTime, ancillarySettings ).transpose( ) *
                Eigen::Vector3d::UnitZ( );
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

bool getPositionAngleUseStellarAberration( const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    if( ancillarySettings != nullptr && ancillarySettings->hasAncillaryIntData( position_angle_use_stellar_aberration ) )
    {
        return ancillarySettings->getAncillaryIntData( position_angle_use_stellar_aberration ) != 0;
    }
    return false;
}

}  // namespace observation_models
}  // namespace tudat
