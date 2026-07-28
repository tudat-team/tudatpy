/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createSp3Ephemeris.h"

#include <algorithm>
#include <cctype>
#include <set>

#include "tudat/astro/basic_astro/physicalConstants.h"
#if TUDAT_BUILD_WITH_SOFA_INTERFACE
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#endif
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/io/readSp3File.h"

namespace tudat
{

namespace simulation_setup
{

enum class Sp3FrameFamily { terrestrial, gcrs, j2000, eclipj2000, unknown };

//! Normalize an SP3 coordinate-system or time-system tag for comparison.
static std::string normalizeSp3FrameName( std::string frameName )
{
    frameName.erase(
            std::remove_if(
                    frameName.begin( ), frameName.end( ), []( const unsigned char character ) { return std::isspace( character ); } ),
            frameName.end( ) );
    std::transform( frameName.begin( ), frameName.end( ), frameName.begin( ), []( const unsigned char character ) {
        return static_cast< char >( std::toupper( character ) );
    } );
    return frameName;
}

//! Convert a two-digit ITRF year to the naming convention used by Tudat's ITRF transformations.
static std::string getCanonicalItrfNameFromTwoDigitYear( const std::string& year )
{
    const int numericalYear = std::stoi( year );
    return numericalYear >= 88 ? "ITRF" + year : "ITRF20" + year;
}

//! Resolve an SP3 coordinate-system tag to a supported frame family and canonical Tudat name.
static std::pair< Sp3FrameFamily, std::string > identifySp3Frame( const std::string& rawFrameName )
{
    const std::string frameName = normalizeSp3FrameName( rawFrameName );
    if( frameName == "GCRF" || frameName == "GCRS" )
    {
        return { Sp3FrameFamily::gcrs, "GCRS" };
    }
    if( frameName == "EME00" || frameName == "EME2K" || frameName == "EME2000" || frameName == "J2000" )
    {
        return { Sp3FrameFamily::j2000, "J2000" };
    }
    if( frameName == "ECLIPJ2000" )
    {
        return { Sp3FrameFamily::eclipj2000, "ECLIPJ2000" };
    }
    if( frameName == "ITRS" || frameName == "ECEF" )
    {
        return { Sp3FrameFamily::terrestrial, "ITRS" };
    }
    if( frameName == "WGS84" )
    {
        return { Sp3FrameFamily::terrestrial, "WGS84" };
    }
    if( frameName.rfind( "ITRF", 0 ) == 0 )
    {
        const std::string year = frameName.substr( 4 );
        if( std::all_of( year.begin( ), year.end( ), []( const unsigned char character ) { return std::isdigit( character ); } ) )
        {
            if( year.size( ) == 2 )
            {
                return { Sp3FrameFamily::terrestrial, getCanonicalItrfNameFromTwoDigitYear( year ) };
            }
            if( year.size( ) == 4 && ( year.rfind( "19", 0 ) == 0 || year.rfind( "20", 0 ) == 0 ) )
            {
                return { Sp3FrameFamily::terrestrial, year.rfind( "19", 0 ) == 0 ? "ITRF" + year.substr( 2, 2 ) : "ITRF" + year };
            }
        }
    }
    if( frameName.size( ) == 5 && std::isdigit( static_cast< unsigned char >( frameName.at( 3 ) ) ) &&
        std::isdigit( static_cast< unsigned char >( frameName.at( 4 ) ) ) )
    {
        const std::string prefix = frameName.substr( 0, 3 );
        if( prefix == "IER" || prefix == "ITR" || prefix == "IGS" || prefix == "IGB" || prefix == "IGC" || prefix == "SLR" )
        {
            return { Sp3FrameFamily::terrestrial, getCanonicalItrfNameFromTwoDigitYear( frameName.substr( 3, 2 ) ) };
        }
    }
    return { Sp3FrameFamily::unknown, frameName };
}

//! Check whether Tudat has Helmert parameters for an ITRF realization.
static bool isSupportedItrfHelmertFrame( const std::string& frameName )
{
    static const std::set< std::string > supportedFrames = { "ITRF88", "ITRF89", "ITRF90",   "ITRF91",   "ITRF92",   "ITRF93",  "ITRF94",
                                                             "ITRF96", "ITRF97", "ITRF2000", "ITRF2005", "ITRF2008", "ITRF2014" };
    return supportedFrames.count( frameName ) > 0;
}

#if TUDAT_BUILD_WITH_SOFA_INTERFACE
//! Map an SP3 time-system tag to the time scale and epoch offset required by Tudat's Earth-rotation model.
static std::pair< basic_astrodynamics::TimeScales, double > getSp3EarthRotationTimeScaleAndOffset( const std::string& timeScale )
{
    const std::string normalizedTimeScale = normalizeSp3FrameName( timeScale );
    if( normalizedTimeScale == "TAI" )
    {
        return { basic_astrodynamics::tai_scale, 0.0 };
    }
    if( normalizedTimeScale == "GPS" || normalizedTimeScale == "GAL" || normalizedTimeScale == "QZS" || normalizedTimeScale == "IRN" )
    {
        return { basic_astrodynamics::tai_scale, 19.0 };
    }
    if( normalizedTimeScale == "BDT" )
    {
        return { basic_astrodynamics::tai_scale, 33.0 };
    }
    if( normalizedTimeScale == "UTC" )
    {
        return { basic_astrodynamics::utc_scale, 0.0 };
    }
    if( normalizedTimeScale == "GLO" )
    {
        return { basic_astrodynamics::utc_scale, -3.0 * 3600.0 };
    }
    throw std::invalid_argument( "Cannot rotate SP3 states expressed in unsupported time system '" + timeScale + "'." );
}

//! Return the base-frame name expected by Tudat's ITRS rotation model for an SP3 inertial frame.
static std::string getSp3EarthRotationBaseFrame( const Sp3FrameFamily family )
{
    switch( family )
    {
        case Sp3FrameFamily::gcrs:
            return "GCRS";
        case Sp3FrameFamily::j2000:
            return "J2000";
        case Sp3FrameFamily::eclipj2000:
            return "ECLIPJ2000";
        default:
            throw std::invalid_argument( "The requested SP3 frame is not a supported inertial frame." );
    }
}
#endif

//! Retrieve one SP3 state history and convert it to the requested environment frame.
static std::map< double, Eigen::Vector6d > getTransformedSp3StateHistory( const input_output::Sp3FileContents& fileContents,
                                                                          const std::string& satelliteId,
                                                                          const std::string& targetFrameOrientation )
{
    const std::map< double, Eigen::Vector6d > sourceStateHistory = fileContents.getSatelliteStateHistory( satelliteId );
    if( fileContents.frameName.empty( ) )
    {
        throw std::invalid_argument( "Cannot create SP3 ephemeris settings because the file does not declare a reference frame." );
    }
    if( normalizeSp3FrameName( targetFrameOrientation ) == normalizeSp3FrameName( fileContents.frameName ) )
    {
        return sourceStateHistory;
    }

    const std::pair< Sp3FrameFamily, std::string > sourceFrame = identifySp3Frame( fileContents.frameName );
    const std::pair< Sp3FrameFamily, std::string > targetFrame = identifySp3Frame( targetFrameOrientation );
    if( sourceFrame.first == Sp3FrameFamily::unknown || targetFrame.first == Sp3FrameFamily::unknown )
    {
        throw std::invalid_argument( "Cannot transform SP3 states from frame '" + fileContents.frameName + "' to frame '" +
                                     targetFrameOrientation + "': one or both frame tags are not recognized." );
    }
    if( sourceFrame == targetFrame )
    {
        return sourceStateHistory;
    }

    const bool sourceIsSupportedTerrestrial = sourceFrame.first != Sp3FrameFamily::terrestrial || sourceFrame.second == "ITRS" ||
            isSupportedItrfHelmertFrame( sourceFrame.second );
    const bool targetIsSupportedTerrestrial = targetFrame.first != Sp3FrameFamily::terrestrial || targetFrame.second == "ITRS" ||
            isSupportedItrfHelmertFrame( targetFrame.second );
    if( !sourceIsSupportedTerrestrial || !targetIsSupportedTerrestrial )
    {
        throw std::invalid_argument( "Cannot transform SP3 states from frame '" + fileContents.frameName + "' to frame '" +
                                     targetFrameOrientation + "': no unambiguous Tudat frame-realization transformation is available." );
    }

    const double referenceEpochOffset =
            ( fileContents.referenceJulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;
    std::map< double, Eigen::Vector6d > transformedStateHistory;

    if( sourceFrame.first == Sp3FrameFamily::terrestrial && targetFrame.first == Sp3FrameFamily::terrestrial )
    {
        if( sourceFrame.second == "ITRS" || targetFrame.second == "ITRS" )
        {
            throw std::invalid_argument( "Cannot transform between generic ITRS and a specific ITRF realization." );
        }
        for( const auto& state : sourceStateHistory )
        {
            transformedStateHistory[ state.first ] = reference_frames::convertStateBetweenItrfFrames(
                    state.second, state.first + referenceEpochOffset, sourceFrame.second, targetFrame.second );
        }
        return transformedStateHistory;
    }

#if TUDAT_BUILD_WITH_SOFA_INTERFACE
    const std::pair< basic_astrodynamics::TimeScales, double > earthRotationTime =
            getSp3EarthRotationTimeScaleAndOffset( fileContents.timeScale );
    const std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            earth_orientation::createStandardEarthOrientationCalculator( );
    std::shared_ptr< ephemerides::RotationalEphemeris > sourceRotationModel;
    std::shared_ptr< ephemerides::RotationalEphemeris > targetRotationModel;
    if( sourceFrame.first != Sp3FrameFamily::terrestrial )
    {
        sourceRotationModel = std::make_shared< ephemerides::GcrsToItrsRotationModel >(
                earthOrientationCalculator, earthRotationTime.first, getSp3EarthRotationBaseFrame( sourceFrame.first ) );
    }
    if( targetFrame.first != Sp3FrameFamily::terrestrial )
    {
        targetRotationModel = std::make_shared< ephemerides::GcrsToItrsRotationModel >(
                earthOrientationCalculator, earthRotationTime.first, getSp3EarthRotationBaseFrame( targetFrame.first ) );
    }

    for( const auto& state : sourceStateHistory )
    {
        const double epochSinceJ2000 = state.first + referenceEpochOffset;
        const double earthRotationEpoch = epochSinceJ2000 + earthRotationTime.second;
        Eigen::Vector6d stateInItrs;
        if( sourceFrame.first == Sp3FrameFamily::terrestrial )
        {
            stateInItrs = sourceFrame.second == "ITRS"
                    ? state.second
                    : reference_frames::convertStateBetweenItrfFrames( state.second, epochSinceJ2000, sourceFrame.second, "ITRF2014" );
        }
        else
        {
            stateInItrs = ephemerides::transformStateToTargetFrame( state.second, earthRotationEpoch, sourceRotationModel );
        }

        if( targetFrame.first == Sp3FrameFamily::terrestrial )
        {
            transformedStateHistory[ state.first ] = targetFrame.second == "ITRS"
                    ? stateInItrs
                    : reference_frames::convertStateBetweenItrfFrames( stateInItrs, epochSinceJ2000, "ITRF2014", targetFrame.second );
        }
        else
        {
            transformedStateHistory[ state.first ] =
                    ephemerides::transformStateToInertialOrientation( stateInItrs, earthRotationEpoch, targetRotationModel );
        }
    }
    return transformedStateHistory;
#else
    throw std::invalid_argument( "Cannot transform SP3 states from frame '" + fileContents.frameName + "' to frame '" +
                                 targetFrameOrientation + "': this Tudat build has no SOFA/IERS Earth-orientation support." );
#endif
}

//! Create tabulated ephemeris settings from parsed SP3 contents.
std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::shared_ptr< input_output::Sp3FileContents >& fileContents,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin,
                                                           const std::string& frameOrientation )
{
    if( fileContents == nullptr )
    {
        throw std::invalid_argument( "Cannot create SP3 ephemeris settings from null file contents." );
    }

    const std::string resolvedFrameOrientation = frameOrientation.empty( ) ? fileContents->frameName : frameOrientation;
    if( resolvedFrameOrientation.empty( ) )
    {
        throw std::invalid_argument( "The SP3 file does not declare a reference frame; provide frameOrientation explicitly." );
    }
    return tabulatedEphemerisSettings(
            getTransformedSp3StateHistory( *fileContents, satelliteId, resolvedFrameOrientation ), frameOrigin, resolvedFrameOrientation );
}

//! Read an SP3 file and create tabulated ephemeris settings for one of its satellites.
std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::string& fileName,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin,
                                                           const std::string& frameOrientation,
                                                           const double referenceJulianDay )
{
    return sp3EphemerisSettings( input_output::readSp3File( fileName, referenceJulianDay ), satelliteId, frameOrigin, frameOrientation );
}

}  // namespace simulation_setup

}  // namespace tudat
