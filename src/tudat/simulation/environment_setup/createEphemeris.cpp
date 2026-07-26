/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <cctype>
#include <mutex>
#include <set>

#include "tudat/astro/basic_astro/physicalConstants.h"
#if TUDAT_BUILD_WITH_SOFA_INTERFACE
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#endif
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/io/readSp3File.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

enum class Sp3FrameFamily { terrestrial, gcrs, j2000, eclipj2000, unknown };

struct Sp3FrameDefinition {
    Sp3FrameFamily family;
    std::string canonicalName;
};

std::string normalizeSp3FrameName( std::string frameName )
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

bool hasTwoDigitYearSuffix( const std::string& frameName )
{
    return frameName.size( ) == 5 && std::isdigit( static_cast< unsigned char >( frameName.at( 3 ) ) ) &&
            std::isdigit( static_cast< unsigned char >( frameName.at( 4 ) ) );
}

std::string canonicalItrfNameFromTwoDigitYear( const std::string& year )
{
    const int numericalYear = std::stoi( year );
    return "ITRF" + std::to_string( numericalYear >= 88 ? 1900 + numericalYear : 2000 + numericalYear );
}

Sp3FrameDefinition identifySp3Frame( const std::string& rawFrameName )
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
        if( ( year.size( ) == 2 || year.size( ) == 4 ) &&
            std::all_of( year.begin( ), year.end( ), []( const unsigned char character ) { return std::isdigit( character ); } ) )
        {
            return { Sp3FrameFamily::terrestrial, year.size( ) == 2 ? canonicalItrfNameFromTwoDigitYear( year ) : "ITRF" + year };
        }
    }
    if( hasTwoDigitYearSuffix( frameName ) )
    {
        const std::string prefix = frameName.substr( 0, 3 );
        if( prefix == "IER" || prefix == "ITR" || prefix == "IGS" || prefix == "IGB" || prefix == "IGC" || prefix == "SLR" )
        {
            return { Sp3FrameFamily::terrestrial, canonicalItrfNameFromTwoDigitYear( frameName.substr( 3, 2 ) ) };
        }
    }
    return { Sp3FrameFamily::unknown, frameName };
}

bool isSupportedItrfHelmertFrame( const std::string& frameName )
{
    static const std::set< std::string > supportedFrames = { "ITRF1988", "ITRF1989", "ITRF1990", "ITRF1991", "ITRF1992",
                                                             "ITRF1993", "ITRF1994", "ITRF1996", "ITRF1997", "ITRF2000",
                                                             "ITRF2005", "ITRF2008", "ITRF2014" };
    return supportedFrames.count( frameName ) > 0;
}

std::string getTudatItrfFrameName( const std::string& canonicalFrameName )
{
    if( canonicalFrameName.size( ) == 8 && canonicalFrameName.rfind( "ITRF19", 0 ) == 0 )
    {
        return "ITRF" + canonicalFrameName.substr( 6, 2 );
    }
    return canonicalFrameName;
}

Eigen::Vector6d transformStateFromItrf2014( const Eigen::Vector6d& state, const double epochSinceJ2000, const std::string& targetFrame )
{
    const std::string tudatTargetFrame = getTudatItrfFrameName( targetFrame );
    const Eigen::Matrix3d matrixAtReferenceEpoch = reference_frames::getItrf2014ToArbitraryItrfRotationMatrix( tudatTargetFrame );
    const Eigen::Matrix3d matrixDerivative = reference_frames::getItrf2014ToArbitraryItrfRotationMatrixDerivative( tudatTargetFrame );
    const Eigen::Vector6d translation = reference_frames::getItrf2014ToArbitraryItrfTranslation( tudatTargetFrame );
    const double timeFromReferenceEpoch = epochSinceJ2000 - 10.0 * physical_constants::JULIAN_YEAR;
    const Eigen::Matrix3d matrixAtEpoch = matrixAtReferenceEpoch + matrixDerivative * timeFromReferenceEpoch;

    Eigen::Vector6d transformedState;
    transformedState.segment< 3 >( 0 ) = translation.segment< 3 >( 0 ) + translation.segment< 3 >( 3 ) * timeFromReferenceEpoch +
            matrixAtEpoch * state.segment< 3 >( 0 );
    transformedState.segment< 3 >( 3 ) =
            translation.segment< 3 >( 3 ) + matrixDerivative * state.segment< 3 >( 0 ) + matrixAtEpoch * state.segment< 3 >( 3 );
    return transformedState;
}

Eigen::Vector6d transformStateToItrf2014( const Eigen::Vector6d& state, const double epochSinceJ2000, const std::string& sourceFrame )
{
    const std::string tudatSourceFrame = getTudatItrfFrameName( sourceFrame );
    const Eigen::Matrix3d matrixAtReferenceEpoch = reference_frames::getItrf2014ToArbitraryItrfRotationMatrix( tudatSourceFrame );
    const Eigen::Matrix3d matrixDerivative = reference_frames::getItrf2014ToArbitraryItrfRotationMatrixDerivative( tudatSourceFrame );
    const Eigen::Vector6d translation = reference_frames::getItrf2014ToArbitraryItrfTranslation( tudatSourceFrame );
    const double timeFromReferenceEpoch = epochSinceJ2000 - 10.0 * physical_constants::JULIAN_YEAR;
    const Eigen::Matrix3d matrixAtEpoch = matrixAtReferenceEpoch + matrixDerivative * timeFromReferenceEpoch;
    const Eigen::Matrix3d inverseMatrixAtEpoch = matrixAtEpoch.inverse( );

    Eigen::Vector6d transformedState;
    transformedState.segment< 3 >( 0 ) = inverseMatrixAtEpoch *
            ( state.segment< 3 >( 0 ) - translation.segment< 3 >( 0 ) - translation.segment< 3 >( 3 ) * timeFromReferenceEpoch );
    transformedState.segment< 3 >( 3 ) = inverseMatrixAtEpoch *
            ( state.segment< 3 >( 3 ) - translation.segment< 3 >( 3 ) - matrixDerivative * transformedState.segment< 3 >( 0 ) );
    return transformedState;
}

Eigen::Vector6d transformStateBetweenItrfFrames( const Eigen::Vector6d& state,
                                                 const double epochSinceJ2000,
                                                 const std::string& sourceFrame,
                                                 const std::string& targetFrame )
{
    if( sourceFrame == targetFrame )
    {
        return state;
    }
    const Eigen::Vector6d stateInItrf2014 =
            sourceFrame == "ITRF2014" ? state : transformStateToItrf2014( state, epochSinceJ2000, sourceFrame );
    return targetFrame == "ITRF2014" ? stateInItrf2014 : transformStateFromItrf2014( stateInItrf2014, epochSinceJ2000, targetFrame );
}

#if TUDAT_BUILD_WITH_SOFA_INTERFACE
std::pair< basic_astrodynamics::TimeScales, double > getEarthRotationTimeScaleAndOffset( const std::string& timeScale )
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

std::string getEarthRotationBaseFrame( const Sp3FrameFamily family )
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

std::shared_ptr< ephemerides::RotationalEphemeris > createSp3EarthRotationModel( const Sp3FrameFamily inertialFrame,
                                                                                 const basic_astrodynamics::TimeScales inputTimeScale )
{
    static const std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            earth_orientation::createStandardEarthOrientationCalculator( );
    return std::make_shared< ephemerides::GcrsToItrsRotationModel >(
            earthOrientationCalculator, inputTimeScale, getEarthRotationBaseFrame( inertialFrame ) );
}

std::mutex& getSp3EarthRotationMutex( )
{
    static std::mutex earthRotationMutex;
    return earthRotationMutex;
}

Eigen::Vector6d rotateStateFromItrs( const Eigen::Vector6d& state,
                                     const double epoch,
                                     const std::shared_ptr< ephemerides::RotationalEphemeris >& rotationModel )
{
    return ephemerides::transformStateToFrameFromRotations(
            state, rotationModel->getRotationToBaseFrame( epoch ), rotationModel->getDerivativeOfRotationToBaseFrame( epoch ) );
}

Eigen::Vector6d rotateStateToItrs( const Eigen::Vector6d& state,
                                   const double epoch,
                                   const std::shared_ptr< ephemerides::RotationalEphemeris >& rotationModel )
{
    return ephemerides::transformStateToFrameFromRotations(
            state, rotationModel->getRotationToTargetFrame( epoch ), rotationModel->getDerivativeOfRotationToTargetFrame( epoch ) );
}
#endif

}  // namespace

std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::shared_ptr< input_output::Sp3FileContents >& fileContents,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin,
                                                           const std::string& frameOrientation )
{
    if( fileContents == nullptr )
    {
        throw std::invalid_argument( "Cannot create SP3 ephemeris settings from null file contents." );
    }

    const auto satelliteIterator = fileContents->satelliteStates.find( satelliteId );
    if( satelliteIterator == fileContents->satelliteStates.end( ) )
    {
        throw std::invalid_argument( "Satellite '" + satelliteId + "' is not present in the SP3 file." );
    }
    if( satelliteIterator->second.empty( ) )
    {
        throw std::invalid_argument( "Satellite '" + satelliteId + "' has no states in the SP3 file." );
    }

    std::map< double, Eigen::Vector6d > sourceStateHistory;
    for( const auto& state : satelliteIterator->second )
    {
        if( state.second.size( ) != 6 )
        {
            throw std::runtime_error( "SP3 state for satellite '" + satelliteId + "' does not have six components." );
        }

        const Eigen::Vector6d fixedSizeState = state.second;
        if( !fixedSizeState.allFinite( ) )
        {
            throw std::runtime_error( "SP3 ephemeris settings require finite position and velocity records for satellite '" + satelliteId +
                                      "'. The selected SP3 file contains a missing value or does not include velocity records." );
        }
        sourceStateHistory[ state.first ] = fixedSizeState;
    }

    const std::string resolvedFrameOrientation = frameOrientation.empty( ) ? fileContents->frameName : frameOrientation;
    if( resolvedFrameOrientation.empty( ) )
    {
        throw std::invalid_argument( "The SP3 file does not declare a reference frame; provide frameOrientation explicitly." );
    }

    if( frameOrientation.empty( ) || normalizeSp3FrameName( resolvedFrameOrientation ) == normalizeSp3FrameName( fileContents->frameName ) )
    {
        return tabulatedEphemerisSettings( sourceStateHistory, frameOrigin, resolvedFrameOrientation );
    }
    if( fileContents->frameName.empty( ) )
    {
        throw std::invalid_argument( "Cannot transform SP3 states because the file does not declare their source frame." );
    }

    const Sp3FrameDefinition sourceFrame = identifySp3Frame( fileContents->frameName );
    const Sp3FrameDefinition targetFrame = identifySp3Frame( resolvedFrameOrientation );
    if( sourceFrame.family == Sp3FrameFamily::unknown || targetFrame.family == Sp3FrameFamily::unknown )
    {
        throw std::invalid_argument( "Cannot transform SP3 states from frame '" + fileContents->frameName + "' to frame '" +
                                     resolvedFrameOrientation + "': one or both frame tags are not recognized." );
    }

    std::map< double, Eigen::Vector6d > transformedStateHistory;
    const double referenceEpochOffset =
            ( fileContents->referenceJulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;

    if( sourceFrame.family == Sp3FrameFamily::terrestrial && targetFrame.family == Sp3FrameFamily::terrestrial )
    {
        if( sourceFrame.canonicalName == targetFrame.canonicalName || sourceFrame.canonicalName == "ITRS" ||
            targetFrame.canonicalName == "ITRS" )
        {
            return tabulatedEphemerisSettings( sourceStateHistory, frameOrigin, resolvedFrameOrientation );
        }
        if( !isSupportedItrfHelmertFrame( sourceFrame.canonicalName ) || !isSupportedItrfHelmertFrame( targetFrame.canonicalName ) )
        {
            throw std::invalid_argument( "Cannot transform SP3 states from terrestrial frame '" + fileContents->frameName + "' to '" +
                                         resolvedFrameOrientation +
                                         "': no unambiguous Tudat frame-realization transformation is available." );
        }
        for( const auto& state : sourceStateHistory )
        {
            transformedStateHistory[ state.first ] = transformStateBetweenItrfFrames(
                    state.second, state.first + referenceEpochOffset, sourceFrame.canonicalName, targetFrame.canonicalName );
        }
        return tabulatedEphemerisSettings( transformedStateHistory, frameOrigin, resolvedFrameOrientation );
    }

#if TUDAT_BUILD_WITH_SOFA_INTERFACE
    const std::lock_guard< std::mutex > earthRotationLock( getSp3EarthRotationMutex( ) );
    const auto earthRotationTime = getEarthRotationTimeScaleAndOffset( fileContents->timeScale );
    std::shared_ptr< ephemerides::RotationalEphemeris > sourceRotationModel;
    std::shared_ptr< ephemerides::RotationalEphemeris > targetRotationModel;
    if( sourceFrame.family != Sp3FrameFamily::terrestrial )
    {
        sourceRotationModel = createSp3EarthRotationModel( sourceFrame.family, earthRotationTime.first );
    }
    if( targetFrame.family != Sp3FrameFamily::terrestrial )
    {
        targetRotationModel = createSp3EarthRotationModel( targetFrame.family, earthRotationTime.first );
    }

    for( const auto& state : sourceStateHistory )
    {
        const double earthRotationEpoch = state.first + referenceEpochOffset + earthRotationTime.second;
        Eigen::Vector6d stateInItrs = sourceFrame.family == Sp3FrameFamily::terrestrial
                ? state.second
                : rotateStateToItrs( state.second, earthRotationEpoch, sourceRotationModel );
        transformedStateHistory[ state.first ] = targetFrame.family == Sp3FrameFamily::terrestrial
                ? stateInItrs
                : rotateStateFromItrs( stateInItrs, earthRotationEpoch, targetRotationModel );
    }

    return tabulatedEphemerisSettings( transformedStateHistory, frameOrigin, resolvedFrameOrientation );
#else
    throw std::invalid_argument( "Cannot transform SP3 states from frame '" + fileContents->frameName + "' to frame '" +
                                 resolvedFrameOrientation + "': this Tudat build has no SOFA/IERS Earth-orientation support." );
#endif
}

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
