/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROCESS_PSF_FILE_H
#define TUDAT_PROCESS_PSF_FILE_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/io/readPsfFile.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createCameras.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"

namespace tudat
{

namespace observation_models
{

//! Settings for converting raw PSF contents to Tudat camera models and pixel-coordinate observation sets.
class PsfFileObservationConversionSettings
{
public:
    PsfFileObservationConversionSettings(
            const std::string& receiverBodyName,
            const std::map< std::string, std::string >& imageNameToBodyName = std::map< std::string, std::string >( ) ):
        receiverBodyName_( receiverBodyName ), imageNameToBodyName_( imageNameToBodyName )
    {}

    //! Body that receives the optical observations, typically the spacecraft named by PSF SCID.
    std::string receiverBodyName_;

    //! Optional conversion from PSF IMG names to Tudat body names.
    std::map< std::string, std::string > imageNameToBodyName_;

    //! If true, unmapped PSF IMG names are used directly as Tudat body names.
    bool useRawImageNameAsBodyNameIfUnmapped_ = true;

    //! Use Z - ZC as the observed pixel/line value.
    bool useCorrectedPixelLine_ = true;

    //! Interpret TOB as end-of-exposure and store the mid-exposure epoch by subtracting EXPTIM/2.
    bool useMidExposureTime_ = true;

    //! Include pictures with nonzero PICDEL.
    bool includeDeletedPictures_ = false;

    //! Include $IM end-marker records as observations.
    bool includeEndMarkerRecords_ = false;

    //! If true, only measurements whose USE flag equals requiredUseFlag_ are converted.
    bool filterByUseFlag_ = false;

    //! USE flag retained when filterByUseFlag_ is true.
    int requiredUseFlag_ = 0;

    //! Body-fixed camera reference point position used when adding PSF cameras to the receiver body.
    Eigen::Vector3d bodyFixedCameraPosition_ = Eigen::Vector3d::Zero( );

    //! Use picture-specific RA/DEC/TWIST as a direct inertial-to-camera pointing source.
    bool usePicturePointing_ = true;
};

inline std::string convertStringToUpperCase( const std::string& input )
{
    std::string output = input;
    std::transform( output.begin( ), output.end( ), output.begin( ), []( unsigned char character ) {
        return static_cast< char >( std::toupper( character ) );
    } );
    return output;
}

inline int getPsfMonthNumber( const std::string& monthString )
{
    static const std::map< std::string, int > monthNumbers = { { "JAN", 1 }, { "FEB", 2 },  { "MAR", 3 },  { "APR", 4 },
                                                               { "MAY", 5 }, { "JUN", 6 },  { "JUL", 7 },  { "AUG", 8 },
                                                               { "SEP", 9 }, { "OCT", 10 }, { "NOV", 11 }, { "DEC", 12 } };

    const std::string upperMonthString = convertStringToUpperCase( monthString );
    if( monthNumbers.count( upperMonthString ) == 0 )
    {
        throw std::runtime_error( "Error when converting PSF epoch: unsupported month string '" + monthString + "'." );
    }
    return monthNumbers.at( upperMonthString );
}

namespace detail
{

template< typename TimeType = double >
TimeType convertPsfUtcStringToSecondsSinceJ2000( const std::string& utcString, const bool convertUtcToTdb = true )
{
    std::stringstream stream( utcString );
    int year = 0;
    std::string monthString;
    int day = 0;
    std::string timeString;
    stream >> year >> monthString >> day >> timeString;
    if( stream.fail( ) || timeString.empty( ) )
    {
        throw std::runtime_error( "Error when converting PSF epoch: malformed UTC string '" + utcString + "'." );
    }

    std::replace( timeString.begin( ), timeString.end( ), ':', ' ' );
    std::stringstream timeStream( timeString );
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;
    timeStream >> hour >> minute >> seconds;
    if( timeStream.fail( ) )
    {
        throw std::runtime_error( "Error when converting PSF epoch: malformed time field in '" + utcString + "'." );
    }

    const TimeType utcSecondsSinceJ2000 =
            basic_astrodynamics::DateTime( year, getPsfMonthNumber( monthString ), day, hour, minute, seconds ).epoch< TimeType >( );
    if( convertUtcToTdb )
    {
        return earth_orientation::TerrestrialTimeScaleConverter( ).getCurrentTime< TimeType >(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcSecondsSinceJ2000, Eigen::Vector3d::Zero( ) );
    }
    else
    {
        return utcSecondsSinceJ2000;
    }
}

template< typename TimeType = double >
TimeType getPsfPictureObservationTime( const input_output::psf::RawPsfFileImageContents& imageContents,
                                       const bool useMidExposureTime = true,
                                       const bool convertUtcToTdb = true )
{
    TimeType observationTime =
            convertPsfUtcStringToSecondsSinceJ2000< TimeType >( imageContents.endOfExposureTimeUtcString_, convertUtcToTdb );
    if( useMidExposureTime )
    {
        observationTime -= static_cast< TimeType >( 0.5 * imageContents.exposureTimeSeconds_ );
    }
    return observationTime;
}

//! Return the direct inertial-to-camera attitude encoded by PSF picture RA, DEC and TWIST.
inline Eigen::Quaterniond getPsfPictureRotationFromInertialToCameraFrame( const input_output::psf::RawPsfFileImageContents& imageContents )
{
    const double rightAscension = imageContents.rightAscensionDegrees_ * mathematical_constants::PI / 180.0;
    const double declination = imageContents.declinationDegrees_ * mathematical_constants::PI / 180.0;
    const double plateAngle = ( imageContents.twistDegrees_ - 90.0 ) * mathematical_constants::PI / 180.0;

    const Eigen::Vector3d boresight = ( Eigen::Vector3d( ) << std::cos( declination ) * std::cos( rightAscension ),
                                        std::cos( declination ) * std::sin( rightAscension ),
                                        std::sin( declination ) )
                                              .finished( );
    const Eigen::Vector3d east = ( Eigen::Vector3d( ) << -std::sin( rightAscension ), std::cos( rightAscension ), 0.0 ).finished( );
    const Eigen::Vector3d north = ( Eigen::Vector3d( ) << -std::sin( declination ) * std::cos( rightAscension ),
                                    -std::sin( declination ) * std::sin( rightAscension ),
                                    std::cos( declination ) )
                                          .finished( );

    Eigen::Matrix3d rotationFromCameraToInertial;
    rotationFromCameraToInertial.col( 0 ) = std::cos( plateAngle ) * east + std::sin( plateAngle ) * north;
    rotationFromCameraToInertial.col( 1 ) = -std::sin( plateAngle ) * east + std::cos( plateAngle ) * north;
    rotationFromCameraToInertial.col( 2 ) = boresight;

    return Eigen::Quaterniond( rotationFromCameraToInertial.transpose( ) );
}

template< typename TimeType = double >
std::function< Eigen::Quaterniond( const double ) > createNearestPsfPicturePointingFunction(
        const input_output::psf::RawPsfFileContents& psfFileContents,
        const std::string& cameraId,
        const bool useMidExposureTime = true,
        const bool convertUtcToTdb = true )
{
    std::map< double, Eigen::Quaterniond > picturePointing;
    for( const input_output::psf::RawPsfFileImageContents& imageContents : psfFileContents.images_ )
    {
        if( imageContents.cameraId_ == cameraId )
        {
            picturePointing[ static_cast< double >(
                    getPsfPictureObservationTime< TimeType >( imageContents, useMidExposureTime, convertUtcToTdb ) ) ] =
                    getPsfPictureRotationFromInertialToCameraFrame( imageContents );
        }
    }

    if( picturePointing.empty( ) )
    {
        throw std::runtime_error( "Error when creating PSF camera pointing function: no pictures found for camera '" + cameraId + "'." );
    }

    return [ = ]( const double time ) {
        auto upperIterator = picturePointing.lower_bound( time );
        if( upperIterator == picturePointing.begin( ) )
        {
            return upperIterator->second;
        }
        if( upperIterator == picturePointing.end( ) )
        {
            return std::prev( upperIterator )->second;
        }

        auto lowerIterator = std::prev( upperIterator );
        if( std::fabs( time - lowerIterator->first ) <= std::fabs( upperIterator->first - time ) )
        {
            return lowerIterator->second;
        }
        return upperIterator->second;
    };
}

}  // namespace detail

template< typename TimeType = double >
TimeType convertPsfUtcStringToSecondsSinceJ2000( const std::string& utcString )
{
    return detail::convertPsfUtcStringToSecondsSinceJ2000< TimeType >( utcString, true );
}

template< typename TimeType = double >
TimeType getPsfPictureObservationTime( const input_output::psf::RawPsfFileImageContents& imageContents,
                                       const bool useMidExposureTime = true )
{
    return detail::getPsfPictureObservationTime< TimeType >( imageContents, useMidExposureTime, true );
}

//! Return the direct inertial-to-camera attitude encoded by PSF picture RA, DEC and TWIST.
inline Eigen::Quaterniond getPsfPictureRotationFromInertialToCameraFrame( const input_output::psf::RawPsfFileImageContents& imageContents )
{
    return detail::getPsfPictureRotationFromInertialToCameraFrame( imageContents );
}

template< typename TimeType = double >
std::function< Eigen::Quaterniond( const double ) > createNearestPsfPicturePointingFunction(
        const input_output::psf::RawPsfFileContents& psfFileContents,
        const std::string& cameraId,
        const bool useMidExposureTime = true )
{
    return detail::createNearestPsfPicturePointingFunction< TimeType >( psfFileContents, cameraId, useMidExposureTime, true );
}

inline void addPsfCamerasToBody( const input_output::psf::RawPsfFileContents& psfFileContents,
                                 const std::shared_ptr< simulation_setup::Body >& receiverBody,
                                 const PsfFileObservationConversionSettings& conversionSettings )
{
    if( receiverBody == nullptr )
    {
        throw std::runtime_error( "Error when adding PSF cameras to body: receiver body is nullptr." );
    }

    for( const auto& cameraEntry : psfFileContents.cameraProperties_ )
    {
        const std::string& cameraId = cameraEntry.first;
        std::function< Eigen::Quaterniond( const double ) > pointingFunction = nullptr;
        if( conversionSettings.usePicturePointing_ )
        {
            pointingFunction =
                    createNearestPsfPicturePointingFunction< double >( psfFileContents, cameraId, conversionSettings.useMidExposureTime_ );
        }

        std::shared_ptr< simulation_setup::CameraSettings > cameraSettings = std::make_shared< simulation_setup::CameraSettings >(
                cameraId,
                Eigen::Vector3d::Zero( ),
                input_output::psf::createPsfCameraProjectionModel( cameraEntry.second ),
                conversionSettings.bodyFixedCameraPosition_,
                pointingFunction );
        simulation_setup::createCamera( receiverBody, cameraSettings );
    }
}

inline void addPsfCamerasToBodies( const input_output::psf::RawPsfFileContents& psfFileContents,
                                   const simulation_setup::SystemOfBodies& bodies,
                                   const PsfFileObservationConversionSettings& conversionSettings )
{
    addPsfCamerasToBody( psfFileContents, bodies.at( conversionSettings.receiverBodyName_ ), conversionSettings );
}

inline bool shouldConvertPsfMeasurement( const input_output::psf::RawPsfMeasurement& measurement,
                                         const PsfFileObservationConversionSettings& conversionSettings )
{
    if( !conversionSettings.includeEndMarkerRecords_ && measurement.opticalImageType_ == input_output::psf::OpticalImageType::end_marker )
    {
        return false;
    }

    if( conversionSettings.filterByUseFlag_ && measurement.useFlag_ != conversionSettings.requiredUseFlag_ )
    {
        return false;
    }

    return true;
}

inline std::string getTudatBodyNameForPsfMeasurement( const input_output::psf::RawPsfMeasurement& measurement,
                                                      const PsfFileObservationConversionSettings& conversionSettings )
{
    if( conversionSettings.imageNameToBodyName_.count( measurement.imageName_ ) != 0 )
    {
        return conversionSettings.imageNameToBodyName_.at( measurement.imageName_ );
    }

    const std::string upperImageName = convertStringToUpperCase( measurement.imageName_ );
    if( conversionSettings.imageNameToBodyName_.count( upperImageName ) != 0 )
    {
        return conversionSettings.imageNameToBodyName_.at( upperImageName );
    }

    if( conversionSettings.useRawImageNameAsBodyNameIfUnmapped_ )
    {
        return measurement.imageName_;
    }

    return "";
}

template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createPsfFileObservationDataset(
        const input_output::psf::RawPsfFileContents& psfFileContents,
        const PsfFileObservationConversionSettings& conversionSettings );

template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createPsfFileObservationDataset(
        const std::string& psfFile,
        const PsfFileObservationConversionSettings& conversionSettings );

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
createPsfFileObservationSets( const input_output::psf::RawPsfFileContents& psfFileContents,
                              const PsfFileObservationConversionSettings& conversionSettings )
{
    return createObservationCollection< ObservationScalarType, TimeType >(
                   createPsfFileObservationDataset< ObservationScalarType, TimeType >( psfFileContents, conversionSettings ) )
            ->getObservationsSets( );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createPsfFileObservationDataset(
        const input_output::psf::RawPsfFileContents& psfFileContents,
        const PsfFileObservationConversionSettings& conversionSettings )
{
    std::map< LinkEnds, std::vector< TimeType > > observationTimesMap;
    std::map< LinkEnds, std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observablesMap;
    std::map< LinkEnds, std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > weightsMap;

    for( const input_output::psf::RawPsfFileImageContents& imageContents : psfFileContents.images_ )
    {
        if( !conversionSettings.includeDeletedPictures_ && imageContents.pictureDeletionFlag_ != 0 )
        {
            continue;
        }

        const TimeType observationTime = getPsfPictureObservationTime< TimeType >( imageContents, conversionSettings.useMidExposureTime_ );

        for( const std::shared_ptr< input_output::psf::RawPsfMeasurement >& measurement : imageContents.measurements_ )
        {
            if( measurement == nullptr || !shouldConvertPsfMeasurement( *measurement, conversionSettings ) )
            {
                continue;
            }

            const std::string targetBodyName = getTudatBodyNameForPsfMeasurement( *measurement, conversionSettings );
            if( targetBodyName.empty( ) )
            {
                continue;
            }

            LinkEnds currentLinkEnds;
            currentLinkEnds[ transmitter ] = LinkEndId( targetBodyName );
            currentLinkEnds[ receiver ] = LinkEndId( conversionSettings.receiverBodyName_, imageContents.cameraId_ );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservable( 2 );
            const Eigen::Vector2d pixelLine =
                    conversionSettings.useCorrectedPixelLine_ ? measurement->getEffectivePixelLine( ) : measurement->observedPixelLine_;
            currentObservable( 0 ) = static_cast< ObservationScalarType >( pixelLine( 0 ) );
            currentObservable( 1 ) = static_cast< ObservationScalarType >( pixelLine( 1 ) );

            Eigen::Matrix< double, Eigen::Dynamic, 1 > currentWeights( 2 );
            for( int i = 0; i < 2; ++i )
            {
                const double sigma = measurement->sigmaPixelLine_( i );
                currentWeights( i ) = ( sigma > 0.0 ) ? 1.0 / ( sigma * sigma ) : 1.0;
            }

            observationTimesMap[ currentLinkEnds ].push_back( observationTime );
            observablesMap[ currentLinkEnds ].push_back( currentObservable );
            weightsMap[ currentLinkEnds ].push_back( currentWeights );
        }
    }

    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > observationDataset =
            std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );
    for( const auto& observationEntry : observablesMap )
    {
        const LinkEnds& currentLinkEnds = observationEntry.first;
        observationDataset->addObservationSet( pixel_coordinates,
                                               currentLinkEnds,
                                               observationEntry.second,
                                               observationTimesMap.at( currentLinkEnds ),
                                               receiver,
                                               std::vector< Eigen::VectorXd >( ),
                                               nullptr,
                                               nullptr,
                                               weightsMap.at( currentLinkEnds ) );
    }

    return observationDataset;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createPsfFileObservationDataset(
        const std::string& psfFile,
        const PsfFileObservationConversionSettings& conversionSettings )
{
    return createPsfFileObservationDataset< ObservationScalarType, TimeType >( input_output::psf::readPsfFile( psfFile ),
                                                                               conversionSettings );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createPsfFileObservationCollection(
        const input_output::psf::RawPsfFileContents& psfFileContents,
        const PsfFileObservationConversionSettings& conversionSettings )
{
    return createObservationCollection< ObservationScalarType, TimeType >(
            createPsfFileObservationDataset< ObservationScalarType, TimeType >( psfFileContents, conversionSettings ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createPsfFileObservationCollection(
        const std::string& psfFile,
        const PsfFileObservationConversionSettings& conversionSettings )
{
    return createObservationCollection< ObservationScalarType, TimeType >(
            createPsfFileObservationDataset< ObservationScalarType, TimeType >( psfFile, conversionSettings ) );
}

inline void addPsfCamerasToBodies( const std::string& psfFile,
                                   const simulation_setup::SystemOfBodies& bodies,
                                   const PsfFileObservationConversionSettings& conversionSettings )
{
    addPsfCamerasToBodies( input_output::psf::readPsfFile( psfFile ), bodies, conversionSettings );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_PROCESS_PSF_FILE_H
