/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READ_PSF_FILE_H
#define TUDAT_READ_PSF_FILE_H

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace tudat
{

namespace input_output
{

namespace psf
{

enum class OpticalImageType { star, planet, satellite, rock, end_marker, unknown };

class RawPsfMeasurement
{
public:
    virtual ~RawPsfMeasurement( ) {}

    OpticalImageType opticalImageType_ = OpticalImageType::unknown;

    std::string imageName_;     // raw IMG (string, without quotes)
    std::int64_t imageId_ = 0;  // raw IMGID
    int useFlag_ = 0;           // raw USE

    Eigen::Vector2d observedPixelLine_ = Eigen::Vector2d::Zero( );  // Z
    Eigen::Vector2d localCorrection_ = Eigen::Vector2d::Zero( );    // ZC
    Eigen::Vector2d sigmaPixelLine_ = Eigen::Vector2d::Zero( );     // SIG

    Eigen::Vector2d getEffectivePixelLine( ) const;
};

class RawPsfStarMeasurement : public RawPsfMeasurement
{
public:
    double rightAscensionDegrees_ = 0.0;  // STRA
    double declinationDegrees_ = 0.0;     // STDEC
};

class RawPsfFileImageContents
{
public:
    RawPsfFileImageContents( ) = default;

    std::string pictureName_;                 // PICNM
    int pictureNumber_ = 0;                   // PICNO
    std::string endOfExposureTimeUtcString_;  // TOB
    std::string cameraId_;                    // CAMERA

    double exposureTimeSeconds_ = 0.0;  // EXPTIM
    int pictureDeletionFlag_ = 0;       // PICDEL

    double rightAscensionDegrees_ = 0.0;  // RA
    double declinationDegrees_ = 0.0;     // DEC
    double twistDegrees_ = 0.0;           // TWIST

    std::vector< std::shared_ptr< RawPsfMeasurement > > measurements_;
};

class RawPsfFileContents
{
public:
    RawPsfFileContents( ) = default;

    // $ID
    std::string spacecraftId_;  // SCID
    int numberOfCameras_ = 0;   // NCAM
    int equinox_ = 0;           // EQUNOX

    std::string psfId_;                       // PSFID
    std::string psfProgram_;                  // PSFPRG
    std::string psfGenerationTimeUtcString_;  // PSFTIM
    std::vector< std::string > psfComments_;  // PSFCOM list

    // $CAM
    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > trackingSupplementaryData_;

    // $PIC/$IM
    std::vector< RawPsfFileImageContents > images_;
};

//! Retrieve the camera supplementary data read from a PSF $CAM block.
std::shared_ptr< const data::CameraInstrumentSupplementaryData > getPsfCameraInstrumentSupplementaryData(
        const RawPsfFileContents& psfFileContents,
        const std::string& cameraId );

double getPsfPictureObservationTime( const RawPsfFileImageContents& imageContents, const bool useMidExposureTime = true );

RawPsfFileContents readRawPsfFile( const std::string& psfFile );

std::string getTudatBodyNameForPsfMeasurement( const RawPsfMeasurement& measurement,
                                               const std::map< std::string, std::string >& imageNameToBodyName,
                                               const bool useRawImageNameAsBodyNameIfUnmapped );

bool shouldConvertPsfMeasurement( const RawPsfMeasurement& measurement,
                                  const bool includeEndMarkerRecords,
                                  const bool filterByUseFlag,
                                  const int requiredUseFlag );

std::vector< std::shared_ptr< data::TrackingSupplementaryData > > getPsfTrackingSupplementaryDataForReceiver(
        const RawPsfFileContents& psfFileContents,
        const std::string& receiverBodyName );

inline data::PlainLinkDefinition getPsfLinkEnds( const std::string& targetBodyName,
                                                 const std::string& receiverBodyName,
                                                 const std::string& cameraId )
{
    return data::PlainLinkDefinition( { std::make_pair( std::make_pair( targetBodyName, "" ), "transmitter" ),
                                        std::make_pair( std::make_pair( receiverBodyName, cameraId ), "receiver" ) } );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertRawPsfFiles( const std::vector< RawPsfFileContents >& rawPsfDataVector,
                    const std::string& receiverBodyName = "",
                    const std::map< std::string, std::string >& imageNameToBodyName = std::map< std::string, std::string >( ),
                    const bool useRawImageNameAsBodyNameIfUnmapped = true,
                    const bool useCorrectedPixelLine = true,
                    const bool useMidExposureTime = true,
                    const bool includeDeletedPictures = false,
                    const bool includeEndMarkerRecords = false,
                    const bool filterByUseFlag = false,
                    const int requiredUseFlag = 0 )
{
    std::map< data::PlainLinkDefinition, std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observationsMap;
    std::map< data::PlainLinkDefinition, std::vector< TimeType > > observationTimesMap;
    std::map< data::PlainLinkDefinition, std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > weightsMap;

    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > trackingSupplementaryDataSets;
    std::set< std::pair< std::string, std::string > > addedSupplementaryData;

    for( const RawPsfFileContents& psfFileContents : rawPsfDataVector )
    {
        const std::string currentReceiverBodyName = receiverBodyName.empty( ) ? psfFileContents.spacecraftId_ : receiverBodyName;

        for( const std::shared_ptr< data::TrackingSupplementaryData >& supplementaryData :
             getPsfTrackingSupplementaryDataForReceiver( psfFileContents, currentReceiverBodyName ) )
        {
            const std::pair< std::string, std::string > key =
                    std::make_pair( supplementaryData->getBodyName( ), supplementaryData->getReferencePointName( ) );
            if( addedSupplementaryData.insert( key ).second )
            {
                trackingSupplementaryDataSets.push_back( supplementaryData );
            }
        }

        for( const RawPsfFileImageContents& imageContents : psfFileContents.images_ )
        {
            if( !includeDeletedPictures && imageContents.pictureDeletionFlag_ != 0 )
            {
                continue;
            }

            const TimeType observationTime = static_cast< TimeType >( getPsfPictureObservationTime( imageContents, useMidExposureTime ) );

            for( const std::shared_ptr< RawPsfMeasurement >& measurement : imageContents.measurements_ )
            {
                if( measurement == nullptr ||
                    !shouldConvertPsfMeasurement( *measurement, includeEndMarkerRecords, filterByUseFlag, requiredUseFlag ) )
                {
                    continue;
                }

                const std::string targetBodyName =
                        getTudatBodyNameForPsfMeasurement( *measurement, imageNameToBodyName, useRawImageNameAsBodyNameIfUnmapped );
                if( targetBodyName.empty( ) )
                {
                    continue;
                }

                const data::PlainLinkDefinition linkEnds =
                        getPsfLinkEnds( targetBodyName, currentReceiverBodyName, imageContents.cameraId_ );
                const Eigen::Vector2d pixelLine =
                        useCorrectedPixelLine ? measurement->getEffectivePixelLine( ) : measurement->observedPixelLine_;

                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservable( 2 );
                currentObservable( 0 ) = static_cast< ObservationScalarType >( pixelLine( 0 ) );
                currentObservable( 1 ) = static_cast< ObservationScalarType >( pixelLine( 1 ) );

                Eigen::Matrix< double, Eigen::Dynamic, 1 > currentWeights( 2 );
                for( int i = 0; i < 2; ++i )
                {
                    const double sigma = measurement->sigmaPixelLine_( i );
                    currentWeights( i ) = ( sigma > 0.0 ) ? 1.0 / ( sigma * sigma ) : 1.0;
                }

                observationsMap[ linkEnds ].push_back( currentObservable );
                observationTimesMap[ linkEnds ].push_back( observationTime );
                weightsMap[ linkEnds ].push_back( currentWeights );
            }
        }
    }

    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;
    for( const auto& observationEntry : observationsMap )
    {
        auto trackingData =
                std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >( "PixelCoordinates",
                                                                                           observationEntry.first,
                                                                                           observationEntry.second,
                                                                                           observationTimesMap.at( observationEntry.first ),
                                                                                           "receiver",
                                                                                           "UTC" );
        trackingData->setObservationWeights( weightsMap.at( observationEntry.first ) );
        trackingDataSets.push_back( trackingData );
    }

    return std::make_pair( trackingDataSets, trackingSupplementaryDataSets );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readPsfFiles( const std::vector< std::string >& psfFileNames,
              const std::string& receiverBodyName = "",
              const std::map< std::string, std::string >& imageNameToBodyName = std::map< std::string, std::string >( ),
              const bool useRawImageNameAsBodyNameIfUnmapped = true,
              const bool useCorrectedPixelLine = true,
              const bool useMidExposureTime = true,
              const bool includeDeletedPictures = false,
              const bool includeEndMarkerRecords = false,
              const bool filterByUseFlag = false,
              const int requiredUseFlag = 0 )
{
    std::vector< RawPsfFileContents > rawPsfDataVector;
    rawPsfDataVector.reserve( psfFileNames.size( ) );
    for( const std::string& psfFileName : psfFileNames )
    {
        rawPsfDataVector.push_back( readRawPsfFile( psfFileName ) );
    }

    return convertRawPsfFiles< ObservationScalarType, TimeType >( rawPsfDataVector,
                                                                  receiverBodyName,
                                                                  imageNameToBodyName,
                                                                  useRawImageNameAsBodyNameIfUnmapped,
                                                                  useCorrectedPixelLine,
                                                                  useMidExposureTime,
                                                                  includeDeletedPictures,
                                                                  includeEndMarkerRecords,
                                                                  filterByUseFlag,
                                                                  requiredUseFlag );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readPsfFile( const std::string& psfFileName,
             const std::string& receiverBodyName = "",
             const std::map< std::string, std::string >& imageNameToBodyName = std::map< std::string, std::string >( ),
             const bool useRawImageNameAsBodyNameIfUnmapped = true,
             const bool useCorrectedPixelLine = true,
             const bool useMidExposureTime = true,
             const bool includeDeletedPictures = false,
             const bool includeEndMarkerRecords = false,
             const bool filterByUseFlag = false,
             const int requiredUseFlag = 0 )
{
    return readPsfFiles< ObservationScalarType, TimeType >( std::vector< std::string >( { psfFileName } ),
                                                            receiverBodyName,
                                                            imageNameToBodyName,
                                                            useRawImageNameAsBodyNameIfUnmapped,
                                                            useCorrectedPixelLine,
                                                            useMidExposureTime,
                                                            includeDeletedPictures,
                                                            includeEndMarkerRecords,
                                                            filterByUseFlag,
                                                            requiredUseFlag );
}

}  // namespace psf

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_PSF_FILE_H
