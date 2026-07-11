#include "tudat/simulation/estimation_setup/createObservationCollection.h"
#include "tudat/simulation/environment_setup/createCameras.h"

namespace tudat
{

namespace observation_models
{

observation_models::ObservableType getObservableTypeFromTrackingDataString( const std::string& observableTypeString )
{
    if( observableTypeString == "DsnNWayAveragedDoppler" )
    {
        return observation_models::dsn_n_way_averaged_doppler;
    }
    else if( observableTypeString == "DopplerMeasuredFrequency" )
    {
        return observation_models::doppler_measured_frequency;
    }

    try
    {
        return observation_models::getObservableType( observableTypeString );
    }
    catch( const std::exception& e )
    {
        throw std::runtime_error( "Error when creating ObservationCollection from TrackingData: observable type '" + observableTypeString +
                                  "' is not recognised. Underlying error: " + e.what( ) );
    }
}

observation_models::LinkEnds getLinkEndsFromTrackingData(
        const std::vector< std::pair< std::pair< std::string, std::string >, std::string > >& rawLinkEnds )
{
    observation_models::LinkEnds linkEnds;
    for( const auto& linkEnd : rawLinkEnds )
    {
        observation_models::LinkEndType type = observation_models::getLinkEndTypeFromString( linkEnd.second );
        observation_models::LinkEndId id = observation_models::LinkEndId( linkEnd.first );
        if( !linkEnds.emplace( type, id ).second )
        {
            throw std::runtime_error( "Duplicate link-end type '" + linkEnd.second + "' in tracking data." );
        }
    }
    return linkEnds;
}

void checkTrackingDataLinkEnds( const observation_models::ObservableType observableType,
                                const observation_models::LinkEnds& linkEnds,
                                const observation_models::LinkEndType referenceLinkEnd )
{
    // Check that link ends are provided
    if( linkEnds.empty( ) )
    {
        throw std::runtime_error(
                "Error when creating ObservationCollection from TrackingData object: no link ends provided for observable " +
                observation_models::getObservableName( observableType ) + "." );
    }

    // Check that the reference link end is defined
    if( linkEnds.count( referenceLinkEnd ) == 0 )
    {
        throw std::runtime_error( "Error when creating ObservationCollection from TrackingData object: the reference link end " +
                                  observation_models::getLinkEndTypeString( referenceLinkEnd ) + " is not defined." );
    }

    // Missing check for link ends size and content depending on observable type?
}

bool shouldSkipObservationCollectionAncillarySetting( const std::string& ancillarySetting )
{
    return ancillarySetting == "Doppler base frequency" || ancillarySetting == "note2" || ancillarySetting == "catalog";
}

void resetTabulatedEphemerisFromTrackingSupplementaryStateHistory( const std::map< double, Eigen::Vector6d >& stateHistory,
                                                                   const std::shared_ptr< ephemerides::Ephemeris > ephemeris,
                                                                   const std::string& bodyName )
{
    if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, double > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory, std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, double > >( ephemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, double > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory, std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, double > >( ephemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory, std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, Time > >( ephemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory, std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) );
    }
    else
    {
        throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                  ", tabulated ephemeris type was not recognized." );
    }
}

void resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
        const std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        const std::string& bodyName )
{
    if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, double > >( rotationalEphemeris ) != nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, double > >( rotationalEphemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, double > >( rotationalEphemeris ) !=
             nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, double > >( rotationalEphemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, Time > >( rotationalEphemeris ) != nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, Time > >( rotationalEphemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, Time > >( rotationalEphemeris ) != nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, Time > >( rotationalEphemeris ) );
    }
    else
    {
        throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                  ", tabulated rotational ephemeris type was not recognized." );
    }
}

std::map< double, Eigen::Vector6d > getTranslationalStateHistoryWithVelocity(
        const data::TranslationalStateSupplementaryData& translationalStateSupplementaryData )
{
    std::map< double, Eigen::Vector6d > stateHistory = translationalStateSupplementaryData.getStateHistory( );

    if( stateHistory.size( ) > 0 && !translationalStateSupplementaryData.isVelocityDefined( ) )
    {
        if( stateHistory.size( ) == 1 )
        {
            stateHistory.begin( )->second.segment( 3, 3 ).setZero( );
        }
        else
        {
            for( auto stateIterator = stateHistory.begin( ); stateIterator != stateHistory.end( ); ++stateIterator )
            {
                Eigen::Vector3d velocity = Eigen::Vector3d::Zero( );

                auto nextStateIterator = stateIterator;
                ++nextStateIterator;

                if( stateIterator == stateHistory.begin( ) )
                {
                    velocity = ( nextStateIterator->second.segment( 0, 3 ) - stateIterator->second.segment( 0, 3 ) ) /
                            ( nextStateIterator->first - stateIterator->first );
                }
                else if( nextStateIterator == stateHistory.end( ) )
                {
                    auto previousStateIterator = stateIterator;
                    --previousStateIterator;

                    velocity = ( stateIterator->second.segment( 0, 3 ) - previousStateIterator->second.segment( 0, 3 ) ) /
                            ( stateIterator->first - previousStateIterator->first );
                }
                else
                {
                    auto previousStateIterator = stateIterator;
                    --previousStateIterator;

                    const double previousTime = previousStateIterator->first;
                    const double currentTime = stateIterator->first;
                    const double nextTime = nextStateIterator->first;

                    const double previousCoefficient =
                            ( currentTime - nextTime ) / ( ( previousTime - currentTime ) * ( previousTime - nextTime ) );
                    const double currentCoefficient = ( 2.0 * currentTime - previousTime - nextTime ) /
                            ( ( currentTime - previousTime ) * ( currentTime - nextTime ) );
                    const double nextCoefficient =
                            ( currentTime - previousTime ) / ( ( nextTime - previousTime ) * ( nextTime - currentTime ) );

                    velocity = previousCoefficient * previousStateIterator->second.segment( 0, 3 ) +
                            currentCoefficient * stateIterator->second.segment( 0, 3 ) +
                            nextCoefficient * nextStateIterator->second.segment( 0, 3 );
                }

                stateIterator->second.segment( 3, 3 ) = velocity;
            }
        }
    }

    return stateHistory;
}

void setTranslationalStateSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< data::TranslationalStateSupplementaryData > >&
                translationalStateSupplementaryData )
{
    for( auto it = translationalStateSupplementaryData.begin( ); it != translationalStateSupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;
        if( referencePointName != "" )
        {
            throw std::runtime_error(
                    "Error, reference point ID for setting ephemeris from tracking supplementary data must be empty, but found " +
                    referencePointName );
        }
        std::map< double, Eigen::Vector6d > stateHistory;
        std::vector< std::pair< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > > > inconsistentDuplicateStateHistoryEntries;
        std::string frameOrigin;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( i == 0 )
            {
                frameOrigin = it->second.at( i ).getFrameOrigin( );
            }
            else if( it->second.at( i ).getFrameOrigin( ) != frameOrigin )
            {
                throw std::runtime_error(
                        "Error, inconsistent frame origins found when setting translational state from tracking "
                        "supplementary data for body " +
                        bodyName + ". Found " + it->second.at( i ).getFrameOrigin( ) + " but expected " + frameOrigin + "." );
            }

            std::map< double, Eigen::Vector6d > currentStateHistory = getTranslationalStateHistoryWithVelocity( it->second.at( i ) );
            for( auto stateIterator = currentStateHistory.begin( ); stateIterator != currentStateHistory.end( ); ++stateIterator )
            {
                if( stateHistory.count( stateIterator->first ) == 0 )
                {
                    stateHistory[ stateIterator->first ] = stateIterator->second;
                }
                else if( !stateHistory.at( stateIterator->first ).isApprox( stateIterator->second, 0.0 ) )
                {
                    inconsistentDuplicateStateHistoryEntries.push_back( std::make_pair(
                            stateIterator->first, std::make_pair( stateHistory.at( stateIterator->first ), stateIterator->second ) ) );
                }
            }
        }

        std::shared_ptr< ephemerides::Ephemeris > ephemeris = bodies.at( bodyName )->getEphemeris( );
        if( ephemeris != nullptr )
        {
            if( ephemeris->getReferenceFrameOrigin( ) != frameOrigin )
            {
                throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                          ", existing ephemeris frame origin is " + ephemeris->getReferenceFrameOrigin( ) +
                                          " but supplementary data frame origin is " + frameOrigin + "." );
            }

            if( !ephemerides::isTabulatedEphemeris( ephemeris ) )
            {
                throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                          ", existing ephemeris is not tabulated." );
            }

            resetTabulatedEphemerisFromTrackingSupplementaryStateHistory( stateHistory, ephemeris, bodyName );
        }
        else
        {
            std::map< Time, Eigen::Vector6d > timeStateHistory;
            utilities::castMatrixMap< double, double, Time, double, 6, 1 >( stateHistory, timeStateHistory );
            bodies.at( bodyName )
                    ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< double, Time > >(
                            interpolators::createOneDimensionalInterpolator( timeStateHistory, interpolators::linearInterpolation( ) ),
                            frameOrigin ) );
        }
    }
}

void setRotationalStateSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< data::RotationalStateSupplementaryData > >&
                rotationalStateSupplementaryData )
{
    for( auto it = rotationalStateSupplementaryData.begin( ); it != rotationalStateSupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;
        if( referencePointName != "" )
        {
            throw std::runtime_error(
                    "Error, reference point ID for setting rotational ephemeris from tracking supplementary data must be empty, but "
                    "found " +
                    referencePointName );
        }
        std::map< double, Eigen::Vector7d > rotationalStateHistory;
        std::vector< std::pair< double, std::pair< Eigen::Vector7d, Eigen::Vector7d > > >
                inconsistentDuplicateRotationalStateHistoryEntries;
        std::string baseFrameOrientation;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( i == 0 )
            {
                baseFrameOrientation = it->second.at( i ).getBaseFrameOrientation( );
            }
            else if( it->second.at( i ).getBaseFrameOrientation( ) != baseFrameOrientation )
            {
                throw std::runtime_error(
                        "Error, inconsistent base frame orientations found when setting rotational state from tracking "
                        "supplementary data for body " +
                        bodyName + ". Found " + it->second.at( i ).getBaseFrameOrientation( ) + " but expected " + baseFrameOrientation +
                        "." );
            }

            const std::map< double, Eigen::Vector7d >& currentRotationalStateHistory = it->second.at( i ).getRotationalStateHistory( );
            for( auto stateIterator = currentRotationalStateHistory.begin( ); stateIterator != currentRotationalStateHistory.end( );
                 ++stateIterator )
            {
                if( rotationalStateHistory.count( stateIterator->first ) == 0 )
                {
                    rotationalStateHistory[ stateIterator->first ] = stateIterator->second;
                }
                else if( !rotationalStateHistory.at( stateIterator->first ).isApprox( stateIterator->second, 0.0 ) )
                {
                    inconsistentDuplicateRotationalStateHistoryEntries.push_back(
                            std::make_pair( stateIterator->first,
                                            std::make_pair( rotationalStateHistory.at( stateIterator->first ), stateIterator->second ) ) );
                }
            }
        }

        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris = bodies.at( bodyName )->getRotationalEphemeris( );
        if( rotationalEphemeris != nullptr )
        {
            if( rotationalEphemeris->getBaseFrameOrientation( ) != baseFrameOrientation )
            {
                throw std::runtime_error( "Error when setting tracking supplementary rotational data in body " + bodyName +
                                          ", existing rotational ephemeris base frame orientation is " +
                                          rotationalEphemeris->getBaseFrameOrientation( ) +
                                          " but supplementary data base frame orientation is " + baseFrameOrientation + "." );
            }

            if( !ephemerides::isTabulatedRotationalEphemeris( rotationalEphemeris ) )
            {
                throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                          ", existing rotational ephemeris is not tabulated." );
            }

            resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory( rotationalStateHistory, rotationalEphemeris, bodyName );
        }
        else
        {
            std::map< Time, Eigen::Vector7d > timeRotationalStateHistory;
            utilities::castMatrixMap< double, double, Time, double, 7, 1 >( rotationalStateHistory, timeRotationalStateHistory );
            bodies.at( bodyName )
                    ->setRotationalEphemeris( std::make_shared< ephemerides::TabulatedRotationalEphemeris< double, Time > >(
                            interpolators::createOneDimensionalInterpolator( timeRotationalStateHistory,
                                                                             interpolators::linearInterpolation( ) ),
                            baseFrameOrientation,
                            bodyName + "_Fixed" ) );
        }
    }
}

void setFrequencySupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > >&
                frequencySupplementaryData )
{
    for( auto it = frequencySupplementaryData.begin( ); it != frequencySupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;

        if( it->second.empty( ) )
        {
            continue;
        }

        std::vector< data::RampedFrequencySupplementaryData::FrequencyRamp > frequencyRamps;
        std::vector< std::map< double, double > > piecewiseConstantFrequencyHistories;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( it->second.at( i ) == nullptr )
            {
                throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName + ", reference point " +
                                          referencePointName + ": frequency data entry is null." );
            }

            if( i > 0 &&
                it->second.at( i )->getFrequencySupplementaryDataType( ) != it->second.at( 0 )->getFrequencySupplementaryDataType( ) )
            {
                throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName + ", reference point " +
                                          referencePointName + ": all frequency supplementary data entries must have the same type." );
            }

            if( it->second.at( i )->getFrequencySupplementaryDataType( ) == data::FrequencySupplementaryDataType::ramped_frequency )
            {
                std::shared_ptr< data::RampedFrequencySupplementaryData > rampedFrequencySupplementaryData =
                        std::dynamic_pointer_cast< data::RampedFrequencySupplementaryData >( it->second.at( i ) );
                if( rampedFrequencySupplementaryData == nullptr )
                {
                    throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName + ", reference point " +
                                              referencePointName +
                                              ": frequency data type is ramped, but derived object type is inconsistent." );
                }
                const std::vector< data::RampedFrequencySupplementaryData::FrequencyRamp >& currentFrequencyRamps =
                        rampedFrequencySupplementaryData->getFrequencyRamps( );
                frequencyRamps.insert( frequencyRamps.end( ), currentFrequencyRamps.begin( ), currentFrequencyRamps.end( ) );
            }
            else if( it->second.at( i )->getFrequencySupplementaryDataType( ) ==
                     data::FrequencySupplementaryDataType::piecewise_constant_frequency )
            {
                std::shared_ptr< data::PiecewiseConstantFrequencySupplementaryData > piecewiseConstantFrequencySupplementaryData =
                        std::dynamic_pointer_cast< data::PiecewiseConstantFrequencySupplementaryData >( it->second.at( i ) );
                if( piecewiseConstantFrequencySupplementaryData == nullptr )
                {
                    throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName + ", reference point " +
                                              referencePointName +
                                              ": frequency data type is piecewise constant, but derived object type is inconsistent." );
                }
                piecewiseConstantFrequencyHistories.push_back( piecewiseConstantFrequencySupplementaryData->getFrequencyHistory( ) );
            }
        }

        if( it->second.at( 0 )->getFrequencySupplementaryDataType( ) == data::FrequencySupplementaryDataType::ramped_frequency )
        {
            if( frequencyRamps.empty( ) )
            {
                throw std::runtime_error( "Error when setting ramped frequency supplementary data in body " + bodyName +
                                          ", reference point " + referencePointName + ": no frequency ramps were found." );
            }

            std::vector< Time > startTimes;
            std::vector< Time > endTimes;
            std::vector< double > rampRates;
            std::vector< double > startFrequencies;

            for( unsigned int i = 0; i < frequencyRamps.size( ); ++i )
            {
                startTimes.push_back( Time( frequencyRamps.at( i ).startTime_ ) );
                endTimes.push_back( Time( frequencyRamps.at( i ).endTime_ ) );
                rampRates.push_back( frequencyRamps.at( i ).frequencyRate_ );
                startFrequencies.push_back( frequencyRamps.at( i ).startFrequency_ );
            }

            std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > frequencyInterpolator =
                    std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                            startTimes, endTimes, rampRates, startFrequencies );

            if( referencePointName != "" )
            {
                if( bodies.at( bodyName )->getGroundStationMap( ).count( referencePointName ) == 0 )
                {
                    throw std::runtime_error( "Error when setting ramped frequency supplementary data in body " + bodyName +
                                              ", ground station " + referencePointName + " was not found." );
                }

                std::shared_ptr< ground_stations::GroundStation > groundStation =
                        bodies.at( bodyName )->getGroundStation( referencePointName );
                if( !groundStation->hasFrequencyCalculator( ) )
                {
                    groundStation->setTransmittingFrequencyCalculator( frequencyInterpolator );
                }
                else
                {
                    std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > existingFrequencyInterpolator =
                            std::dynamic_pointer_cast< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                                    groundStation->getTransmittingFrequencyCalculator( ) );
                    if( existingFrequencyInterpolator == nullptr )
                    {
                        throw std::runtime_error( "Error when setting ramped frequency supplementary data in body " + bodyName +
                                                  ", ground station " + referencePointName +
                                                  " already has a non-ramped frequency calculator." );
                    }
                    existingFrequencyInterpolator->addFrequencyInterpolator( frequencyInterpolator );
                }
            }
            else
            {
                if( bodies.at( bodyName )->getVehicleSystems( ) == nullptr )
                {
                    bodies.at( bodyName )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
                }

                std::shared_ptr< ground_stations::StationFrequencyInterpolator > existingFrequencyCalculator =
                        bodies.at( bodyName )->getVehicleSystems( )->getTransmittedFrequencyCalculator( );
                if( existingFrequencyCalculator == nullptr )
                {
                    bodies.at( bodyName )->getVehicleSystems( )->setTransmittedFrequencyCalculator( frequencyInterpolator );
                }
                else
                {
                    std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > existingFrequencyInterpolator =
                            std::dynamic_pointer_cast< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                                    existingFrequencyCalculator );
                    if( existingFrequencyInterpolator == nullptr )
                    {
                        throw std::runtime_error( "Error when setting ramped frequency supplementary data in body " + bodyName +
                                                  ", vehicle systems already contain a non-ramped frequency calculator." );
                    }
                    existingFrequencyInterpolator->addFrequencyInterpolator( frequencyInterpolator );
                }
            }
        }
    }
}

void setInstrumentSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::InstrumentSupplementaryData > > >&
                instrumentSupplementaryData )
{
    for( auto it = instrumentSupplementaryData.begin( ); it != instrumentSupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( it->second.at( i ) == nullptr )
            {
                throw std::runtime_error( "Error when setting instrument supplementary data in body " + bodyName + ", reference point " +
                                          referencePointName + ": instrument data entry is null." );
            }

            if( it->second.at( i )->getInstrumentSupplementaryDataType( ) == data::InstrumentSupplementaryDataType::camera_settings )
            {
                std::shared_ptr< data::CameraInstrumentSupplementaryData > cameraSupplementaryData =
                        std::dynamic_pointer_cast< data::CameraInstrumentSupplementaryData >( it->second.at( i ) );
                if( cameraSupplementaryData == nullptr )
                {
                    throw std::runtime_error( "Error when setting camera supplementary data in body " + bodyName + ", reference point " +
                                              referencePointName +
                                              ": instrument data type is camera, but derived object type is inconsistent." );
                }

                if( referencePointName != "" && referencePointName != cameraSupplementaryData->getCameraId( ) )
                {
                    throw std::runtime_error( "Error when setting camera supplementary data in body " + bodyName + ", reference point " +
                                              referencePointName + " does not match camera id " + cameraSupplementaryData->getCameraId( ) +
                                              "." );
                }

                if( bodies.at( bodyName )->getVehicleSystems( ) == nullptr )
                {
                    bodies.at( bodyName )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
                }

                if( bodies.at( bodyName )->getVehicleSystems( )->getCameraMap( ).count( cameraSupplementaryData->getCameraId( ) ) != 0 )
                {
                    throw std::runtime_error( "Error when setting camera supplementary data in body " + bodyName + ", camera " +
                                              cameraSupplementaryData->getCameraId( ) +
                                              " already exists. Overriding camera settings is not allowed." );
                }

                std::shared_ptr< system_models::PsfCameraProjectionModel > projectionModel =
                        std::make_shared< system_models::PsfCameraProjectionModel >( cameraSupplementaryData->getFocalLength( ),
                                                                                     cameraSupplementaryData->getPrincipalPoint( ),
                                                                                     cameraSupplementaryData->getKMatrix( ),
                                                                                     cameraSupplementaryData->getDistortionCoefficients( ),
                                                                                     cameraSupplementaryData->getMountingOffsets( ),
                                                                                     cameraSupplementaryData->getFieldOfViewBounds( ) );

                std::shared_ptr< simulation_setup::CameraSettings > cameraSettings = std::make_shared< simulation_setup::CameraSettings >(
                        cameraSupplementaryData->getCameraId( ), Eigen::Vector3d::Zero( ), projectionModel );
                simulation_setup::createCamera( bodies.at( bodyName ), cameraSettings );
            }
            else
            {
                throw std::runtime_error( "Error when setting instrument supplementary data in body " + bodyName + ", reference point " +
                                          referencePointName + ": unsupported instrument data type." );
            }
        }
    }
}

void setTrackingSupplementaryDataInBodies( simulation_setup::SystemOfBodies& bodies,
                                           const std::vector< data::TrackingSupplementaryData >& supplementaryData )
{
    std::map< std::pair< std::string, std::string >, std::vector< data::TranslationalStateSupplementaryData > >
            translationalStateSupplementaryData;
    std::map< std::pair< std::string, std::string >, std::vector< data::RotationalStateSupplementaryData > >
            rotationalStateSupplementaryData;
    std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > >
            frequencySupplementaryData;
    std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::InstrumentSupplementaryData > > >
            instrumentSupplementaryData;

    for( const data::TrackingSupplementaryData& currentSupplementaryData : supplementaryData )
    {
        const std::pair< std::string, std::string > bodyReferencePoint =
                std::make_pair( currentSupplementaryData.getBodyName( ), currentSupplementaryData.getReferencePointName( ) );

        if( !currentSupplementaryData.getTranslationalStateSupplementaryData( ).getStateHistory( ).empty( ) )
        {
            translationalStateSupplementaryData[ bodyReferencePoint ].push_back(
                    currentSupplementaryData.getTranslationalStateSupplementaryData( ) );
        }
        if( !currentSupplementaryData.getRotationalStateSupplementaryData( ).getRotationalStateHistory( ).empty( ) )
        {
            rotationalStateSupplementaryData[ bodyReferencePoint ].push_back(
                    currentSupplementaryData.getRotationalStateSupplementaryData( ) );
        }

        const std::vector< std::shared_ptr< data::FrequencySupplementaryData > >& currentFrequencySupplementaryData =
                currentSupplementaryData.getFrequencySupplementaryData( );
        if( !currentFrequencySupplementaryData.empty( ) )
        {
            frequencySupplementaryData[ bodyReferencePoint ].insert( frequencySupplementaryData[ bodyReferencePoint ].end( ),
                                                                     currentFrequencySupplementaryData.begin( ),
                                                                     currentFrequencySupplementaryData.end( ) );
        }

        const std::vector< std::shared_ptr< data::InstrumentSupplementaryData > >& currentInstrumentSupplementaryData =
                currentSupplementaryData.getInstrumentSupplementaryData( );
        if( !currentInstrumentSupplementaryData.empty( ) )
        {
            instrumentSupplementaryData[ bodyReferencePoint ].insert( instrumentSupplementaryData[ bodyReferencePoint ].end( ),
                                                                      currentInstrumentSupplementaryData.begin( ),
                                                                      currentInstrumentSupplementaryData.end( ) );
        }
    }

    setTranslationalStateSupplementaryDataInBodies( bodies, translationalStateSupplementaryData );
    setRotationalStateSupplementaryDataInBodies( bodies, rotationalStateSupplementaryData );
    setFrequencySupplementaryDataInBodies( bodies, frequencySupplementaryData );
    setInstrumentSupplementaryDataInBodies( bodies, instrumentSupplementaryData );
}

void setTrackingSupplementaryDataInBodies( simulation_setup::SystemOfBodies& bodies,
                                           const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData )
{
    std::vector< data::TrackingSupplementaryData > supplementaryDataValues;
    supplementaryDataValues.reserve( supplementaryData.size( ) );
    for( const std::shared_ptr< data::TrackingSupplementaryData >& currentSupplementaryData : supplementaryData )
    {
        if( currentSupplementaryData == nullptr )
        {
            throw std::runtime_error( "Error when setting tracking supplementary data in bodies, supplementary data entry is null." );
        }
        supplementaryDataValues.push_back( *currentSupplementaryData );
    }

    setTrackingSupplementaryDataInBodies( bodies, supplementaryDataValues );
}

}  // namespace observation_models

}  // namespace tudat
