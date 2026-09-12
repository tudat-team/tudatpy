/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ANCILLARYSETTINGS_H
#define TUDAT_ANCILLARYSETTINGS_H

#include <Eigen/Core>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/serialization/core.h"
#include "tudat/io/serialization/file_io_declarations.h"

namespace tudat
{

namespace observation_models
{

enum ObservationAncillarySimulationVariable {
    link_ends_delays,
    frequency_bands,
    reception_reference_frequency_band,
    doppler_integration_time,
    doppler_reference_frequency,
    sequential_range_lowest_ranging_component,
    range_conversion_factor,
};

enum ObservationIntermediateSimulationVariable { transmitter_frequency_intermediate, received_frequency_intermediate };

struct ObservationAncillarySimulationSettings {
public:
    ObservationAncillarySimulationSettings( ) {}

    virtual ~ObservationAncillarySimulationSettings( ) {}

    //! Save ancillary settings to a JSON file
    TUDAT_DECLARE_FILE_IO( ObservationAncillarySimulationSettings )

    void setAncillaryDoubleData( const ObservationAncillarySimulationVariable& variableType, const double variable )
    {
        switch( variableType )
        {
            case doppler_integration_time:
            case doppler_reference_frequency:
            case reception_reference_frequency_band:
            case sequential_range_lowest_ranging_component:
            case range_conversion_factor:
                doubleData_[ variableType ] = variable;
                break;
            default:
                throw std::runtime_error(
                        "Error when setting double ancillary observation "
                        "data; could not set type " +
                        getAncillaryDataName( variableType ) );
        }
    }

    void setAncillaryDoubleVectorData( const ObservationAncillarySimulationVariable& variableType, const std::vector< double >& variable )
    {
        switch( variableType )
        {
            case link_ends_delays:
            case frequency_bands:
                doubleVectorData_[ variableType ] = variable;
                break;
            default:
                throw std::runtime_error(
                        "Error when setting double vector ancillary "
                        "observation data; could not set type " +
                        getAncillaryDataName( variableType ) );
        }
    }

    double getAncillaryDoubleData( const ObservationAncillarySimulationVariable& variableType, const bool throwException = true )
    {
        double returnVariable = TUDAT_NAN;
        try
        {
            switch( variableType )
            {
                case doppler_integration_time:
                case doppler_reference_frequency:
                case reception_reference_frequency_band:
                case sequential_range_lowest_ranging_component:
                case range_conversion_factor:
                    returnVariable = doubleData_.at( variableType );
                    break;
                default:
                    if( throwException )
                    {
                        throw std::runtime_error(
                                "Error when getting double ancillary observation "
                                "data; could not retrieve type " +
                                getAncillaryDataName( variableType ) );
                    }
                    break;
            }
        }
        catch( ... )
        {
            if( throwException )
            {
                throw std::runtime_error(
                        "Error when getting double ancillary observation "
                        "data; could not retrieve type " +
                        getAncillaryDataName( variableType ) );
            }
        }
        return returnVariable;
    }

    std::vector< double > getAncillaryDoubleVectorData( const ObservationAncillarySimulationVariable& variableType,
                                                        const bool throwException = true )
    {
        std::vector< double > returnVariable;
        try
        {
            switch( variableType )
            {
                case link_ends_delays:
                case frequency_bands:
                    returnVariable = doubleVectorData_.at( variableType );
                    break;
                default:
                    if( throwException )
                    {
                        throw std::runtime_error(
                                "Error when getting double vector ancillary "
                                "observation data; could not retrieve type " +
                                getAncillaryDataName( variableType ) );
                    }
                    break;
            }
        }
        catch( ... )
        {
            if( throwException )
            {
                throw std::runtime_error(
                        "Error when getting double vector ancillary observation "
                        "data; could not retrieve type " +
                        getAncillaryDataName( variableType ) );
            }
        }
        return returnVariable;
    }

    std::string getAncillaryDataName( const ObservationAncillarySimulationVariable& variableType )
    {
        std::string name;

        switch( variableType )
        {
            case link_ends_delays:
                name = "link ends time delays";
                break;
            case frequency_bands:
                name = "frequency bands";
                break;
            case doppler_integration_time:
                name = "Doppler observable integration time";
                break;
            case doppler_reference_frequency:
                name = "DSN Doppler reference frequency";
                break;
            case reception_reference_frequency_band:
                name = "DSN reference frequency band at reception";
                break;
            case sequential_range_lowest_ranging_component:
                name = "DSN sequential range lowest ranging component";
                break;
            case range_conversion_factor:
                name = "DSN range conversion factor from RU to meter";
                break;
            default:
                throw std::runtime_error(
                        "Error when getting ancillary observation data name; "
                        "variable type not recognized." );
                break;
        }

        return name;
    }

    ObservationAncillarySimulationVariable getAncillaryVariableFromString( const std::string& name )
    {
        std::map< std::string, ObservationAncillarySimulationVariable > ancillaryVariableFromStringMap = {
            { "link ends time delays", link_ends_delays },
            { "frequency bands", frequency_bands },
            { "Doppler observable integration time", doppler_integration_time },
            { "DSN Doppler reference frequency", doppler_reference_frequency },
            { "DSN reference frequency band at reception", reception_reference_frequency_band },
            { "DSN sequential range lowest ranging component", sequential_range_lowest_ranging_component },
            { "DSN range conversion factor from RU to meter", range_conversion_factor }
        };

        const auto it = ancillaryVariableFromStringMap.find( name );
        if( it == ancillaryVariableFromStringMap.end( ) )
        {
            std::string validAncillaryVariableStrings;
            for( const auto& entry : ancillaryVariableFromStringMap )
                validAncillaryVariableStrings += ( validAncillaryVariableStrings.empty( ) ? "" : ", " ) + entry.first;
            throw std::runtime_error( "Error when converting ancillary setting key " + name +
                                      ", ancillary variable not recognised. Valid options are: " + validAncillaryVariableStrings + "." );
        }
        return it->second;
    }

    void setIntermediateDoubleData( const ObservationIntermediateSimulationVariable& variableType, const double variable )
    {
        switch( variableType )
        {
            case transmitter_frequency_intermediate:
            case received_frequency_intermediate:
                doubleIntermediateData_[ variableType ] = variable;
                break;
            default:
                throw std::runtime_error(
                        "Error when setting double intermediate observation "
                        "data; could not set type " +
                        std::to_string( static_cast< int >( variableType ) ) );
        }
    }

    double getIntermediateDoubleData( const ObservationIntermediateSimulationVariable& variableType, const bool throwException = true )
    {
        double returnVariable = TUDAT_NAN;
        try
        {
            switch( variableType )
            {
                case transmitter_frequency_intermediate:
                case received_frequency_intermediate:
                    returnVariable = doubleIntermediateData_.at( variableType );
                    break;
                default:
                    if( throwException )
                    {
                        throw std::runtime_error(
                                "Error when getting double intermediate observation "
                                "data; could not retrieve type " +
                                std::to_string( static_cast< int >( variableType ) ) );
                    }
                    break;
            }
        }
        catch( ... )
        {
            if( throwException )
            {
                throw std::runtime_error(
                        "Error when getting double intermediate observation "
                        "data; could not retrieve type " +
                        std::to_string( static_cast< int >( variableType ) ) );
            }
        }
        return returnVariable;
    }

    bool operator==( const ObservationAncillarySimulationSettings& rightSettings ) const
    {
        return equals( rightSettings );
    }

    bool operator!=( const ObservationAncillarySimulationSettings& rightSettings ) const
    {
        return !( *this == rightSettings );
    }

    //! Equality comparison via equals method
    bool equals( const ObservationAncillarySimulationSettings& rhs ) const
    {
        return doubleData_ == rhs.doubleData_ && doubleVectorData_ == rhs.doubleVectorData_ &&
                doubleIntermediateData_ == rhs.doubleIntermediateData_;
    }

    std::map< ObservationAncillarySimulationVariable, double > getDoubleData( ) const
    {
        return doubleData_;
    }

    std::map< ObservationAncillarySimulationVariable, std::vector< double > > getDoubleVectorData( ) const
    {
        return doubleVectorData_;
    }

protected:
    std::map< ObservationAncillarySimulationVariable, double > doubleData_;
    std::map< ObservationAncillarySimulationVariable, std::vector< double > > doubleVectorData_;

    std::map< ObservationIntermediateSimulationVariable, double > doubleIntermediateData_;

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( doubleData_ ) );
        ar( CEREAL_NVP( doubleVectorData_ ) );
        ar( CEREAL_NVP( doubleIntermediateData_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( doubleData_ ) );
        ar( CEREAL_NVP( doubleVectorData_ ) );
        ar( CEREAL_NVP( doubleIntermediateData_ ) );
    }
};

inline std::shared_ptr< ObservationAncillarySimulationSettings > getAveragedDopplerAncillarySettings( const double integrationTime = 60.0 )
{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, integrationTime );
    return ancillarySettings;
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getNWayRangeAncillarySettings(
        const std::vector< double > linkEndsDelays = std::vector< double >( ),
        const std::vector< FrequencyBands >& frequencyBands = std::vector< FrequencyBands >( ) )
{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleVectorData( link_ends_delays, linkEndsDelays );
    ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, convertFrequencyBandsToDoubleVector( frequencyBands ) );
    return ancillarySettings;
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getNWayAveragedDopplerAncillarySettings(
        const double integrationTime = 60.0,
        const std::vector< double > linkEndsDelays = std::vector< double >( ),
        const std::vector< FrequencyBands >& frequencyBands = std::vector< FrequencyBands >( ) )
{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, integrationTime );
    ancillarySettings->setAncillaryDoubleVectorData( link_ends_delays, linkEndsDelays );
    ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, convertFrequencyBandsToDoubleVector( frequencyBands ) );
    return ancillarySettings;
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getTwoWayRangeAncillarySettings( const double retransmissionTime )
{
    return getNWayRangeAncillarySettings( std::vector< double >( { retransmissionTime } ) );
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getTwoWayAveragedDopplerAncillarySettings(
        const double integrationTime = 60.0,
        const double retransmissionTime = 0.0 )
{
    return getNWayAveragedDopplerAncillarySettings( integrationTime, std::vector< double >( { retransmissionTime } ) );
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getDsnNWayAveragedDopplerAncillarySettings(
        const std::vector< FrequencyBands >& frequencyBands,
        const FrequencyBands receptionReferenceFrequencyBand,
        const double referenceFrequency,
        const double integrationTime = 60.0,
        const std::vector< double > linkEndsDelays = std::vector< double >( ) )
{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );

    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, integrationTime );
    ancillarySettings->setAncillaryDoubleData( doppler_reference_frequency, referenceFrequency );
    ancillarySettings->setAncillaryDoubleData( reception_reference_frequency_band,
                                               convertFrequencyBandToDouble( receptionReferenceFrequencyBand ) );

    ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, convertFrequencyBandsToDoubleVector( frequencyBands ) );
    ancillarySettings->setAncillaryDoubleVectorData( link_ends_delays, linkEndsDelays );

    return ancillarySettings;
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getDsnNWayRangeAncillarySettings(
        const std::vector< FrequencyBands >& frequencyBands,
        const double lowestRangingComponent,
        const std::vector< double > linkEndsDelays = std::vector< double >( ) )

{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );

    ancillarySettings->setAncillaryDoubleData( sequential_range_lowest_ranging_component, lowestRangingComponent );

    ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, convertFrequencyBandsToDoubleVector( frequencyBands ) );
    ancillarySettings->setAncillaryDoubleVectorData( link_ends_delays, linkEndsDelays );

    return ancillarySettings;
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getDopplerMeasuredFrequencyAncillarySettings(
        const std::vector< FrequencyBands >& frequencyBands )
{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );

    ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, convertFrequencyBandsToDoubleVector( frequencyBands ) );

    return ancillarySettings;
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > getDefaultAncillaryObservationSettings(
        const observation_models::ObservableType observableType )
{
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = nullptr;
    switch( observableType )
    {
        case observation_models::one_way_differenced_range:
            ancillarySettings = getAveragedDopplerAncillarySettings( 60.0 );
            break;
        case observation_models::n_way_differenced_range:
            ancillarySettings = getAveragedDopplerAncillarySettings( 60.0 );
            break;
        default:
            break;
    }
    return ancillarySettings;
}

inline std::vector< ObservationAncillarySimulationVariable > getListRequiredAncillaryVariables(
        const observation_models::ObservableType observableType )
{
    switch( observableType )
    {
        case observation_models::one_way_differenced_range:
        case observation_models::n_way_differenced_range:
            return { doppler_integration_time };
        case observation_models::dsn_n_way_averaged_doppler:
            return { doppler_integration_time, doppler_reference_frequency, reception_reference_frequency_band, frequency_bands };
        case observation_models::dsn_n_way_range:
            return { sequential_range_lowest_ranging_component, frequency_bands };
        case observation_models::doppler_measured_frequency:
        case observation_models::one_way_doppler_measured_frequency:
            return { frequency_bands };
        default:
            return {};
    }
}

inline void checkTrackingDataAncillarySettings( const observation_models::ObservableType observableType,
                                                const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    // Extract ancillary data
    const auto doubleData = ancillarySettings->getDoubleData( );
    const auto doubleVectorData = ancillarySettings->getDoubleVectorData( );

    std::vector< ObservationAncillarySimulationVariable > requiredAncillaryVariables = getListRequiredAncillaryVariables( observableType );

    for( const ObservationAncillarySimulationVariable requiredAncillary : requiredAncillaryVariables )
    {
        bool ancillaryVariableDetected = doubleData.count( requiredAncillary ) > 0 || doubleVectorData.count( requiredAncillary ) > 0;
        if( !ancillaryVariableDetected )
        {
            throw std::runtime_error( "Error when checking ancillary settings consistency: observable '" +
                                      observation_models::getObservableName( observableType ) + "' requires ancillary setting '" +
                                      ancillarySettings->getAncillaryDataName( requiredAncillary ) + "', which was not provided." );
        }
    }
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_ANCILLARYSETTINGS_H
