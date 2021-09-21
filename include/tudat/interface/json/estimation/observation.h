/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_OBSERVATION_H
#define TUDAT_JSONINTERFACE_OBSERVATION_H

#include <boost/lexical_cast.hpp>

#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace observation_models
{

//! Map of `ObservableType` string representations.
static std::map< ObservableType, std::string > observationTypes =
{
    { one_way_range, "oneWayRange" },
    { angular_position, "angularPosition" },
    { position_observable, "positionObservable" },
    { one_way_doppler, "oneWayDoppler" },
    { one_way_differenced_range, "oneWayDifferencedRange" },
    { n_way_range, "nWayRange" },
    { two_way_doppler, "twoWayDoppler" }
};

//! Map of `ObservableType` string representations.
static std::map< std::string, ObservableType > observationTypesInverse =
{
    { "oneWayRange", one_way_range },
    { "angularPosition", angular_position },
    { "positionObservable", position_observable },
    { "oneWayDoppler", one_way_doppler },
    { "oneWayDifferencedRange", one_way_differenced_range },
    { "nWayRange", n_way_range },
    { "twoWayDoppler", two_way_doppler }
};

//! Convert `ObservableType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const ObservableType& observableType )
{
    jsonObject = json_interface::stringFromEnum( observableType, observationTypes );
}

//! Convert `json` to `ObservableType`.
inline void from_json( const nlohmann::json& jsonObject, ObservableType& observableType )
{
    observableType = json_interface::enumFromString( jsonObject, observationTypes );
}

//! Map of `ObservableType` string representations.
static std::map< ObservationBiasTypes, std::string > observationBiasTypes =
{

    { multiple_observation_biases, "multipleObservationBiases" },
    { constant_absolute_bias, "constantAbsoluteBias" },
    { constant_relative_bias, "constantRelativeBias" },
    { arc_wise_constant_absolute_bias, "arcWiseConstantAbsoluteBias" },
    { arc_wise_constant_relative_bias, "arcWiseConstantRelativeBias" }
};

//! Convert `ObservableType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const ObservationBiasTypes& observableBiasType )
{
    jsonObject = json_interface::stringFromEnum( observableBiasType, observationBiasTypes );
}

//! Convert `json` to `ObservableType`.
inline void from_json( const nlohmann::json& jsonObject, ObservationBiasTypes& observableBiasType )
{
    observableBiasType = json_interface::enumFromString( jsonObject, observationBiasTypes );
}


//! Map of `LinkEndType` string representations.
static std::map< LinkEndType, std::string > linkEndTypes =
{
    { transmitter, "transmitter" },
    { reflector1, "reflector1" },
    { reflector, "reflector" },
    { reflector2, "reflector2" },
    { reflector3, "reflector3" },
    { reflector4, "reflector4" },
    { receiver, "receiver" },
    { observed_body, "observedBody" }
};

//! Map of `LinkEndType` string representations.
static std::map< std::string, LinkEndType > linkEndTypesInverse =
{
    { "transmitter", transmitter,  },
    { "reflector1", reflector1 },
    { "reflector", reflector },
    { "reflector2", reflector2 },
    { "reflector3", reflector3 },
    { "reflector4", reflector4 },
    { "receiver", receiver },
    { "observedBody", observed_body }
};

//! Convert `LinkEndType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const LinkEndType& observableType )
{
    jsonObject = json_interface::stringFromEnum( observableType, linkEndTypes );
}

//! Convert `json` to `LinkEndType`.
inline void from_json( const nlohmann::json& jsonObject, LinkEndType& observableType )
{
    observableType = json_interface::enumFromString( jsonObject, linkEndTypes );
}


//! Map of `LinkEndType` string representations.
static std::map< DopplerProperTimeRateType, std::string > dopplerProperTimeRateTypes =
{
    { direct_first_order_doppler_proper_time_rate, "firsOrderProperTimeRate" }
};

//! Convert `LinkEndType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const DopplerProperTimeRateType& rateType )
{
    jsonObject = json_interface::stringFromEnum( rateType, dopplerProperTimeRateTypes );
}

//! Convert `json` to `LinkEndType`.
inline void from_json( const nlohmann::json& jsonObject, DopplerProperTimeRateType& rateType )
{
    rateType = json_interface::enumFromString( jsonObject, dopplerProperTimeRateTypes );
}

//! Map of `LinkEndType` string representations.
static std::map< LightTimeCorrectionType, std::string > lichtTimeCorrectionTypes =
{
    { first_order_relativistic, "firsOrderRelativistic" }
};


//! Convert `LinkEndType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const LightTimeCorrectionType& correctionType )
{
    jsonObject = json_interface::stringFromEnum( correctionType, lichtTimeCorrectionTypes );
}

//! Convert `json` to `LinkEndType`.
inline void from_json( const nlohmann::json& jsonObject, LightTimeCorrectionType& correctionType )
{
    correctionType = json_interface::enumFromString( jsonObject, lichtTimeCorrectionTypes );
}


//! Create a `json` object from a shared pointer to a `ObservationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ObservationModelSettings >& parameterSettings );

//! Create a shared pointer to a `ObservationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ObservationModelSettings >& parameterSettings );




//! Create a `json` object from a shared pointer to a `ObservationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ObservationBiasSettings >& parameterSettings );

//! Create a shared pointer to a `ObservationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ObservationBiasSettings >& parameterSettings );



//! Create a `json` object from a shared pointer to a `ObservationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< DopplerProperTimeRateSettings >& properTimeRateSettings );

//! Create a shared pointer to a `ObservationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< DopplerProperTimeRateSettings >& properTimeRateSettings );



//! Create a `json` object from a shared pointer to a `ObservationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< LightTimeCorrectionSettings >& lightTimeCorrectionSettings );

//! Create a shared pointer to a `ObservationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< LightTimeCorrectionSettings >& lightTimeCorrectionSettings );


//! Create a `json` object from a shared pointer to a `ObservationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ObservationSimulationSettings< double > >& ObservationSimulationSettings );

//! Create a shared pointer to a `ObservationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ObservationSimulationSettings< double > >& ObservationSimulationSettings );


//! Create a `json` object from a shared pointer to a `ObservationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ObservationViabilitySettings >& observationViabilitySettings );

//! Create a shared pointer to a `ObservationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ObservationViabilitySettings >& observationViabilitySettings );



} // namespace propagators

} // namespace tudat

namespace boost
{

template<>
tudat::observation_models::ObservableType lexical_cast( const std::string & s );


template<>
std::string lexical_cast(const tudat::observation_models::ObservableType & component );


template<>
tudat::observation_models::LinkEndType lexical_cast( const std::string & s );


template<>
std::string lexical_cast(const tudat::observation_models::LinkEndType & component );


template<>
tudat::observation_models::LinkEnds lexical_cast( const std::string & s );

template<>
std::string lexical_cast(const tudat::observation_models::LinkEnds & component );

//}

}


#endif // TUDAT_JSONINTERFACE_OBSERVATION_H
