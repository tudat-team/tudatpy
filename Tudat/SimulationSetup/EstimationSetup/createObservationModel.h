/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOBSERVATIONMODEL_H
#define TUDAT_CREATEOBSERVATIONMODEL_H

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>


#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/angularPositionObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/positionObservationModel.h"
#include "Tudat/SimulationSetup/body.h"


namespace tudat
{

namespace observation_models
{

//! Interface class for creating observation models
/*!
 *  Interface class for creating observation models. This class is used instead of a single templated free function to
 *  allow ObservationModel deroved classed with different ObservationSize template arguments to be created using the same
 *  interface. This class has template specializations for each value of ObservationSize, and contains a single
 *  createObservationModel function that performs the required operation.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double,
          typename TimeType = double, typename StateScalarType = ObservationScalarType >
class ObservationModelCreator
{
public:

    //! Function to create an observation model.
    /*!
     * Function to create an observation model.
     * \param observableType Type of observation model that is to be created.
     * \param linkEnds Link ends for observation model that is to be created
     * \param bodyMap List of body objects that comprises the environment
     * \param singleObservableCorrections List of light time corrections that are used when computing the observable.
     * \param observationBiasCalculator
     * \return Observation model of required settings.
     */
    static boost::shared_ptr< observation_models::ObservationModel<
    ObservationSize, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBias< ObservationSize > > observationBiasCalculator = NULL );
};

//! Interface class for creating observation models of size 1.
template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModelCreator< 1, ObservationScalarType, TimeType, StateScalarType >
{
public:

    //! Function to create an observation model of size 1.
    /*!
     * Function to create an observation model of size 1.
     * \param observableType Type of observation model that is to be created.
     * \param linkEnds Link ends for observation model that is to be created
     * \param bodyMap List of body objects that comprises the environment
     * \param singleObservableCorrections List of light time corrections that are used when computing the observable.
     * \param observationBiasCalculator
     * \return Observation model of required settings.
     */
    static boost::shared_ptr< observation_models::ObservationModel<
    1, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL )
    {
        using namespace observation_models;

        boost::shared_ptr< observation_models::ObservationModel<
                1, ObservationScalarType, TimeType, StateScalarType > > observationModel;

        // Check type of observation model.
        switch( observableType )
        {
        case oneWayRange:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage =
                        "Error when making 1 way range model, " +
                        boost::lexical_cast< std::string >( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way range model, no receiver found" );
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way range model, no transmitter found" );
            }

            // Create observation model
            observationModel = boost::make_shared< OneWayRangeObservationModel<
                    ObservationScalarType, TimeType, StateScalarType > >(
                        createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                        linkEnds.at( transmitter ), linkEnds.at( receiver ),
                        bodyMap, singleObservableCorrections ), observationBiasCalculator );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + boost::lexical_cast< std::string >( observableType ) +
                    "  not recognized when making size 1 observation model.";
            throw std::runtime_error( errorMessage );
        }
        return observationModel;
    }
};

//! Interface class for creating observation models of size 2.
template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModelCreator< 2, ObservationScalarType, TimeType, StateScalarType >
{
public:

    //! Function to create an observation model of size 2.
    /*!
     * Function to create an observation model of size 2.
     * \param observableType Type of observation model that is to be created.
     * \param linkEnds Link ends for observation model that is to be created
     * \param bodyMap List of body objects that comprises the environment
     * \param singleObservableCorrections List of light time corrections that are used when computing the observable.
     * \param observationBiasCalculator
     * \return Observation model of required settings.
     */
    static boost::shared_ptr< observation_models::ObservationModel<
    2, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = NULL )
    {
        using namespace observation_models;
        boost::shared_ptr< observation_models::ObservationModel<
                2, ObservationScalarType, TimeType, StateScalarType > > observationModel;

        // Check type of observation model.
        switch( observableType )
        {
        case angular_position:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage =
                        "Error when making angular position model, " +
                        boost::lexical_cast< std::string >( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making angular position model, no receiver found" );
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making angular position model, no transmitter found" );
            }

            // Create observation model
            observationModel = boost::make_shared< AngularPositionObservationModel<
                    ObservationScalarType, TimeType, StateScalarType > >(
                        createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                            linkEnds.at( transmitter ), linkEnds.at( receiver ),
                            bodyMap, singleObservableCorrections ), observationBiasCalculator );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + boost::lexical_cast< std::string >( observableType ) +
                    "  not recognized when making size 2 observation model.";
            throw std::runtime_error( errorMessage );
            break;
        }

        return observationModel;
    }
};

//! Interface class for creating observation models of size 3.
template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModelCreator< 3, ObservationScalarType, TimeType, StateScalarType >
{
public:

    //! Function to create an observation model of size 3.
    /*!
     * Function to create an observation model of size 3.
     * \param observableType Type of observation model that is to be created.
     * \param linkEnds Link ends for observation model that is to be created
     * \param bodyMap List of body objects that comprises the environment
     * \param singleObservableCorrections List of light time corrections that are used when computing the observable.
     * \param observationBiasCalculator
     * \return Observation model of required settings.
     */
    static boost::shared_ptr< observation_models::ObservationModel<
    3, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
           const boost::shared_ptr< ObservationBias < 3 > > observationBiasCalculator = NULL )
    {
        using namespace observation_models;
        boost::shared_ptr< observation_models::ObservationModel<
                3, ObservationScalarType, TimeType, StateScalarType > > observationModel;

        // Check type of observation model.
        switch( observableType )
        {
        case position_observable:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 1 )
            {
                std::string errorMessage =
                        "Error when making position observable model, " +
                        boost::lexical_cast< std::string >( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }

            if( linkEnds.count( observed_body ) == 0 )
            {
                throw std::runtime_error( "Error when making position observable model, no observed_body found" );
            }

            if( singleObservableCorrections.size( ) > 0 )
            {
                throw std::runtime_error( "Error when making position observable model, found light time corrections" );
            }
            if( linkEnds.at( observed_body ).second != "" )
            {
                throw std::runtime_error( "Error, cannot yet create position function for reference point" );
            }

            // Create observation model
            observationModel = boost::make_shared< PositionObservationModel<
                    ObservationScalarType, TimeType, StateScalarType > >(
                        boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris<
                                     StateScalarType, TimeType >,
                                     bodyMap.at( linkEnds.at( observed_body ).first ), _1 ), observationBiasCalculator );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + boost::lexical_cast< std::string >( observableType ) +
                    "  not recognized when making size 3 observation model.";
            throw std::runtime_error( errorMessage );
            break;
        }
        return observationModel;
    }
};

}

}

#endif // TUDAT_CREATEOBSERVATIONMODEL_H
