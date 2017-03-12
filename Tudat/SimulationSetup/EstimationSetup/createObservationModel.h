/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include "Tudat/Astrodynamics/ObservationModels/oneWayDopplerObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/angularPositionObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/positionObservationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCalculator.h"


namespace tudat
{

namespace observation_models
{

class ObservationBiasSettings
{
public:

    ObservationBiasSettings(
            const observation_models::ObservationBiasTypes observationBiasType ):
    observationBiasType_( observationBiasType ){ }

    virtual ~ObservationBiasSettings( ){ }

    observation_models::ObservationBiasTypes observationBiasType_;
};

class ConstantObservationBiasSettings: public ObservationBiasSettings
{
public:

    ConstantObservationBiasSettings(
            const Eigen::VectorXd& observationBias ):
    ObservationBiasSettings( constant_additive_bias ), observationBias_( observationBias )
    { }

    ~ConstantObservationBiasSettings( ){ }

    Eigen::VectorXd observationBias_;

};

class ObservationSettings
{
public:


    ObservationSettings(
            const observation_models::ObservableType observableType,
            const boost::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const boost::shared_ptr< ObservationBiasSettings > biasSettings = NULL ):
        observableType_( observableType ),
        biasSettings_( biasSettings )
    {
        if( lightTimeCorrections != NULL )
        {
            lightTimeCorrectionsList_.push_back( lightTimeCorrections );
        }
    }

    ObservationSettings(
            const observation_models::ObservableType observableType,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBiasSettings > biasSettings = NULL ):
        observableType_( observableType ),lightTimeCorrectionsList_( lightTimeCorrectionsList ),
        biasSettings_( biasSettings ){ }

    ~ObservationSettings( ){ }

    observation_models::ObservableType observableType_;

    std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList_;

    boost::shared_ptr< ObservationBiasSettings > biasSettings_;
};

typedef std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSettings > > > SortedObservationSettingsMap;
typedef std::multimap< LinkEnds, boost::shared_ptr< ObservationSettings > > ObservationSettingsMap;

SortedObservationSettingsMap convertUnsortedToSortedObservationSettingsMap(
        const ObservationSettingsMap& unsortedObservationSettingsMap );

template< int ObservationSize = 1 >
boost::shared_ptr< ObservationBias< ObservationSize > > createObservationBiasCalculator(
            const LinkEnds linkEnds,
            const boost::shared_ptr< ObservationBiasSettings > biasSettings,
            const simulation_setup::NamedBodyMap &bodyMap )
{
    boost::shared_ptr< ObservationBias< ObservationSize > > observationBias;
    switch( biasSettings->observationBiasType_ )
    {
    case constant_additive_bias:
    {
        boost::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings = boost::dynamic_pointer_cast<
                ConstantObservationBiasSettings >( biasSettings );
        if( constantBiasSettings == NULL )
        {
            throw std::runtime_error( "Error when making constant observation bias, settings are inconsistent" );
        }

        if( constantBiasSettings->observationBias_.rows( ) != ObservationSize )
        {
            throw std::runtime_error( "Error when making constant observation bias, bias size is inconsistent" );
        }
        observationBias = boost::make_shared< ConstantObservationBias< ObservationSize > >(
                    constantBiasSettings->observationBias_ );
        break;
    }
    default:
    {
        std::string errorMessage = "Error when making observation bias, bias type " +
                boost::lexical_cast< std::string >( biasSettings->observationBiasType_  ) + " not recognized";
        throw std::runtime_error( errorMessage );
    }
    }
    return observationBias;
}

//! Interface class for creating observation models
/*!
 *  Interface class for creating observation models. This class is used instead of a single templated free function to
 *  allow ObservationModel deroved classed with different ObservationSize template arguments to be created using the same
 *  interface. This class has template specializations for each value of ObservationSize, and contains a single
 *  createObservationModel function that performs the required operation.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
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
    ObservationSize, ObservationScalarType, TimeType > > createObservationModel(
            const LinkEnds linkEnds,
            const boost::shared_ptr< ObservationSettings > observationSettings,
            const simulation_setup::NamedBodyMap &bodyMap );

   static boost::shared_ptr< ObservationBias< ObservationSize > > createObservationBiasCalculator(
           const LinkEnds linkEnds,
           const boost::shared_ptr< ObservationBiasSettings > observationSettings,
           const simulation_setup::NamedBodyMap &bodyMap );
};

//! Interface class for creating observation models of size 1.
template< typename ObservationScalarType, typename TimeType >
class ObservationModelCreator< 1, ObservationScalarType, TimeType >
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
    1, ObservationScalarType, TimeType > > createObservationModel(
            const LinkEnds linkEnds,
            const boost::shared_ptr< ObservationSettings > observationSettings,
            const simulation_setup::NamedBodyMap &bodyMap )
    {
        using namespace observation_models;

        boost::shared_ptr< observation_models::ObservationModel<
                1, ObservationScalarType, TimeType > > observationModel;

        // Check type of observation model.
        switch( observationSettings->observableType_ )
        {
        case one_way_range:
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

            boost::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != NULL )
            {
                observationBias =
                    createObservationBiasCalculator(
                       linkEnds, observationSettings->biasSettings_,bodyMap );
            }

            // Create observation model
            observationModel = boost::make_shared< OneWayRangeObservationModel<
                    ObservationScalarType, TimeType > >(
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                        linkEnds.at( transmitter ), linkEnds.at( receiver ),
                        bodyMap, observationSettings->lightTimeCorrectionsList_ ),
                         observationBias );

            break;
        }
        case one_way_doppler:
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

             boost::shared_ptr< ObservationBias< 1 > > observationBias;
             if( observationSettings->biasSettings_ != NULL )
             {
                 observationBias =
                     createObservationBiasCalculator(
                        linkEnds, observationSettings->biasSettings_,bodyMap );
             }

             // Create observation model
             observationModel = boost::make_shared< OneWayDopplerObservationModel<
                     ObservationScalarType, TimeType > >(
                         createLightTimeCalculator< ObservationScalarType, TimeType >(
                         linkEnds.at( transmitter ), linkEnds.at( receiver ),
                         bodyMap, observationSettings->lightTimeCorrectionsList_ ),
                         observationBias );

             break;
         }
        default:
            std::string errorMessage = "Error, observable " + boost::lexical_cast< std::string >(
                        observationSettings->observableType_ ) +
                    "  not recognized when making size 1 observation model.";
            throw std::runtime_error( errorMessage );
        }
        return observationModel;
    }

};

//! Interface class for creating observation models of size 2.
template< typename ObservationScalarType, typename TimeType >
class ObservationModelCreator< 2, ObservationScalarType, TimeType >
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
    2, ObservationScalarType, TimeType > > createObservationModel(
            const LinkEnds linkEnds,
            const boost::shared_ptr< ObservationSettings > observationSettings,
            const simulation_setup::NamedBodyMap &bodyMap )
    {
        using namespace observation_models;
        boost::shared_ptr< observation_models::ObservationModel<
                2, ObservationScalarType, TimeType > > observationModel;

        // Check type of observation model.
        switch( observationSettings->observableType_ )
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


            boost::shared_ptr< ObservationBias< 2 > > observationBias;
            if( observationSettings->biasSettings_ != NULL )
            {
                observationBias =
                    createObservationBiasCalculator< 2 >(
                       linkEnds, observationSettings->biasSettings_,bodyMap );
            }

            // Create observation model
            observationModel = boost::make_shared< AngularPositionObservationModel<
                    ObservationScalarType, TimeType > >(
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds.at( transmitter ), linkEnds.at( receiver ),
                            bodyMap, observationSettings->lightTimeCorrectionsList_ ),
                        observationBias );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + boost::lexical_cast< std::string >(
                        observationSettings->observableType_ ) +
                    "  not recognized when making size 2 observation model.";
            throw std::runtime_error( errorMessage );
            break;
        }

        return observationModel;
    }

};

//! Interface class for creating observation models of size 3.
template< typename ObservationScalarType, typename TimeType >
class ObservationModelCreator< 3, ObservationScalarType, TimeType >
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
    3, ObservationScalarType, TimeType > > createObservationModel(
            const LinkEnds linkEnds,
            const boost::shared_ptr< ObservationSettings > observationSettings,
            const simulation_setup::NamedBodyMap &bodyMap )
    {
        using namespace observation_models;
        boost::shared_ptr< observation_models::ObservationModel<
                3, ObservationScalarType, TimeType > > observationModel;

        // Check type of observation model.
        switch( observationSettings->observableType_ )
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

            if( observationSettings->lightTimeCorrectionsList_.size( ) > 0 )
            {
                throw std::runtime_error( "Error when making position observable model, found light time corrections" );
            }
            if( linkEnds.at( observed_body ).second != "" )
            {
                throw std::runtime_error( "Error, cannot yet create position function for reference point" );
            }

            boost::shared_ptr< ObservationBias< 3 > > observationBias;
            if( observationSettings->biasSettings_ != NULL )
            {
                observationBias =
                    createObservationBiasCalculator< 3 >(
                       linkEnds, observationSettings->biasSettings_,bodyMap );
            }


            // Create observation model
            observationModel = boost::make_shared< PositionObservationModel<
                    ObservationScalarType, TimeType > >(
                        boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris<
                                     ObservationScalarType, TimeType >,
                                     bodyMap.at( linkEnds.at( observed_body ).first ), _1 ),
                        observationBias );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + boost::lexical_cast< std::string >(
                        observationSettings->observableType_ ) +
                    "  not recognized when making size 3 observation model.";
            throw std::runtime_error( errorMessage );
            break;
        }
        return observationModel;
    }
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_CREATEOBSERVATIONMODEL_H
