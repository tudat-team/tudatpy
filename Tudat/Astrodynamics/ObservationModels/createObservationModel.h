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
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/createLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/angularPositionObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/positionObservationModel.h"
#include "Tudat/SimulationSetup/body.h"


namespace tudat
{

namespace observation_models
{

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class ObservationModelCreator
{
public:
    static boost::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBias< ObservationSize > > observationBiasCalculator = NULL );
};

template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModelCreator< 1, ObservationScalarType, TimeType, StateScalarType >
{
public:
    static boost::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL )
    {
        using namespace observation_models;

        boost::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType > > observationModel;

        switch( observableType )
        {
        case oneWayRange:
        {
            if( linkEnds.size( ) != 2 )
            {
                std::cerr<<"Error when making 1 way range model, "<< linkEnds.size( ) <<" link ends found"<<std::endl;
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                std::cerr<<"Error when making 1 way range model, no receiver found"<<std::endl;
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                std::cerr<<"Error when making 1 way range model, no transmitter found"<<std::endl;
            }
            observationModel = boost::make_shared< OneWayRangeObservationModel<
                    ObservationScalarType, TimeType, StateScalarType > >(
                        createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                        linkEnds.at( transmitter ), linkEnds.at( receiver ),
                        bodyMap, singleObservableCorrections ), observationBiasCalculator );

            break;
        }
        default:
            std::cerr<<"Error, observable "<<observableType<<" not recognized when making size 1 observation model."<<std::endl;
            break;
        }
        return observationModel;
    }
};

template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModelCreator< 2, ObservationScalarType, TimeType, StateScalarType >
{
public:
    static boost::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = NULL )
    {
        using namespace observation_models;
        boost::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType > > observationModel;

        switch( observableType )
        {
        case angular_position:
        {
            if( linkEnds.size( ) != 2 )
            {
                std::cerr<<"Error when making angular position model, "<< linkEnds.size( ) <<" link ends found"<<std::endl;
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                std::cerr<<"Error when making angular position model, no receiver found"<<std::endl;
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                std::cerr<<"Error when making angular position model, no transmitter found"<<std::endl;
            }
            observationModel = boost::make_shared< AngularPositionObservationModel<
                    ObservationScalarType, TimeType, StateScalarType > >(
                        createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                        linkEnds.at( transmitter ), linkEnds.at( receiver ),
                        bodyMap, singleObservableCorrections ), observationBiasCalculator );

            break;
        }
        default:
            std::cerr<<"Error, observable "<<observableType<<" not recognized when making size 2 observation model."<<std::endl;
            break;
        }

        return observationModel;
    }
};

template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModelCreator< 3, ObservationScalarType, TimeType, StateScalarType >
{
public:
    static boost::shared_ptr< observation_models::ObservationModel< 3, ObservationScalarType, TimeType, StateScalarType > > createObservationModel(
            const ObservableType observableType,
            const LinkEnds& linkEnds,
            const simulation_setup::NamedBodyMap &bodyMap,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& singleObservableCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
           const boost::shared_ptr< ObservationBias < 3 > > observationBiasCalculator = NULL )
    {
        using namespace observation_models;
        boost::shared_ptr< observation_models::ObservationModel< 3, ObservationScalarType, TimeType, StateScalarType > > observationModel;
        switch( observableType )
        {
        case position_observable:
        {
            if( linkEnds.size( ) != 1 )
            {
                std::cerr<<"Error position observable model, "<< linkEnds.size( ) <<" link ends found"<<std::endl;
            }

            if( linkEnds.count( observed_body ) == 0 )
            {
                std::cerr<<"Error when making position observable model, no observed_body found"<<std::endl;
            }

            if( singleObservableCorrections.size( ) > 0 )
            {
                std::cerr<<"Error when making position observable model, found light time corrections"<<std::endl;
            }

            if( linkEnds.at( observed_body ).second != "" )
            {
                std::cerr<<"Error, cannot yet create position function for reference point"<<std::endl;
            }
            observationModel = boost::make_shared< PositionObservationModel<
                    ObservationScalarType, TimeType, StateScalarType > >(
                        boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris<
                                     StateScalarType, TimeType >,
                                     bodyMap.at( linkEnds.at( observed_body ).first ), _1 ), observationBiasCalculator );

            break;
        }
        default:
            std::cerr<<"Error, observable "<<observableType<<" not recognized when making size 3 observation model."<<std::endl;
            break;
        }
        return observationModel;
    }
};

}

}

#endif // TUDAT_CREATEOBSERVATIONMODEL_H
