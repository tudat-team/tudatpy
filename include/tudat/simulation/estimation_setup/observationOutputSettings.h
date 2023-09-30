/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONOUTPUTSETTINGS
#define TUDAT_OBSERVATIONOUTPUTSETTINGS

#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace simulation_setup
{

enum ObservationDependentVariables
{
    station_elevation_angle,
    station_azimuth_angle,
    target_range,
    body_avoidance_angle_variable,
    link_body_center_distance,
    link_limb_distance,
    integration_time
};

class ObservationDependentVariableSettings
{
public:
    ObservationDependentVariableSettings(
            const ObservationDependentVariables variableType ):
        variableType_( variableType )
    {

    }

    virtual ~ObservationDependentVariableSettings( ){ }

    ObservationDependentVariables variableType_;

//    virtual std::string getIdentifier( ) = 0;

};

enum IntegratedObservationPropertyHandling
{
    interval_start,
    interval_end,
    interval_undefined
};

class StationAngleObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    StationAngleObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const observation_models::LinkEndId relevantLinkEnd,
            const observation_models::LinkEndType linkEndRole = observation_models::unidentified_link_end,
            const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined,
            const observation_models::LinkEndType originatingLinkEndRole = observation_models::unidentified_link_end ):
        ObservationDependentVariableSettings( variableType ),
        relevantLinkEnd_( relevantLinkEnd ), linkEndRole_( linkEndRole ),
        integratedObservableHandling_( integratedObservableHandling ),
        originatingLinkEndRole_( originatingLinkEndRole )
    {

    }

//    std::string getIdentifier( )
//    {
//        std::string identifier = "Station: (" + relevantLinkEnd_.bodyName_ + ", " + relevantLinkEnd_.stationName_ + ")";
//        if( linkEndRole_ != observation_models::unidentified_link_end )
//        {
//            throw std::runtime_error( "Error, StationAngleObservationDependentVariableSettings ID not yet implemented for link end roles" );
//        }
//
//        if( integratedObservableHandling_ != interval_undefined )
//        {
//            throw std::runtime_error( "Error, StationAngleObservationDependentVariableSettings ID not yet implemented for integrated obs. handling" );
//        }
//        return identifier;
//    }

    observation_models::LinkEndId relevantLinkEnd_;

    observation_models::LinkEndType linkEndRole_;

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    observation_models::LinkEndType originatingLinkEndRole_;

};

class InterlinkObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    InterlinkObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const observation_models::LinkEndType startLinkEnd,
            const observation_models::LinkEndType endLinkEnd,
            const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined,
            const std::string relativeBody = "" ):
    ObservationDependentVariableSettings( variableType ),
        startLinkEnd_( startLinkEnd ), endLinkEnd_( endLinkEnd ),
        integratedObservableHandling_( integratedObservableHandling ){ }

    ~InterlinkObservationDependentVariableSettings( ){ }

//    std::string getIdentifier( )
//    {
//        std::string identifier = "Link: (" + startLinkEnd_.bodyName_ + ", " + startLinkEnd_.stationName_ + ")";
//        if( linkEndRole_ != observation_models::unidentified_link_end )
//        {
//            throw std::runtime_error( "Error, StationAngleObservationDependentVariableSettings ID not yet implemented for link end roles" );
//        }
//
//        if( integratedObservableHandling_ != interval_undefined )
//        {
//            throw std::runtime_error( "Error, StationAngleObservationDependentVariableSettings ID not yet implemented for integrated obs. handling" );
//        }
//        return identifier;
//    }

    observation_models::LinkEndType startLinkEnd_;

    observation_models::LinkEndType endLinkEnd_;

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    std::string relativeBody_;
};



std::string getObservationDependentVariableName(
        const ObservationDependentVariables variableType );

std::string getObservationDependentVariableId(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableVectorial(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableGroundStationProperty(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableSimpleLinkProperty(
    const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

//bool doesStationAngleVariableExistForGivenLink(
//        const observation_models::ObservableType observableType,
//        const observation_models::LinkEnds& linkEnds,
//        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings );
//
//bool doesObservationDependentVariableExistForGivenLink(
//        const observation_models::ObservableType observableType,
//        const observation_models::LinkEnds& linkEnds,
//        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

}

}
#endif // TUDAT_OBSERVATIONOUTPUTSETTINGS
