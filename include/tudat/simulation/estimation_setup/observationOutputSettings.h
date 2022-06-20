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
    integration_time,
    station_local_time,
    limb_separation_angle
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

    virtual std::string getIdentifier( ) = 0;

};

enum IntegratedObservationPropertyHandling
{
    interval_start,
    interval_end,
    interval_start_and_end,
    interval_center,
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
            const bool useIncomingLink = true ):
        ObservationDependentVariableSettings( variableType ),
        relevantLinkEnd_( relevantLinkEnd ), linkEndRole_( linkEndRole ),
        integratedObservableHandling_( integratedObservableHandling ),
        useIncomingLink_( useIncomingLink )
    {

    }

    std::string getIdentifier( )
    {
        std::string identifier = "Station: (" + relevantLinkEnd_.first + ", " + relevantLinkEnd_.second + ")";
        if( linkEndRole_ != observation_models::unidentified_link_end )
        {
            throw std::runtime_error( "Error, StationAngleObservationDependentVariableSettings ID not yet implemented for link end roles" );
        }

        if( integratedObservableHandling_ != interval_undefined )
        {
            throw std::runtime_error( "Error, StationAngleObservationDependentVariableSettings ID not yet implemented for integrated obs. handling" );
        }
        return identifier;
    }

    observation_models::LinkEndId relevantLinkEnd_;

    observation_models::LinkEndType linkEndRole_;

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    bool useIncomingLink_;
};

class InterlinkObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    InterlinkObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const observation_models::LinkEndType startLinkEnd = observation_models::unidentified_link_end,
            const observation_models::LinkEndType endLinkEnd = observation_models::unidentified_link_end ):
        ObservationDependentVariableSettings( variableType ),
        startLinkEnd_( startLinkEnd ), endLinkEnd_( endLinkEnd ){ }

    observation_models::LinkEndType startLinkEnd_;

    observation_models::LinkEndType endLinkEnd_;
};

class BodyAvoidanceObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    BodyAvoidanceObservationDependentVariableSettings(
            const std::string& bodyAvoidance,
            const observation_models::LinkEndType startLinkEnd = observation_models::unidentified_link_end,
            const observation_models::LinkEndType endLinkEnd = observation_models::unidentified_link_end ):
        ObservationDependentVariableSettings( body_avoidance_angle_variable ),
        startLinkEnd_( startLinkEnd ), endLinkEnd_( endLinkEnd ), bodyAvoidance_( bodyAvoidance ){ }

    observation_models::LinkEndType startLinkEnd_;

    observation_models::LinkEndType endLinkEnd_;

    std::string bodyAvoidance_;
};



std::string getObservationDependentVariableName(
        const ObservationDependentVariables variableType );

std::string getObservationDependentVariableId(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableVectorial(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableGroundStationProperty(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool checkStationAngleVariableForGivenLink(
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings );

bool checkObservationDependentVariableForGivenLink(
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

}

}
#endif // TUDAT_OBSERVATIONOUTPUTSETTINGS
