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
    body_avoidance_angle,
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
        ObservationDependentVariableSettings( body_avoidance_angle ),
        startLinkEnd_( startLinkEnd ), endLinkEnd_( endLinkEnd ), bodyAvoidance_( bodyAvoidance ){ }

    observation_models::LinkEndType startLinkEnd_;

    observation_models::LinkEndType endLinkEnd_;

    std::string bodyAvoidance_;
};


std::string getObservationDependentVariableName(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableGroundStationProperty(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

void checkObservationDependentVariableEnvironment(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

typedef std::function< double( const std::vector< double >& ,
                               const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                               const Eigen::VectorXd& ) > DoubleObservationDependentVariableFunction;

typedef std::function< Eigen::VectorXd( const std::vector< double >& ,
                               const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                               const Eigen::VectorXd& ) > VectorObservationDependentVariableFunction;

DoubleObservationDependentVariableFunction getBodyAvoidanceFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< BodyAvoidanceObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

DoubleObservationDependentVariableFunction getTargetRangeFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

DoubleObservationDependentVariableFunction getStationObservationAngleFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

DoubleObservationDependentVariableFunction getObservationDoubleDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );


}

}
#endif // TUDAT_OBSERVATIONOUTPUTSETTINGS
