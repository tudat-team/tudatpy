#ifndef TUDAT_PROPAGATIONOUTPUTSETTINGS_H
#define TUDAT_PROPAGATIONOUTPUTSETTINGS_H

#include <string>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace propagators
{


//! Enum listing the dependent variables that can be saved during the propagation
enum PropagationDependentVariables
{
    mach_number_dependent_variable = 0,
    altitude_dependent_variable = 1,
    airspeed_dependent_variable = 2,
    local_density_dependent_variable = 3,
    relative_speed_dependent_variable = 4,
    relative_position_dependent_variable = 5,
    relative_distance_dependent_variable = 6,
    relative_velocity_dependent_variable = 7,
    radiation_pressure_dependent_variable = 8,
    total_acceleration_norm_dependent_variable = 9,
    single_acceleration_norm_dependent_variable = 10,
    total_acceleration_dependent_variable = 11,
    single_acceleration_dependent_variable = 12
};


class SingleDependentVariableSaveSettings
{
public:
    SingleDependentVariableSaveSettings(
            const PropagationDependentVariables variableType,
            const std::string& associatedBody,
            const std::string& secondaryBody = "" ):
        variableType_( variableType ), associatedBody_( associatedBody ), secondaryBody_( secondaryBody ){ }

    virtual ~SingleDependentVariableSaveSettings( ){ }

    PropagationDependentVariables variableType_;

    std::string associatedBody_;

    std::string secondaryBody_;

};

class SingleAccelerationNormDependentVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    SingleAccelerationNormDependentVariableSaveSettings(
            const basic_astrodynamics::AvailableAcceleration accelerationModeType,
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration ):
        SingleDependentVariableSaveSettings(
            single_acceleration_norm_dependent_variable, bodyUndergoingAcceleration, bodyExertingAcceleration ),
        accelerationModeType_( accelerationModeType ){ }

    basic_astrodynamics::AvailableAcceleration accelerationModeType_;

};

class SingleAccelerationDependentVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    SingleAccelerationDependentVariableSaveSettings(
            const basic_astrodynamics::AvailableAcceleration accelerationModeType,
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration ):
        SingleDependentVariableSaveSettings(
            single_acceleration_dependent_variable, bodyUndergoingAcceleration, bodyExertingAcceleration ),
        accelerationModeType_( accelerationModeType ){ }

    basic_astrodynamics::AvailableAcceleration accelerationModeType_;

};

class DependentVariableSaveSettings
{
public:
    DependentVariableSaveSettings(
            const std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables ):
    dependentVariables_( dependentVariables ){ }

    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables_;
};

}

}

#endif // TUDAT_PROPAGATIONOUTPUTSETTINGS_H
