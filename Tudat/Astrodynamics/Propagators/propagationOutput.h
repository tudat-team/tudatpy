#ifndef TUDAT_PROPAGATIONOUTPUT_H
#define TUDAT_PROPAGATIONOUTPUT_H

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace propagators
{

template< typename OutputType, typename InputType >
OutputType evaluateReferenceFunction(
        const boost::function< OutputType( const InputType&, const InputType& ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

template< typename OutputType, typename InputType >
OutputType evaluateFunction(
        const boost::function< OutputType( const InputType, const InputType ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

int getDependentVariableSize(
        const PropagationDependentVariables dependentVariableSettings );

template< typename TimeType = double, typename StateScalarType = double >
boost::function< double( ) > getDoubleDependentVariableFunction(
        const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    boost::function< double( ) > variableFunction;
    PropagationDependentVariables dependentVariable = dependentVariableSettings->variableType_;
    const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

    switch( dependentVariable )
    {
    case mach_number_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }

        boost::function< double( const double, const double ) > functionToEvaluate =
                boost::bind( &aerodynamics::computeMachNumber, _1, _2 );
        boost::function< double( ) > firstInput =
                boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeed,
                             bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        boost::function< double( ) > secondInput =
                boost::bind( &aerodynamics::FlightConditions::getCurrentSpeedOfSound,
                             bodyMap.at( bodyWithProperty )->getFlightConditions( ) );


        variableFunction = boost::bind( &evaluateFunction< double, double >,
                                        functionToEvaluate, firstInput, secondInput );
        break;
    }
    case altitude_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAltitude,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case airspeed_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeed,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case local_density_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentDensity,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case radiation_pressure_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).count( secondaryBody ) == 0 )
        {

        }
        variableFunction = boost::bind( &electro_magnetism::RadiationPressureInterface::getCurrentRadiationPressure,
                                        bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).at( secondaryBody ) );
        break;
    case relative_distance_dependent_variable:
    {

        boost::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeNormOfVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );
        break;
    }
    case relative_speed_dependent_variable:
    {
        boost::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeNormOfVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );

        break;
    }
    case single_acceleration_norm_dependent_variable:
    {
        boost::shared_ptr< SingleAccelerationNormDependentVariableSaveSettings > accelerationDependentVariableSettings =
                boost::dynamic_pointer_cast< SingleAccelerationNormDependentVariableSaveSettings >(
                    dependentVariableSettings );
        if( accelerationDependentVariableSettings == NULL )
        {
            std::string errorMessage;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                        accelerationDependentVariableSettings->associatedBody_,
                        accelerationDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, accelerationDependentVariableSettings->accelerationModeType_ );
            if( listOfSuitableAccelerationModels.size( ) != 1 )
            {
                std::string errorMessage;
                throw std::runtime_error( errorMessage );
            }
            else
            {
                boost::function< Eigen::Vector3d( ) > vectorFunction =
                        boost::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                     listOfSuitableAccelerationModels.at( 0 ) );
                variableFunction = boost::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );
            }
        }
        break;
    }
    case total_acceleration_norm_dependent_variable:
    {
        boost::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        boost::function< Eigen::Vector3d( ) > vectorFunction =
                boost::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody,
                             nBodyModel, bodyWithProperty );
        variableFunction = boost::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );

        break;
    }
    default:
        std::string errorMessage =
                "Error, did not recognize double dependent variable type when making variable function: " +
                boost::lexical_cast< std::string >( dependentVariableSettings->variableType_ );
        throw std::runtime_error( errorMessage );
    }
    return variableFunction;
}

template< typename TimeType = double, typename StateScalarType = double >
std::pair< boost::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction(
        const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    boost::function< Eigen::VectorXd( ) > variableFunction;
    int parameterSize;
    PropagationDependentVariables dependentVariable = dependentVariableSettings->variableType_;
    const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

    switch( dependentVariable )
    {
    case relative_position_dependent_variable:
    {
        boost::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        break;
    }
    case relative_velocity_dependent_variable:
    {
        boost::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );

        break;
    }
    case total_acceleration_dependent_variable:
    {
        boost::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        variableFunction =
                boost::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody, nBodyModel,
                             bodyWithProperty );

        break;
    }
    case single_acceleration_dependent_variable:
    {
        boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                boost::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationDependentVariableSettings == NULL )
        {
            std::string errorMessage;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                        accelerationDependentVariableSettings->associatedBody_,
                        accelerationDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, accelerationDependentVariableSettings->accelerationModeType_ );
            if( listOfSuitableAccelerationModels.size( ) != 1 )
            {
                std::string errorMessage;
                throw std::runtime_error( errorMessage );
            }\
            else
            {
                //boost::function< Eigen::Vector3d( ) > vectorFunction =
                variableFunction = boost::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                                listOfSuitableAccelerationModels.at( 0 ) );
                //variableFunction = boost::bind(
                //            &linear_algebra::getVectorNorm, vectorFunction );
            }
        }
        break;
    }
    default:
        std::string errorMessage =
                "Error, did not recognize vector dependent variable type when making variable function: " +
                boost::lexical_cast< std::string >( dependentVariableSettings->variableType_ );
        throw std::runtime_error( errorMessage );
    }
    return std::make_pair( variableFunction, parameterSize );
}

Eigen::VectorXd evaluateListOfFunctions(
        const std::vector< boost::function< double( ) > >& doubleFunctionList,
        const std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize );

template< typename TimeType = double, typename StateScalarType = double >
boost::function< Eigen::VectorXd( ) > createDependentVariableListFunction(
        const boost::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables =
            saveSettings->dependentVariables_;

    std::vector< boost::function< double( ) > > doubleFunctionList;
    std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList;

    int totalVariableSize = 0;
    for( unsigned int i = 0; i < dependentVariables.size( ); i++ )
    {
        if( getDependentVariableSize( dependentVariables.at( i )->variableType_ ) == 1 )
        {
            doubleFunctionList.push_back( getDoubleDependentVariableFunction(
                                              dependentVariables.at( i ),
                                              bodyMap ) );
            totalVariableSize++;
        }
        else
        {
            vectorFunctionList.push_back( getVectorDependentVariableFunction(
                                              dependentVariables.at( i ),
                                              bodyMap ) );
            totalVariableSize += vectorFunctionList.at( vectorFunctionList.size( ) - 1 ).second;
        }
    }

    return boost::bind( &evaluateListOfFunctions, doubleFunctionList, vectorFunctionList, totalVariableSize );
}

}

}
#endif // TUDAT_PROPAGATIONOUTPUT_H
