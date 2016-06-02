#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace propagators
{

boost::function< double( ) > getDependentVariableFunction(
        const PropagationDependentVariables dependentVariable,
        const std::string& bodyWithProperty,
        const std::string& secondaryBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelList )
{
    boost::function< double( ) > variableFunction;
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
    case total_acceleration_norm_dependent_variable:

        break;
    }
    return variableFunction;
}


bool FixedTimeStoppingCondition::checkStopCondition( const double time  )
{
    bool stopPropagation = 0;

    if( propagationDirectionIsPositive_ && ( time >= stopTime_ ) )
    {
        stopPropagation = 1;
    }
    else if( !propagationDirectionIsPositive_ && ( time <= stopTime_ ) )
    {
        stopPropagation = 1;
    }
    return stopPropagation;
}


bool SingleVariableLimitStoppingCondition::checkStopCondition( const double time  )
{
    bool stopPropagation = 0;
    double currentVariable = variableRetrievalFuntion_( );
    if( useAsLowerBound_ && ( currentVariable < limitingValue_ ) )
    {
        stopPropagation = 1;
    }
    else if( !useAsLowerBound_ && ( currentVariable > limitingValue_ ) )
    {
        stopPropagation = 1;
    }
    return stopPropagation;
}


bool HybridStoppingCondition::checkStopCondition( const double time )
{
    if( fulFillSingleCondition_ )
    {
        bool stopPropagation = 0;
        for( unsigned int i = 0; i < stoppingCondition_.size( ); i++ )
        {
            if( stoppingCondition_.at( i )->checkStopCondition( time ) )
            {
                stopPropagation = 1;
                break;
            }
        }
        return stopPropagation;
    }
    else
    {
        bool stopPropagation = 1;
        for( unsigned int i = 0; i < stoppingCondition_.size( ); i++ )
        {
            if( !stoppingCondition_.at( i )->checkStopCondition( time ) )
            {
                stopPropagation = 0;
                break;
            }
        }
        return stopPropagation;
    }
}



boost::shared_ptr< PropagationStoppingCondition > createPropagationStoppingConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep )
{
    boost::shared_ptr< PropagationStoppingCondition > stoppingCondition;
    switch( terminationSettings->terminationType_ )
    {
    case time_stopping_condition:
    {
        boost::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
                boost::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
        stoppingCondition = boost::make_shared< FixedTimeStoppingCondition >(
                    timeTerminationSettings->terminationTime_, ( initialTimeStep > 0 ) ? true : false );
        break;
    }
    case dependent_variable_stopping_condition:
    {
        boost::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                boost::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettings );
        boost::function< double( ) > dependentVariableFunction =
                getDependentVariableFunction( dependentVariableTerminationSettings->variableType_,
                                              dependentVariableTerminationSettings->associatedBody_,
                                              dependentVariableTerminationSettings->secondaryBody_,
                                              bodyMap );
        stoppingCondition = boost::make_shared< SingleVariableLimitStoppingCondition >(
                    std::make_pair( dependentVariableTerminationSettings->variableType_,
                                    dependentVariableTerminationSettings->associatedBody_ ),
                    dependentVariableFunction, dependentVariableTerminationSettings->limitValue_,
                    dependentVariableTerminationSettings->useAsLowerLimit_ );
        break;
    }
    case hybrid_stopping_condition:
    {
        boost::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
                boost::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );

        std::vector< boost::shared_ptr< PropagationStoppingCondition > > stoppingConditionList;
        for( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); i++ )
        {
            stoppingConditionList.push_back(
                        createPropagationStoppingConditions(
                            hybridTerminationSettings->terminationSettings_.at( i ),
                            bodyMap, initialTimeStep ) );
        }
        stoppingCondition = boost::make_shared< HybridStoppingCondition >(
                    stoppingConditionList, hybridTerminationSettings->fulFillSingleCondition_ );
        break;
    }
    default:

        break;
    }
    return stoppingCondition;
}

}

}

