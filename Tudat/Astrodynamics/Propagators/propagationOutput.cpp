#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Propagators/propagationOutput.h"

namespace tudat
{

namespace propagators
{

Eigen::VectorXd evaluateListOfFunctions(
        const std::vector< boost::function< double( ) > >& functionList  )
{
    Eigen::VectorXd variableList = Eigen::VectorXd::Zero( functionList.size( ) );
    for( unsigned int i = 0; i < functionList.size( ); i++ )
    {
        variableList( i ) = functionList.at( i )( );
    }

    return variableList;
}

boost::function< Eigen::VectorXd( ) > createDependentVariableListFunction(
        const boost::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    std::map< PropagationDependentVariables, std::vector< std::pair< std::string, std::string > > > dependentVariables =
            saveSettings->dependentVariables_;

    std::vector< boost::function< double( ) > > doubleFunctionList;
    std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList;

    for( std::map< PropagationDependentVariables, std::vector< std::pair< std::string, std::string > > >::const_iterator variableIterator =
         dependentVariables.begin( ); variableIterator != dependentVariables.end( ); variableIterator++ )
    {
        for( unsigned int i = 0; i < variableIterator->second.size( ); i++ )
        {
            if( true )
            {
            doubleFunctionList.push_back( getDependentVariableFunction(
                                        variableIterator->first, variableIterator->second.at( i ).first,
                                        variableIterator->second.at( i ).second, bodyMap ) );
            }
        }
    }

    return boost::bind( &evaluateListOfFunctions, doubleFunctionList );
}

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

}

}
