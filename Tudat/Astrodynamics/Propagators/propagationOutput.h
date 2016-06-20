#ifndef TUDAT_PROPAGATIONOUTPUT_H
#define TUDAT_PROPAGATIONOUTPUT_H

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
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

boost::function< double( ) > getDependentVariableFunction(
        const PropagationDependentVariables dependentVariable,
        const std::string& bodyWithProperty,
        const std::string& secondaryBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelList = basic_astrodynamics::AccelerationMap( ) );

Eigen::VectorXd evaluateListOfFunctions(
        const std::vector< boost::function< double( ) > >& functionList );

boost::function< Eigen::VectorXd( ) > createDependentVariableListFunction(
        const boost::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap );

}

}
#endif // TUDAT_PROPAGATIONOUTPUT_H
