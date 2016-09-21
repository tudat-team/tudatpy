#ifndef TUDAT_CREATENUMERICALSIMULATOR_H
#define TUDAT_CREATENUMERICALSIMULATOR_H

#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"

namespace tudat
{

namespace simulation_setup
{


//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodyMap List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param centralBodies Map of central bodies for each body undergoing acceleration.
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies );

//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodyMap List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param propagatedBodies List of bodies that are to be propagated
 *  \param centralBodies List of central bodies for each body undergoing acceleration (in same order as propagatedBodies).
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::vector< std::string >& propagatedBodies,
        const std::vector< std::string >& centralBodies );



template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > > createSingleArcDynamicsSimulator(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const boost::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated = true,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true )
{
    return boost::make_shared< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodyMap,
                integratorSettings, propagatorSettings, areEquationsOfMotionToBeIntegrated, clearNumericalSolutions,
                setIntegratedResult );
}

template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
boost::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType, ParameterType > >
createSingleArcVariationalEquationsSolver(
                    const simulation_setup::NamedBodyMap& bodyMap,
                    const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
                    const boost::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
                    const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
                    const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
                    const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
                    = boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                    const bool clearNumericalSolution = 1,
                    const bool integrateEquationsOnCreation = 1 )
{
    return boost::make_shared< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType, ParameterType > >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                integrateDynamicalAndVariationalEquationsConcurrently, variationalOnlyIntegratorSettings,
                clearNumericalSolution, integrateEquationsOnCreation );
}

}

}

#endif // CREATENUMERICALSIMULATOR_H
