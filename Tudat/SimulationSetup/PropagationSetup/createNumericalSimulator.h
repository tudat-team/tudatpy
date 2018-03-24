/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_CREATENUMERICALSIMULATOR_H
#define TUDAT_CREATENUMERICALSIMULATOR_H

#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"

namespace tudat
{

namespace simulation_setup
{

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
boost::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >
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
    return boost::make_shared< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                integrateDynamicalAndVariationalEquationsConcurrently, variationalOnlyIntegratorSettings,
                clearNumericalSolution, integrateEquationsOnCreation );
}

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATENUMERICALSIMULATOR_H
