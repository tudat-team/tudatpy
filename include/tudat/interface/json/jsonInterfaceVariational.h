/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACEVARIATIONAL_H
#define TUDAT_JSONINTERFACEVARIATIONAL_H


#include "tudat/interface/json/jsonInterface.h"
#include "tudat/interface/json/estimation/parameter.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

namespace tudat
{

namespace json_interface
{

template< typename TimeType = double, typename StateScalarType = double >
class JsonVariationalEquationsSimulationManager: public JsonSimulationManager< TimeType, StateScalarType >
{
public:

    using JsonSimulationManager< TimeType, StateScalarType >::jsonObject_;
    using JsonSimulationManager< TimeType, StateScalarType >::initialClockTime_;
    using JsonSimulationManager< TimeType, StateScalarType >::dynamicsSimulator_;
    using JsonSimulationManager< TimeType, StateScalarType >::bodies_;
    using JsonSimulationManager< TimeType, StateScalarType >::integratorSettings_;
    using JsonSimulationManager< TimeType, StateScalarType >::propagatorSettings_;
    using JsonSimulationManager< TimeType, StateScalarType >::exportAsJson;
    using JsonSimulationManager< TimeType, StateScalarType >::profiling;
    using JsonSimulationManager< TimeType, StateScalarType >::exportSettingsVector_;

    JsonVariationalEquationsSimulationManager(
            const std::string& inputFilePath,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
        : JsonSimulationManager< TimeType, StateScalarType >( inputFilePath, initialClockTime ){ }

    //! Constructor from JSON object.
    /*!
     * Constructor.
     * \param jsonObject The root JSON object.
     * \param initialClockTime Initial clock time from which the cumulative CPU time during the propagation will be
     * computed. Default is the moment at which the constructor was called.
     */
    JsonVariationalEquationsSimulationManager(
            const nlohmann::json& jsonObject,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
        : JsonSimulationManager< TimeType, StateScalarType >( jsonObject, initialClockTime ){ }

    virtual ~JsonVariationalEquationsSimulationManager( ){ }

    std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > getVariationalEquationsSolver( ) const
    {
        return variationalEquationsSolver_;
    }

    void resetVariationalEquationsSolver( )
    {
        variationalEquationsSolver_ =
                std::make_shared< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodies_, integratorSettings_, propagatorSettings_, parametersToEstimate_, true,
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, false, false );
        dynamicsSimulator_ = variationalEquationsSolver_->getDynamicsSimulator( );

        if ( profiling )
        {
            std::cout << "resetVariationalEquationsSolver: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }


    void resetParameterSettings( )
    {
        updateFromJSON( parameterSettings_, jsonObject_, Keys::parametersToEstimate );
        parametersToEstimate_ = simulation_setup::createParametersToEstimate< StateScalarType >(
                    parameterSettings_, bodies_, propagatorSettings_  );

        if ( profiling )
        {

            std::cout << "resetParameterSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    virtual void runJsonSimulation( )
    {
        variationalEquationsSolver_->integrateVariationalAndDynamicalEquations(
                    propagatorSettings_->getInitialStates( ) , true );
    }

    virtual void createSimulationObjects( )
    {
        resetVariationalEquationsSolver( );
    }


    virtual void updateSettings( )
    {
        JsonSimulationManager< TimeType, StateScalarType >::updateSettings( );
        resetParameterSettings( );
        if( this->simulationType_ == variational_equations_propagation )
        {
            createSimulationObjects( );
        }
        std::cout<<"S4"<<std::endl;
    }

    virtual void exportResults( )
    {
        JsonSimulationManager< TimeType, StateScalarType >::exportResults( );
        exportResultsOfVariationalEquations( variationalEquationsSolver_, exportSettingsVector_ );
    }


protected:
    std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > variationalEquationsSolver_;

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterSettings_;

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate_;

};

//extern template class JsonVariationalEquationsSimulationManager< double, double >;


} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACEVARIATIONAL_H
