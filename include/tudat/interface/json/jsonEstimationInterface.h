/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONESTIMATIONINTERFACE_H
#define TUDAT_JSONESTIMATIONINTERFACE_H

#include "tudat/interface/json/unitTestSupport.h"
#include "tudat/interface/json/propagation/acceleration.h"
#include "tudat/interface/json/estimation/observation.h"
#include "tudat/interface/json/estimation/parameter.h"
#include "tudat/interface/json/estimation/orbitDetermination.h"
#include "tudat/interface/json/support/valueConversions.h"
#include "tudat/interface/json/jsonInterfaceVariational.h"

#include "tudat/simulation/estimation.h"


namespace tudat
{

namespace json_interface
{



//! Class for managing JSON-based simulations.
template< typename TimeType = double, typename StateScalarType = double >
class JsonEstimationManager: public JsonVariationalEquationsSimulationManager< TimeType, StateScalarType >
{

public:

    using JsonSimulationManager< TimeType, StateScalarType >::jsonObject_;
    using JsonSimulationManager< TimeType, StateScalarType >::applicationOptions_;
    using JsonSimulationManager< TimeType, StateScalarType >::profiling;
    using JsonSimulationManager< TimeType, StateScalarType >::initialClockTime_;
    using JsonSimulationManager< TimeType, StateScalarType >::inputFilePath_;
    using JsonSimulationManager< TimeType, StateScalarType >::dynamicsSimulator_;
    using JsonSimulationManager< TimeType, StateScalarType >::bodies_;
    using JsonSimulationManager< TimeType, StateScalarType >::integratorSettings_;
    using JsonSimulationManager< TimeType, StateScalarType >::propagatorSettings_;
    using JsonSimulationManager< TimeType, StateScalarType >::exportAsJson;

    using JsonVariationalEquationsSimulationManager< TimeType, StateScalarType >::variationalEquationsSolver_;
    using JsonVariationalEquationsSimulationManager< TimeType, StateScalarType >::parametersToEstimate_;

    //! Constructor from JSON file.
    /*!
     * Constructor.
     * \param inputFilePath Path to the root JSON input file. Can be absolute or relative (to the working directory).
     * \param initialClockTime Initial clock time from which the cummulative CPU time during the propagation will be
     * computed. Default is the moment at which the constructor was called.
     */
    JsonEstimationManager(
            const std::string& inputFilePath,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) ):
        JsonVariationalEquationsSimulationManager< TimeType, StateScalarType >( inputFilePath, initialClockTime )
    { }

    //! Constructor from JSON object.
    /*!
     * Constructor.
     * \param jsonObject The root JSON object.
     * \param initialClockTime Initial clock time from which the cummulative CPU time during the propagation will be
     * computed. Default is the moment at which the constructor was called.
     */
    JsonEstimationManager(
            const nlohmann::json& jsonObject,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) ):
        JsonVariationalEquationsSimulationManager< TimeType, StateScalarType >( jsonObject, initialClockTime ){ }

    virtual ~JsonEstimationManager( ){ }

    virtual void updateSettings( )
    {
        std::cout<<"T1"<<std::endl;
        JsonVariationalEquationsSimulationManager< TimeType, StateScalarType >::updateSettings( );
        std::cout<<"T2"<<std::endl;

        resetObservationSettings( );
        std::cout<<"T3"<<std::endl;

        resetEstimationSettings( );
        std::cout<<"T4"<<std::endl;
        createSimulationObjects( );
        std::cout<<"T5"<<std::endl;

    }

    virtual void runJsonSimulation( )
    {
        estimationOutput_ = orbitDeterminationManager_->estimateParameters( estimationInput_, convergenceChecker_ );
    }


    //! Export the results of the dynamics simulation according to the export settings.
    /*!
     * @copybrief exportResults
     */
    void exportEstimationResults( )
    {

    }

protected:

    void resetObservationSettings( )
    {
        observationSettingsMap_ = getValue< observation_models::ObservationSettingsListPerLinkEnd >(
            jsonObject_, Keys::observations );
        std::cout<<"Size of obs. settings "<<observationSettingsMap_.size( )<<std::endl;

        if ( profiling )
        {

            std::cout << "resetObservationSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    void resetEstimationSettings( )
    {
        updatePodSettingsFromJSON(
                    jsonObject_[ Keys::estimationSettings ], estimationInput_, convergenceChecker_, parametersToEstimate_->getParameterSetSize( )  );

        if ( profiling )
        {
            std::cout << "resetEstimationSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }


    //! Reset dynamicsSimulator_ for the current bodies_, integratorSettings_ and propagatorSettings_.
    /*!
     * @copybrief resetDynamicsSimulator
     */
    virtual void resetDynamicsSimulator( )
    {
        this->createSimulationObjects( );
    }

    virtual void resetVariationalEquationsSolver( )
    {
        this->createSimulationObjects( );
    }

    virtual void createSimulationObjects( )
    {
        orbitDeterminationManager_ =
                std::make_shared< simulation_setup::OrbitDeterminationManager< StateScalarType, TimeType > >(
                    bodies_, parametersToEstimate_, observation_models::convertUnsortedToSortedObservationSettingsMap(
                        observationSettingsMap_ ), integratorSettings_, propagatorSettings_,
                    false );
        variationalEquationsSolver_ =
                std::dynamic_pointer_cast< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    orbitDeterminationManager_->getVariationalEquationsSolver( ) );
        dynamicsSimulator_ = variationalEquationsSolver_->getDynamicsSimulator( );
        std::cout<<"check NULL "<<( variationalEquationsSolver_ == NULL )<<" "<<
                   ( orbitDeterminationManager_->getVariationalEquationsSolver( ) == NULL )<<std::endl;
    }

private:

    observation_models::ObservationSettingsListPerLinkEnd observationSettingsMap_;

    std::shared_ptr< simulation_setup::OrbitDeterminationManager< StateScalarType, TimeType > > orbitDeterminationManager_;


    std::shared_ptr< simulation_setup::EstimationConvergenceChecker > convergenceChecker_;

    std::shared_ptr< simulation_setup::EstimationInput< StateScalarType, TimeType > > estimationInput_;


    std::shared_ptr< simulation_setup::EstimationOutput< StateScalarType, TimeType > > estimationOutput_;

};

//extern template class JsonEstimationManager< double, double >;

//#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//extern template class JsonEstimationManager< Time, long double >;
//extern template class JsonEstimationManager< double, double >;
//extern template class JsonEstimationManager< Time, long double >;
//#endif


} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONESTIMATIONINTERFACE_H
