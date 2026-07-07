/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/basic/leastSquaresTraits.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"

// Tunable sparsity/size parameters
const int numberOfArcs = 200;
const int totalNumberOfObservations = 20000;

int main( )
{
    using namespace tudat;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::numerical_integrators;
    using namespace tudat::observation_models;
    using namespace tudat::propagators;
    using namespace tudat::simulation_setup;

    spice_interface::loadStandardSpiceKernels( );

    // Define time period
    const double startTime = 0.0;
    const double stepSize = 3600.0;
    const double finalTime = startTime + 1000.0 * physical_constants::JULIAN_DAY;
    const double arcDuration = ( finalTime - startTime ) / static_cast< double >( numberOfArcs );

    // Create bodies
    BodyListSettings bodySettings = getDefaultBodySettings( { "Sun", "Mars" }, startTime - stepSize, finalTime + stepSize );
    bodySettings.at( "Mars" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define accelerations and state vector
    const std::vector< std::string > bodiesToIntegrate = { "Mars" };
    const std::vector< std::string > centralBodies = { "SSB" };
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Mars" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    AccelerationMap accelerations = createAccelerationModelsMap( bodies, accelerationSettings, bodiesToIntegrate, centralBodies );

    // Define integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKuttaFixedStepSettings( stepSize, CoefficientSets::rungeKuttaFehlberg78 );

    // Create per-arc propagator settings
    std::vector< double > arcStartTimes;
    std::vector< double > observationTimes;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > singleArcSettings;
    for( int arc = 0; arc < numberOfArcs; ++arc )
    {
        // Define single arc propagator settings
        const double arcStart = startTime + arc * arcDuration;
        const double arcEnd = arcStart + arcDuration;
        arcStartTimes.push_back( arcStart );
        singleArcSettings.push_back( std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                centralBodies,
                accelerations,
                bodiesToIntegrate,
                getInitialStateOfBody( "Mars", "SSB", bodies, arcStart ),
                arcStart,
                integratorSettings,
                propagationTimeTerminationSettings( arcEnd ),
                cowell ) );

        // Define single arc observation settings
        const int arcObservations = totalNumberOfObservations / numberOfArcs + ( arc < totalNumberOfObservations % numberOfArcs );
        for( int i = 0; i < arcObservations; ++i )
        {
            const double fraction = ( arcObservations == 1 ) ? 0.5 : static_cast< double >( i ) / ( arcObservations - 1 );
            observationTimes.push_back( arcStart + stepSize + fraction * ( arcDuration - 2.0 * stepSize ) );
        }
    }

    // Define multi-arc propagator settings
    std::shared_ptr< MultiArcPropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double, double > >( singleArcSettings );

    // Define estimated parameters (initial state only)
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings =
            getInitialMultiArcParameterSettings( propagatorSettings, bodies, arcStartTimes );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterSettings, bodies, propagatorSettings );

    // Create obsevation model settings
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Mars", "" );
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings = {
        std::make_shared< ObservationModelSettings >( position_observable, linkEnds )
    };

    // Create estimator
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettings, propagatorSettings );

    // Simuate observations
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings = {
        std::make_shared< TabulatedObservationSimulationSettings< double > >(
                position_observable, linkEnds, observationTimes, observed_body )  };
    std::shared_ptr< ObservationCollection< double, double > > observations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );

    typedef linear_algebra::MatrixTraits< double, linear_algebra::Sparse >::matrix_type SparseDesignMatrix;

    // Perform covariance calculation with sparse estimated-parameter design-matrix storage.
    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    covarianceInput->defineCovarianceSettings( true, true, true, false );
    covarianceInput->setUseSparseDesignMatrix( true );
    std::shared_ptr< CovarianceAnalysisOutput< double, double > > output = orbitDeterminationManager.computeCovariance( covarianceInput );

    const int numberOfEstimatedParameters = parametersToEstimate->getFullParameterValues< double >( ).rows( );
    const int numberOfDesignMatrixRows = 3 * totalNumberOfObservations;
    const int numberOfDesignMatrixColumns = numberOfEstimatedParameters;
    const long long designMatrixEntries = static_cast< long long >( numberOfDesignMatrixRows ) * numberOfDesignMatrixColumns;
    const long long expectedNonZeroEntries = static_cast< long long >( numberOfDesignMatrixRows ) * 6;

    const Eigen::MatrixXd inverseCovariance = output->getNormalizedInverseCovarianceMatrix( );
    if( inverseCovariance.rows( ) != numberOfEstimatedParameters || inverseCovariance.cols( ) != numberOfEstimatedParameters )
    {
        throw std::runtime_error( "Sparse design covariance test produced an inverse covariance matrix with inconsistent dimensions." );
    }
    const SparseDesignMatrix covarianceDesignMatrix = output->getNormalizedDesignMatrix< SparseDesignMatrix >( );
    if( covarianceDesignMatrix.rows( ) != numberOfDesignMatrixRows || covarianceDesignMatrix.cols( ) != numberOfDesignMatrixColumns ||
        covarianceDesignMatrix.nonZeros( ) == designMatrixEntries )
    {
        throw std::runtime_error( "Sparse design covariance test did not store the covariance design matrix sparsely." );
    }

    std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
            observations, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ) );
    estimationInput->defineEstimationSettings( true, true, true, false, false, false );
    estimationInput->setUseSparseDesignMatrix( true );
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    const Eigen::MatrixXd estimationInverseCovariance = estimationOutput->getNormalizedInverseCovarianceMatrix( );
    if( estimationInverseCovariance.rows( ) != numberOfEstimatedParameters ||
        estimationInverseCovariance.cols( ) != numberOfEstimatedParameters )
    {
        throw std::runtime_error( "Sparse design estimation test produced an inverse covariance matrix with inconsistent dimensions." );
    }
    const SparseDesignMatrix estimationDesignMatrix = estimationOutput->getNormalizedDesignMatrix< SparseDesignMatrix >( );
    if( estimationDesignMatrix.rows( ) != numberOfDesignMatrixRows || estimationDesignMatrix.cols( ) != numberOfDesignMatrixColumns ||
        estimationDesignMatrix.nonZeros( ) == designMatrixEntries )
    {
        throw std::runtime_error( "Sparse design estimation test did not store the estimation design matrix sparsely." );
    }
    if( estimationOutput->residuals_.rows( ) != 3 * totalNumberOfObservations )
    {
        throw std::runtime_error( "Sparse design estimation test produced a residual vector with inconsistent dimensions." );
    }

    std::cout << std::setprecision( 16 ) << std::scientific;
    std::cout << "Logical design matrix rows: " << numberOfDesignMatrixRows << "\n"
              << "Logical design matrix columns: " << numberOfDesignMatrixColumns << "\n"
              << "Logical design matrix entries: " << designMatrixEntries << "\n"
              << "Expected structural nonzeros: " << expectedNonZeroEntries << "\n"
              << "Expected structural density: "
              << static_cast< double >( expectedNonZeroEntries ) / static_cast< double >( designMatrixEntries ) << "\n"
              << "Dense covariance matrix entries: " << inverseCovariance.size( ) << "\n";

}
