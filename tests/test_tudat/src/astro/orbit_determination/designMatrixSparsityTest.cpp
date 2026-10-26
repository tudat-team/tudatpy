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
#include <vector>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/interface/spice/spiceInterface.h"
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

#include "tudat/math/basic/leastSquaresTraits.h"

#if defined( __unix__ ) || defined( __APPLE__ )
#include <sys/resource.h>
#endif

//! Peak resident memory of this process so far, in MB (zero where not available)
double getPeakMemoryUsageInMegabytes( )
{
#if defined( __unix__ ) || defined( __APPLE__ )
    struct rusage usage;
    getrusage( RUSAGE_SELF, &usage );
#if defined( __APPLE__ )
    return static_cast< double >( usage.ru_maxrss ) / ( 1024.0 * 1024.0 );  // bytes
#else
    return static_cast< double >( usage.ru_maxrss ) / 1024.0;  // kilobytes
#endif
#else
    return 0.0;
#endif
}

using DesignMatrixStorageClass = tudat::linear_algebra::Sparse;
using DesignMatrixPrecision = double;

using DesignMatrixTraits = tudat::linear_algebra::MatrixTraits<DesignMatrixPrecision, DesignMatrixStorageClass>;

// Tunable sparsity/size parameters
const int numberOfArcs = 100;
const int totalNumberOfObservations = 10000;

// Set to 4 to test proper convergence behaviour
const int numberOfIterations = 1;

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

    // Perturb parameter estimation
    Eigen::VectorXd initialEstimate = parametersToEstimate->getFullParameterValues< double >( );
    for( int i = 0; i < numberOfArcs; ++i )
    {
        initialEstimate.segment( 6 * i, 3 ).array( ) += 1.0;
        initialEstimate.segment( 6 * i + 3, 3 ).array( ) += 1.0E-5;
    }
    parametersToEstimate->resetParameterValues( initialEstimate );

    // Perform estimation, starting from the perturbed initial estimate, with sparse or dense linear algebra.
    // Returns the estimation output, and the parameter values after the final correction has been applied.
    auto runEstimation = [ & ]( const bool useSparseDesignMatrix ) {
        parametersToEstimate->resetParameterValues( initialEstimate );
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( observations );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );
        estimationInput->defineEstimationSettings( true, true, true, true );
        estimationInput->setSparseDesignMatrix( useSparseDesignMatrix );
        std::shared_ptr< EstimationOutput< double, double > > output = orbitDeterminationManager.estimateParameters( estimationInput );
        return std::make_pair( output, Eigen::VectorXd( parametersToEstimate->getFullParameterValues< double >( ) ) );
    };

    // The sparse run is done first, so that the peak memory after it is that of the sparse run only
    const double peakMemoryBeforeEstimation = getPeakMemoryUsageInMegabytes( );
    auto [ sparseOutput, sparseCorrectedParameters ] = runEstimation( true );
    const double peakMemoryAfterSparseEstimation = getPeakMemoryUsageInMegabytes( );
    auto [ denseOutput, denseCorrectedParameters ] = runEstimation( false );
    const double peakMemoryAfterDenseEstimation = getPeakMemoryUsageInMegabytes( );
    std::cout << "Peak memory (MB) before estimation: " << peakMemoryBeforeEstimation
              << ", after sparse estimation: " << peakMemoryAfterSparseEstimation
              << ", after dense estimation: " << peakMemoryAfterDenseEstimation << "\n";

    // Test sparsity
    const DesignMatrixTraits::matrix_type designMatrix =
        sparseOutput->getUnnormalizedDesignMatrixGeneric<DesignMatrixStorageClass>( );
    const double tolerance = 1.0E-20;
    const int nonZeroCount = ( designMatrix.coeffs( ).abs( ) > tolerance ).count( );
    std::cout << std::setprecision( 16 ) << std::scientific;
    std::cout << "Design matrix entries: " << designMatrix.rows( ) * designMatrix.cols( ) << "\n"
              << "Nonzeros > " << tolerance << ": " << nonZeroCount << "\n"
              << "Density: " << static_cast< double >( nonZeroCount ) / static_cast< double >( designMatrix.size( ) ) << "\n";

    // Compare sparse and dense solutions. The corrected parameters are compared, since the parameter estimate in the output
    // is the one from the best iteration, before its correction was applied.
    const Eigen::VectorXd parameterDifference = sparseCorrectedParameters - denseCorrectedParameters;
    const double maximumParameterCorrection = ( denseCorrectedParameters - initialEstimate ).cwiseAbs( ).maxCoeff( );
    const Eigen::VectorXd relativeFormalErrorDifference =
            ( sparseOutput->getFormalErrorVector( ) - denseOutput->getFormalErrorVector( ) )
                    .cwiseQuotient( denseOutput->getFormalErrorVector( ) );
    const double maximumDesignMatrixDifference =
            ( Eigen::MatrixXd( designMatrix ) - denseOutput->getUnnormalizedDesignMatrix( ) ).cwiseAbs( ).maxCoeff( );
    std::cout << "Max. parameter correction (dense): " << maximumParameterCorrection << "\n"
              << "Max. parameter difference (sparse - dense): " << parameterDifference.cwiseAbs( ).maxCoeff( ) << "\n"
              << "Max. relative formal error difference: " << relativeFormalErrorDifference.cwiseAbs( ).maxCoeff( ) << "\n"
              << "Max. design matrix difference: " << maximumDesignMatrixDifference << "\n";

    // Position components are perturbed by 1 m, velocities by 1E-5 m/s; the solutions should agree far below that
    // Check that each run stored the design matrix in the storage it was asked to use
    const bool storageAsRequested =
            std::holds_alternative< Eigen::SparseMatrix< double > >( sparseOutput->normalizedDesignMatrix_ ) &&
            std::holds_alternative< Eigen::MatrixXd >( denseOutput->normalizedDesignMatrix_ );
    std::cout << "Design matrix stored as requested (sparse/dense): " << std::boolalpha << storageAsRequested << "\n";

    const bool solutionsAgree = storageAsRequested && maximumParameterCorrection > 1.0E-3 &&
            parameterDifference.cwiseAbs( ).maxCoeff( ) < 1.0E-6 &&
            relativeFormalErrorDifference.cwiseAbs( ).maxCoeff( ) < 1.0E-8 && maximumDesignMatrixDifference == 0.0;
    std::cout << ( solutionsAgree ? "Sparse and dense solutions agree" : "Sparse and dense solutions DIFFER" ) << std::endl;
    return solutionsAgree ? EXIT_SUCCESS : EXIT_FAILURE;
}
