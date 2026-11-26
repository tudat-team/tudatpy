#include "marsArcBenchmark.h"

int main( int argc, char** argv )
{
    try
    {
        int arcs = 8, degree = 50;
        double days = 2.0, step = 30.0;
        std::vector< int > workers{ 1, 2, 4, 8 };
        for( int i = 1; i < argc; ++i )
        {
            const std::string option = argv[ i ];
            if( option == "--self-test" )
            {
                testExecutor( );
                std::cout << "executor checks passed\n";
                return 0;
            }
            if( option == "--help" )
            {
                std::cout << "--arcs N --degree L --days-per-arc D --step S --workers 1,2,4,8\n";
                return 0;
            }
            if( i + 1 == argc )
            {
                throw std::invalid_argument( "missing option value" );
            }
            const std::string value = argv[ ++i ];
            if( option == "--arcs" )
            {
                arcs = std::stoi( value );
            }
            else if( option == "--degree" )
            {
                degree = std::stoi( value );
            }
            else if( option == "--days-per-arc" )
            {
                days = std::stod( value );
            }
            else if( option == "--step" )
            {
                step = std::stod( value );
            }
            else if( option == "--workers" )
            {
                workers.clear( );
                std::istringstream input( value );
                std::string item;
                while( std::getline( input, item, ',' ) )
                {
                    workers.push_back( std::stoi( item ) );
                }
            }
            else
            {
                throw std::invalid_argument( "unknown option " + option );
            }
        }
        if( arcs < 1 || degree < 2 || degree > 120 || !std::isfinite( days ) || !std::isfinite( step ) || days <= 0 || step <= 0 ||
            days * 86400 < step * 8 || workers.empty( ) )
        {
            throw std::invalid_argument( "invalid benchmark configuration" );
        }
        for( int count : workers )
        {
            if( count < 1 )
            {
                throw std::invalid_argument( "invalid workers" );
            }
        }
        testExecutor( );
        spice_interface::loadStandardSpiceKernels( );
        const auto referenceStart = Clock::now( );
        auto referenceBodies = makeBodies( degree, false );
        std::vector< std::shared_ptr< Settings > > settings;
        std::vector< double > starts;
        for( int arc = 0; arc < arcs; ++arc )
        {
            settings.push_back( makeSettings( referenceBodies, arc, degree, days, step ) );
            starts.push_back( settings.back( )->getInitialTime( ) );
        }
        auto multiSettings = std::make_shared< MultiArcPropagatorSettings< double, double > >( settings, false );
        auto parameterSettings = getInitialMultiArcParameterSettings< double, double >( multiSettings, referenceBodies, starts );
        parameterSettings.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                "Mars", estimatable_parameters::gravitational_parameter ) );
        MultiArcVariationalEquationsSolver< double, double > reference(
                referenceBodies, multiSettings, createParametersToEstimate< double, double >( parameterSettings, referenceBodies ), true );
        const double referenceSeconds = seconds( referenceStart );
        const auto referenceMatrices = reference.getNumericalVariationalEquationsSolution( );
        const auto referenceStates = reference.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
        std::cout << std::setprecision( 12 ) << "{\"reference_seconds\":" << referenceSeconds << ",\"arcs\":" << arcs
                  << ",\"degree\":" << degree << ",\"days_per_arc\":" << days << "}" << std::endl;
        for( int count : workers )
        {
            const auto totalStart = Clock::now( );
            auto contexts = makeContexts( arcs, degree, days, step );
            const double setupSeconds = seconds( totalStart );
            const auto integrationStart = Clock::now( );
            executeIsolatedArcs( arcs, count, [ & ]( std::size_t arc ) {
                auto& context = contexts.at( arc );
                context.solver->integrateVariationalAndDynamicalEquations( context.settings->getInitialStates( ), true );
                if( !context.solver->getDynamicsSimulator( )->integrationCompletedSuccessfully( ) )
                {
                    throw std::runtime_error( "integration did not complete" );
                }
            } );
            const double integrationSeconds = seconds( integrationStart );
            const auto mergeStart = Clock::now( );
            std::size_t epochs = 0;
            // Materialize the full global matrix history as a prototype merge.
            // Keep arc histories separate (no map-key collisions at boundaries).
            std::vector< std::map< double, Eigen::MatrixXd > > combined( arcs );
            for( int arc = 0; arc < arcs; ++arc )
            {
                auto& solver = contexts.at( arc ).solver;
                for( const auto& entry : solver->getStateTransitionMatrixSolution( ) )
                {
                    combined.at( arc ).emplace(
                            entry.first, combine( entry.second, solver->getSensitivityMatrixSolution( ).at( entry.first ), arc, arcs ) );
                    ++epochs;
                }
            }
            const double mergeSeconds = seconds( mergeStart );
            const double totalSeconds = seconds( totalStart );
            // Validation is excluded from timed work, and checks every saved
            // epoch, including the actual native global interface at endpoints.
            for( int arc = 0; arc < arcs; ++arc )
            {
                auto& solver = contexts.at( arc ).solver;
                const auto& states = solver->getEquationsOfMotionSolution( );
                if( states.size( ) != referenceStates.at( arc ).size( ) ||
                    combined.at( arc ).size( ) != referenceMatrices.at( arc ).at( 0 ).size( ) )
                {
                    throw std::runtime_error( "history length mismatch" );
                }
                for( const auto& entry : states )
                {
                    checkMatrix( entry.second, referenceStates.at( arc ).at( entry.first ), "state" );
                }
                for( const auto& entry : combined.at( arc ) )
                {
                    checkMatrix( entry.second,
                                 combine( referenceMatrices.at( arc ).at( 0 ).at( entry.first ),
                                          referenceMatrices.at( arc ).at( 1 ).at( entry.first ),
                                          arc,
                                          arcs ),
                                 "merged variational matrix" );
                }
                const double end = states.rbegin( )->first;
                if( std::abs( end - ( starts.at( arc ) + days * 86400.0 ) ) > 1.0e-6 )
                {
                    throw std::runtime_error( "incorrect arc end" );
                }
                checkMatrix( combined.at( arc ).rbegin( )->second,
                             reference.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( end ),
                             "native global interface" );
            }
            std::cout << "{\"workers\":" << count << ",\"setup_seconds\":" << setupSeconds
                      << ",\"integration_seconds\":" << integrationSeconds << ",\"merge_seconds\":" << mergeSeconds
                      << ",\"total_seconds\":" << totalSeconds << ",\"checked_epochs\":" << epochs << ",\"validation\":\"passed\"}"
                      << std::endl;
        }
    }
    catch( const std::exception& error )
    {
        std::cerr << error.what( ) << std::endl;
        return 1;
    }
}
