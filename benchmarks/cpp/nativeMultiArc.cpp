// Integration/acceptance test for the native independent-arc C++ API.
#include "marsArcBenchmark.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"

using MultiSolver = MultiArcVariationalEquationsSolver< double, double >;
struct NativeRun {
    SystemOfBodies bodies;
    std::shared_ptr< MultiArcPropagatorSettings< double, double > > settings;
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameters;
    std::shared_ptr< MultiSolver > solver;
};

NativeRun makeNative( int arcs, int degree, double days, bool linked = false )
{
    NativeRun result;
    result.bodies = makeBodies( degree, false );
    result.bodies.at( "Satellite" )
            ->setEphemeris( std::make_shared< ephemerides::MultiArcEphemeris >(
                    std::map< double, std::shared_ptr< ephemerides::Ephemeris > >{
                            { 0.0, std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "Mars", "J2000" ) } },
                    "Mars",
                    "J2000" ) );
    result.bodies.processBodyFrameDefinitions< double, double >( );
    std::vector< std::shared_ptr< Settings > > settings;
    std::vector< double > starts;
    for( int arc = 0; arc < arcs; ++arc )
    {
        auto setting = makeSettings( result.bodies, arc, degree, days, 30.0 );
        if( linked && arc > 0 )
        {
            setting->resetInitialTime( 778118400.0 + arc * days * 43200.0 );
        }
        setting->resetDependentVariablesToSave( { relativePositionDependentVariable( "Satellite", "Mars" ) } );
        starts.push_back( setting->getInitialTime( ) );
        settings.push_back( setting );
    }
    result.settings = std::make_shared< MultiArcPropagatorSettings< double, double > >( settings, linked );
    result.settings->getOutputSettings( )->setIntegratedResult( true );
    result.settings->getOutputSettings( )->setUpdateDependentVariableInterpolator( true );
    auto parameters = getInitialMultiArcParameterSettings< double, double >( result.settings, result.bodies, starts );
    parameters.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
            "Mars", estimatable_parameters::gravitational_parameter ) );
    result.parameters = createParametersToEstimate< double, double >( parameters, result.bodies );
    result.solver = std::make_shared< MultiSolver >( result.bodies, result.settings, result.parameters, false );
    return result;
}

MultiSolver::IsolatedArcVariationalSolverFactory makeFactory( int degree, double days )
{
    return [ = ]( std::size_t arc, const std::shared_ptr< Settings >& original, const Eigen::VectorXd& values ) {
        auto bodies = makeBodies( degree, true );
        auto settings = makeSettings( bodies, arc, degree, days, 30.0 );
        settings->resetDependentVariablesToSave( original->getDependentVariablesToSave( ) );
        auto parameters = getInitialStateParameterSettings< double, double >( settings, bodies );
        parameters.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                "Mars", estimatable_parameters::gravitational_parameter ) );
        auto set = createParametersToEstimate< double, double >( parameters, bodies );
        set->resetParameterValues( values );
        return std::make_shared< Solver >( bodies, settings, set, true, false );
    };
}

void integrate( NativeRun& run )
{
    run.solver->integrateVariationalAndDynamicalEquations( run.settings->getInitialStateList( ), true );
    if( !run.solver->getDynamicsSimulator( )->integrationCompletedSuccessfully( ) )
    {
        throw std::runtime_error( "native integration incomplete" );
    }
}

std::size_t compare( NativeRun& actual, NativeRun& reference )
{
    auto first = actual.solver->getDynamicsSimulator( )->getSingleArcDynamicsSimulators( );
    auto second = reference.solver->getDynamicsSimulator( )->getSingleArcDynamicsSimulators( );
    const auto a = actual.solver->getNumericalVariationalEquationsSolution( );
    const auto b = reference.solver->getNumericalVariationalEquationsSolution( );
    std::size_t epochs = 0;
    if( first.size( ) != second.size( ) )
    {
        throw std::runtime_error( "arc count mismatch" );
    }
    for( std::size_t arc = 0; arc < first.size( ); ++arc )
    {
        const auto& states = first.at( arc )->getEquationsOfMotionNumericalSolution( );
        if( states.size( ) != second.at( arc )->getEquationsOfMotionNumericalSolution( ).size( ) ||
            a.at( arc ).at( 0 ).size( ) != b.at( arc ).at( 0 ).size( ) || a.at( arc ).at( 1 ).size( ) != b.at( arc ).at( 1 ).size( ) )
        {
            throw std::runtime_error( "native history length mismatch" );
        }
        for( const auto& entry : states )
        {
            const auto epoch = entry.first;
            checkMatrix( entry.second, second.at( arc )->getEquationsOfMotionNumericalSolution( ).at( epoch ), "native state" );
            for( int kind = 0; kind < 2; ++kind )
            {
                checkMatrix( a.at( arc ).at( kind ).at( epoch ), b.at( arc ).at( kind ).at( epoch ), "native variational history" );
            }
            checkMatrix( first.at( arc )->getDependentVariableHistory( ).at( epoch ),
                         second.at( arc )->getDependentVariableHistory( ).at( epoch ),
                         "dependent variable history" );
            ++epochs;
        }
        const double epoch = states.rbegin( )->first;
        const double middle = 0.5 * ( states.begin( )->first + epoch );
        checkMatrix( actual.solver->getDynamicsSimulator( )
                             ->getMultiArcPropagationResults( )
                             ->getMultiArcDependentVariablesInterface( )
                             ->getDependentVariables( middle ),
                     reference.solver->getDynamicsSimulator( )
                             ->getMultiArcPropagationResults( )
                             ->getMultiArcDependentVariablesInterface( )
                             ->getDependentVariables( middle ),
                     "dependent variable interface" );
        if( first.at( arc )->getSingleArcPropagationResults( )->getTotalNumberOfFunctionEvaluations( ) !=
            second.at( arc )->getSingleArcPropagationResults( )->getTotalNumberOfFunctionEvaluations( ) )
        {
            throw std::runtime_error( "function evaluation metadata mismatch" );
        }
        checkMatrix( actual.solver->getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( epoch ),
                     reference.solver->getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( epoch ),
                     "native combined interface" );
        if( !actual.bodies.at( "Satellite" )->getEphemeris( ) || !reference.bodies.at( "Satellite" )->getEphemeris( ) )
        {
            throw std::runtime_error( "requested integrated ephemeris is missing" );
        }
        checkMatrix( actual.bodies.at( "Satellite" )->getEphemeris( )->getCartesianState( epoch ),
                     reference.bodies.at( "Satellite" )->getEphemeris( )->getCartesianState( epoch ),
                     "updated global ephemeris" );
        if( first.at( arc )->getPropagationTerminationReason( )->getPropagationTerminationReason( ) !=
            second.at( arc )->getPropagationTerminationReason( )->getPropagationTerminationReason( ) )
        {
            throw std::runtime_error( "termination metadata mismatch" );
        }
    }
    return epochs;
}

template< typename Function >
void expectFailure( Function function, const std::string& text )
{
    bool failed = false;
    try
    {
        function( );
    }
    catch( const std::exception& error )
    {
        failed = std::string( error.what( ) ).find( text ) != std::string::npos;
    }
    if( !failed )
    {
        throw std::runtime_error( "expected failure containing: " + text );
    }
}

void regressionChecks( )
{
    // Run before any propagation: concurrently initialize/copy the DP87 table.
    // The old getter tested higherOrder != 7 while initializing it to 8, so it
    // rewrote a shared static table on every call, even after a serial warm-up.
    executeIsolatedArcs( 8, 8, []( std::size_t ) {
        for( int iteration = 0; iteration < 1000; ++iteration )
        {
            const auto coefficients =
                    numerical_integrators::RungeKuttaCoefficients::get( numerical_integrators::rungeKutta87DormandPrince );
            if( coefficients.higherOrder != 8 || coefficients.lowerOrder != 7 || coefficients.aCoefficients.rows( ) != 13 ||
                coefficients.aCoefficients.cols( ) != 12 || !coefficients.aCoefficients.allFinite( ) ||
                coefficients.aCoefficients( 1, 0 ) != 1.0 / 18.0 )
            {
                throw std::runtime_error( "concurrent RK coefficient initialization/copy failed" );
            }
        }
    } );
    const auto factory = makeFactory( 20, 0.02 );
    auto original = makeNative( 2, 20, 0.02 );
    auto isolated = makeNative( 2, 20, 0.02 );
    integrate( original );
    auto results = isolated.solver->getVariationalPropagationResults( );
    auto normalResults = isolated.solver->getDynamicsSimulator( )->getMultiArcPropagationResults( );
    auto interface = isolated.solver->getStateTransitionMatrixInterface( );
    isolated.solver->setIndependentArcExecution( 4, factory );
    integrate( isolated );
    compare( isolated, original );
    if( results != isolated.solver->getVariationalPropagationResults( ) ||
        normalResults != isolated.solver->getDynamicsSimulator( )->getMultiArcPropagationResults( ) ||
        interface != isolated.solver->getStateTransitionMatrixInterface( ) )
    {
        throw std::runtime_error( "native result/interface identity changed" );
    }
    Eigen::VectorXd parameters = original.parameters->getFullParameterValues< double >( );
    parameters( 0 ) += 10.0;
    parameters( parameters.size( ) - 1 ) *= 1.00001;
    original.solver->resetParameterEstimate( parameters, true );
    isolated.solver->resetParameterEstimate( parameters, true );
    compare( isolated, original );
    isolated.solver->clearIndependentArcExecution( );
    integrate( isolated );
    compare( isolated, original );
    if( isolated.solver->getDynamicsSimulator( )->getLastIntegrationUsedIndependentArcs( ) )
    {
        throw std::runtime_error( "clear did not restore sequential execution" );
    }
    expectFailure( [ & ] { isolated.solver->setIndependentArcExecution( 0, factory ); }, "positive workers" );
    // A shallow-copy/aliased body factory must be rejected before threads start.
    isolated.solver->getDynamicsSimulator( )->setIndependentArcExecution( 2, [ & ]( std::size_t arc, const Eigen::MatrixXd& ) {
        return isolated.solver->getDynamicsSimulator( )->getSingleArcDynamicsSimulators( ).at( arc );
    } );
    expectFailure( [ & ] { integrate( isolated ); }, "shared simulator" );
    // Wrong current shared-parameter values must not be silently propagated.
    isolated.solver->setIndependentArcExecution( 2, [ & ]( std::size_t arc, const auto& settings, const auto& values ) {
        Eigen::VectorXd wrong = values;
        wrong( wrong.size( ) - 1 ) *= 1.01;
        return factory( arc, settings, wrong );
    } );
    expectFailure( [ & ] { integrate( isolated ); }, "parameter values" );
    isolated.solver->setIndependentArcExecution( 2, [ & ]( std::size_t arc, const auto& settings, const auto& values ) {
        auto solver = factory( arc, settings, values );
        if( arc == 1 )
        {
            solver->getDynamicsSimulator( )->getDynamicsStateDerivative( )->setStateDerivativeModifierFunction(
                    []( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) -> Eigen::MatrixXd {
                        throw std::runtime_error( "injected worker failure" );
                    } );
        }
        return solver;
    } );
    expectFailure( [ & ] { integrate( isolated ); }, "arc 1" );
    isolated.solver->setIndependentArcExecution( 2, factory );
    integrate( isolated );
    compare( isolated, original );
    // Linked arcs must never invoke the factory; they use the existing state-transfer path.
    auto linked = makeNative( 2, 20, 0.02, true );
    auto linkedReference = makeNative( 2, 20, 0.02, true );
    linked.solver->setIndependentArcExecution( 4, []( std::size_t, const auto&, const auto& ) -> std::shared_ptr< Solver > {
        throw std::runtime_error( "factory must not run for dependent arcs" );
    } );
    integrate( linkedReference );
    integrate( linked );
    compare( linked, linkedReference );
    if( linked.solver->getDynamicsSimulator( )->getLastIntegrationUsedIndependentArcs( ) )
    {
        throw std::runtime_error( "dependent arcs used independent execution" );
    }
    // The first propagation resolves transfer NaNs. Reintegration must still
    // retain the conservative sequential policy for this linked configuration.
    integrate( linkedReference );
    integrate( linked );
    compare( linked, linkedReference );
    // The same native executor must also support dynamics-only output containers.
    auto dynamics = makeNative( 2, 20, 0.02 );
    auto dynamicsReference = makeNative( 2, 20, 0.02 );
    dynamics.solver->setIndependentArcExecution( 2, factory );
    dynamics.solver->getDynamicsSimulator( )->integrateEquationsOfMotion( dynamics.settings->getInitialStateList( ) );
    dynamicsReference.solver->getDynamicsSimulator( )->integrateEquationsOfMotion( dynamicsReference.settings->getInitialStateList( ) );
    const auto dynamicStates = dynamics.solver->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    const auto referenceStates = dynamicsReference.solver->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    for( std::size_t arc = 0; arc < dynamicStates.size( ); ++arc )
    {
        if( dynamicStates.at( arc ).size( ) != referenceStates.at( arc ).size( ) )
        {
            throw std::runtime_error( "dynamics-only history length mismatch" );
        }
        for( const auto& state : dynamicStates.at( arc ) )
        {
            checkMatrix( state.second, referenceStates.at( arc ).at( state.first ), "dynamics-only state" );
        }
    }
    std::cout << "{\"regression_checks\":\"passed\"}\n";
}

int main( int argc, char** argv )
{
    try
    {
        int arcs = 8, degree = 50;
        double days = 2.0;
        bool checksOnly = false;
        for( int i = 1; i < argc; ++i )
        {
            const std::string option = argv[ i ];
            if( option == "--checks-only" )
            {
                checksOnly = true;
                continue;
            }
            if( i + 1 == argc )
            {
                throw std::invalid_argument( "missing value" );
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
            else
            {
                throw std::invalid_argument( "unknown option" );
            }
        }
        if( arcs < 1 || degree < 2 || degree > 120 || !std::isfinite( days ) || days * 86400 < 240 )
        {
            throw std::invalid_argument( "invalid workload" );
        }
        spice_interface::loadStandardSpiceKernels( );
        regressionChecks( );
        if( checksOnly )
        {
            return 0;
        }
        auto reference = makeNative( arcs, degree, days );
        integrate( reference );
        for( int workers : { 1, 2, 4, 8 } )
        {
            const auto start = Clock::now( );
            auto run = makeNative( arcs, degree, days );
            run.solver->setIndependentArcExecution( workers, makeFactory( degree, days ) );
            const auto integrationStart = Clock::now( );
            integrate( run );
            const double integrationSeconds = seconds( integrationStart );
            const double totalSeconds = seconds( start );
            const auto epochs = compare( run, reference );
            if( !run.solver->getDynamicsSimulator( )->getLastIntegrationUsedIndependentArcs( ) )
            {
                throw std::runtime_error( "independent path was not used" );
            }
            std::cout << std::setprecision( 12 ) << "{\"workers\":" << workers << ",\"native_integration_seconds\":" << integrationSeconds
                      << ",\"total_seconds\":" << totalSeconds << ",\"checked_epochs\":" << epochs << ",\"arcs\":" << arcs
                      << ",\"degree\":" << degree << ",\"days_per_arc\":" << days << ",\"validation\":\"passed\"}" << std::endl;
        }
    }
    catch( const std::exception& error )
    {
        std::cerr << error.what( ) << '\n';
        return 1;
    }
}
