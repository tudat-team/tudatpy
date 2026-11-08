// Restricted C++ prototype: independent Mars arcs, isolated environments,
// shared mu column, serialized access to the SPICE-backed input models.
#include "isolatedArcExecutor.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using Clock = std::chrono::steady_clock;
double seconds( Clock::time_point start )
{
    return std::chrono::duration< double >( Clock::now( ) - start ).count( );
}
std::mutex spiceMutex;

class LockedRotation : public ephemerides::RotationalEphemeris
{
    std::shared_ptr< ephemerides::RotationalEphemeris > source_;

public:
    explicit LockedRotation( std::shared_ptr< ephemerides::RotationalEphemeris > source ):
        RotationalEphemeris( source->getBaseFrameOrientation( ), source->getTargetFrameOrientation( ) ), source_( source )
    {}
    Eigen::Quaterniond getRotationToBaseFrame( double t ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex );
        return source_->getRotationToBaseFrame( t );
    }
    Eigen::Quaterniond getRotationToTargetFrame( double t ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex );
        return source_->getRotationToTargetFrame( t );
    }
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( double t ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex );
        return source_->getDerivativeOfRotationToBaseFrame( t );
    }
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( double t ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex );
        return source_->getDerivativeOfRotationToTargetFrame( t );
    }
    void getFullRotationalQuantitiesToTargetFrame( Eigen::Quaterniond& q, Eigen::Matrix3d& d, Eigen::Vector3d& w, double t ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex );
        source_->getFullRotationalQuantitiesToTargetFrame( q, d, w, t );
    }
};

class LockedEphemeris : public ephemerides::Ephemeris
{
    std::shared_ptr< ephemerides::Ephemeris > source_;

public:
    explicit LockedEphemeris( std::shared_ptr< ephemerides::Ephemeris > source ):
        Ephemeris( source->getReferenceFrameOrigin( ), source->getReferenceFrameOrientation( ) ), source_( source )
    {}
    Eigen::Vector6d getCartesianState( double t ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex );
        return source_->getCartesianState( t );
    }
};

using Settings = SingleArcPropagatorSettings< double, double >;
using Solver = SingleArcVariationalEquationsSolver< double, double >;
struct Context {
    SystemOfBodies bodies;
    std::shared_ptr< Settings > settings;
    std::shared_ptr< Solver > solver;
};

SystemOfBodies makeBodies( int degree, bool locked )
{
    auto settings = getDefaultBodySettings( { "Mars" }, "Mars", "J2000" );
    settings.at( "Mars" )->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d, degree );
    auto bodies = createSystemOfBodies< double, double >( settings );
    if( locked )
    {
        auto mars = bodies.at( "Mars" );
        mars->setRotationalEphemeris( std::make_shared< LockedRotation >( mars->getRotationalEphemeris( ) ) );
        mars->setEphemeris( std::make_shared< LockedEphemeris >( mars->getEphemeris( ) ) );
    }
    bodies.createEmptyBody( "Satellite" );
    return bodies;
}

std::shared_ptr< Settings > makeSettings( SystemOfBodies& bodies, int arc, int degree, double days, double step )
{
    const double epoch = 778118400.0 + arc * ( days * 86400.0 + 3600.0 );
    const auto gravity =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( bodies.at( "Mars" )->getGravityFieldModel( ) );
    Eigen::Vector6d elements;
    elements << gravity->getReferenceRadius( ) + 300000.0, 0.001, 93.0 * mathematical_constants::PI / 180.0, 0.2, 0.3,
            std::fmod( 0.4 + 0.7 * arc, 2.0 * mathematical_constants::PI );
    Eigen::Vector6d state =
            orbital_element_conversions::convertKeplerianToCartesianElements( elements, gravity->getGravitationalParameter( ) );
    SelectedAccelerationMap selected;
    selected[ "Satellite" ][ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( degree, degree ) );
    auto accelerations = createAccelerationModelsMap( bodies, selected, { "Satellite" }, { "Mars" } );
    return std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
            std::vector< std::string >{ "Mars" },
            accelerations,
            std::vector< std::string >{ "Satellite" },
            state,
            epoch,
            numerical_integrators::rungeKuttaFixedStepSettings< double >( step, numerical_integrators::rungeKutta87DormandPrince ),
            std::make_shared< PropagationTimeTerminationSettings >( epoch + days * 86400.0, true ) );
}

std::vector< Context > makeContexts( int arcs, int degree, double days, double step )
{
    std::vector< Context > contexts;
    // All model construction and file/kernel access happens before workers run.
    for( int arc = 0; arc < arcs; ++arc )
    {
        Context context;
        context.bodies = makeBodies( degree, true );
        context.settings = makeSettings( context.bodies, arc, degree, days, step );
        auto parameters = getInitialStateParameterSettings< double, double >( context.settings, context.bodies );
        parameters.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                "Mars", estimatable_parameters::gravitational_parameter ) );
        context.solver = std::make_shared< Solver >(
                context.bodies, context.settings, createParametersToEstimate< double, double >( parameters, context.bodies ), true, false );
        contexts.push_back( std::move( context ) );
    }
    return contexts;
}

// Scattering gives every local STM its own global initial-state columns, while
// every local sensitivity maps to the same final global mu column.
Eigen::MatrixXd combine( const Eigen::MatrixXd& stm, const Eigen::MatrixXd& sensitivity, int arc, int arcs )
{
    if( stm.rows( ) != 6 || stm.cols( ) != 6 || sensitivity.rows( ) != 6 || sensitivity.cols( ) != 1 )
    {
        throw std::runtime_error( "unexpected local matrix dimensions" );
    }
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero( 6, 6 * arcs + 1 );
    result.block( 0, 6 * arc, 6, 6 ) = stm;
    result.rightCols( 1 ) = sensitivity;
    return result;
}

void checkMatrix( const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected, const std::string& label )
{
    if( actual.rows( ) != expected.rows( ) || actual.cols( ) != expected.cols( ) || !actual.allFinite( ) || !expected.allFinite( ) )
    {
        throw std::runtime_error( label + ": invalid dimensions or values" );
    }
    if( !( ( actual - expected ).array( ).abs( ) <= ( 1.0e-14 + 1.0e-11 * expected.array( ).abs( ) ) ).all( ) )
    {
        throw std::runtime_error( label + ": numerical mismatch" );
    }
}

void testExecutor( )
{
    std::vector< int > counts( 9, 0 );
    executeIsolatedArcs( 9, 4, [ & ]( std::size_t arc ) { ++counts.at( arc ); } );
    for( int count : counts )
    {
        if( count != 1 )
        {
            throw std::runtime_error( "task coverage failed" );
        }
    }
    try
    {
        executeIsolatedArcs( 9, 4, []( std::size_t arc ) {
            if( arc == 3 )
            {
                throw std::runtime_error( "injected failure" );
            }
        } );
        throw std::logic_error( "missing task exception" );
    }
    catch( const std::runtime_error& error )
    {
        if( std::string( error.what( ) ) != "arc 3: injected failure" )
        {
            throw;
        }
    }
}

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
