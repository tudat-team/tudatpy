// Second inner-loop candidate: acceleration summation and Legendre cache.
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include <chrono>
#include <iostream>
#include <iomanip>
using namespace tudat;
using namespace simulation_setup;
using namespace propagators;
using Clock = std::chrono::steady_clock;
using Solver = SingleArcVariationalEquationsSolver< double, Time >;
using Dynamics = SingleArcDynamicsSimulator< double, Time >;
class ExperimentalAcceleration final : public gravitation::SphericalHarmonicsGravitationalAccelerationModel
{
    std::vector< Eigen::Vector3d > termBuffer_;
    bool reuseTermBuffer_ = false;

public:
    using SphericalHarmonicsGravitationalAccelerationModel::SphericalHarmonicsGravitationalAccelerationModel;
    void setReuseTermBuffer( bool reuse )
    {
        reuseTermBuffer_ = reuse;
    }
    void updateMembers( double time = TUDAT_NAN ) override
    {
        updateMembersWithTermBuffer( time, reuseTermBuffer_ ? &termBuffer_ : nullptr );
    }
};
struct Run {
    std::shared_ptr< Solver > solver;
    std::shared_ptr< Dynamics > dynamics;
    Eigen::Vector6d state;
};
Run makeRun( int degree, double duration, bool variational, bool original = false, bool reuseTermBuffer = false )
{
    auto bs = getDefaultBodySettings( { "Mars" }, "Mars", "J2000" );
    bs.at( "Mars" )->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d, degree );
    auto bodies = createSystemOfBodies< double, Time >( bs );
    bodies.createEmptyBody< double, Time >( "Satellite" );
    bodies.processBodyFrameDefinitions< double, Time >( );
    auto gravity = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( bodies.at( "Mars" )->getGravityFieldModel( ) );
    Eigen::Vector6d elements;
    elements << gravity->getReferenceRadius( ) + 300000., .001, 93. * mathematical_constants::PI / 180., .2, .3, .4;
    auto state = orbital_element_conversions::convertKeplerianToCartesianElements( elements, gravity->getGravitationalParameter( ) );
    basic_astrodynamics::AccelerationMap acceleration;
    if( original )
    {
        SelectedAccelerationMap selected;
        selected[ "Satellite" ][ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( degree, degree ) );
        acceleration = createAccelerationModelsMap( bodies, selected, { "Satellite" }, { "Mars" } );
    }
    else
    {
        auto satellite = bodies.at( "Satellite" ), mars = bodies.at( "Mars" );
        gravitation::SphericalHarmonicsBlock cosine(
                [ gravity ]( ) -> Eigen::MatrixXd& { return gravity->getCosineCoefficientsReference( ); }, degree + 1, degree + 1 );
        gravitation::SphericalHarmonicsBlock sine(
                [ gravity ]( ) -> Eigen::MatrixXd& { return gravity->getSineCoefficientsReference( ); }, degree + 1, degree + 1 );
        auto model =
                std::make_shared< ExperimentalAcceleration >( [ satellite ]( Eigen::Vector3d& out ) { out = satellite->getPosition( ); },
                                                              [ gravity ]( ) { return gravity->getGravitationalParameter( ); },
                                                              gravity->getReferenceRadius( ),
                                                              cosine,
                                                              sine,
                                                              [ mars ]( Eigen::Vector3d& out ) { out = mars->getPosition( ); },
                                                              [ mars ]( ) { return mars->getCurrentRotationToGlobalFrame( ); } );
        model->setReuseTermBuffer( reuseTermBuffer );
        acceleration[ "Satellite" ][ "Mars" ].push_back( model );
    }
    auto settings = std::make_shared< TranslationalStatePropagatorSettings< double, Time > >(
            std::vector< std::string >{ "Mars" },
            acceleration,
            std::vector< std::string >{ "Satellite" },
            state,
            Time( 778118400. ),
            numerical_integrators::rungeKuttaFixedStepSettings< Time >( Time( 30. ), numerical_integrators::rungeKutta87DormandPrince ),
            std::make_shared< PropagationTimeTerminationSettings >( 778118400. + duration, true ) );
    Run result;
    result.state = state;
    if( variational )
    {
        auto ps = getInitialStateParameterSettings< double, Time >( settings, bodies );
        ps.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                "Mars", estimatable_parameters::gravitational_parameter ) );
        auto parameters = createParametersToEstimate< double, Time >( ps, bodies );
        result.solver = std::make_shared< Solver >( bodies, settings, parameters, true, false );
        result.dynamics = result.solver->getDynamicsSimulator( );
    }
    else
        result.dynamics = std::make_shared< Dynamics >( bodies, settings, false );
    return result;
}
void integrate( Run& run )
{
    if( run.solver )
        run.solver->integrateVariationalAndDynamicalEquations( run.state, true );
    else
        run.dynamics->integrateEquationsOfMotion( run.state );
}
void require( bool x, const char* msg )
{
    if( !x ) throw std::runtime_error( msg );
}
double maxAbs = 0, maxRatio = 0;
template< class M >
void compare( const M& a, const M& b )
{
    require( a.size( ) == b.size( ), "history count mismatch" );
    auto j = b.begin( );
    for( const auto& x : a )
    {
        require( x.first == j->first, "epoch mismatch" );
        const auto& u = x.second;
        const auto& v = j->second;
        require( u.rows( ) == v.rows( ) && u.cols( ) == v.cols( ) && u.allFinite( ) && v.allFinite( ), "invalid matrix" );
        for( int i = 0; i < u.size( ); ++i )
        {
            double d = std::abs( u( i ) - v( i ) );
            maxAbs = std::max( maxAbs, d );
            maxRatio = std::max( maxRatio, d / ( 1e-14 + 1e-11 * std::abs( v( i ) ) ) );
        }
        ++j;
    }
}

void cacheChecks( )
{
    using namespace basic_mathematics;
    for( bool normalized : { false, true } )
        for( int degree : { 1, 20, 100 } )
            for( int order : { 0, 5, 100 } )
            {
                if( order > degree ) continue;
                for( double latitude : { -.99999, -.3, 0., .5, .99999 } )
                {
                    LegendreCache serial( degree, order, normalized ), parallel( degree, order, normalized );
                    serial.setComputeFirstDerivatives( true );
                    serial.setComputeSecondDerivatives( true );
                    parallel.setComputeFirstDerivatives( true );
                    parallel.setComputeSecondDerivatives( true );
                    experimentParallelLegendre = false;
                    serial.update( latitude );
                    experimentParallelLegendre = true;
                    omp_set_num_threads( 4 );
                    parallel.update( latitude );
                    for( int n = 0; n <= degree; ++n )
                        for( int m = 0; m <= std::min( n, order ); ++m )
                        {
                            require( serial.getLegendrePolynomial( n, m ) == parallel.getLegendrePolynomial( n, m ),
                                     "Legendre value mismatch" );
                            require( serial.getLegendrePolynomialDerivative( n, m ) == parallel.getLegendrePolynomialDerivative( n, m ),
                                     "Legendre derivative mismatch" );
                            require( serial.getLegendrePolynomialSecondDerivative( n, m ) ==
                                             parallel.getLegendrePolynomialSecondDerivative( n, m ),
                                     "Legendre second derivative mismatch" );
                        }
                }
            }
    experimentParallelLegendre = false;
    omp_set_num_threads( 1 );
}
void bufferChecks( )
{
    using namespace gravitation;
    omp_set_num_threads( 4 );
    experimentParallelAcceleration = true;
    std::vector< Eigen::Vector3d > scratch;
    experimentTermBufferResizes = 0;
    for( int degree : { 3, 8, 3 } )
    {
        Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( degree + 1, degree + 1 );
        Eigen::MatrixXd sine = cosine;
        cosine( 0, 0 ) = 1.0;
        cosine( 2, 0 ) = -0.001;
        sine( 2, 2 ) = 0.00001;
        basic_mathematics::SphericalHarmonicsCache cache( degree, degree );
        std::map< std::pair< int, int >, Eigen::Vector3d > unused;
        for( int sample = 0; sample < 2; ++sample )
        {
            Eigen::Vector3d position( 4.0e6, 1.0e5 * ( sample + 1 ), 2.0e5 );
            const auto expected = computeGeodesyNormalizedGravitationalAccelerationSum(
                    position, 4.2828e13, 3.396e6, cosine, sine, cache, unused, false, Eigen::Matrix3d::Identity( ), true, 1.0e-6, nullptr );
            const auto actual = computeGeodesyNormalizedGravitationalAccelerationSum( position,
                                                                                      4.2828e13,
                                                                                      3.396e6,
                                                                                      cosine,
                                                                                      sine,
                                                                                      cache,
                                                                                      unused,
                                                                                      false,
                                                                                      Eigen::Matrix3d::Identity( ),
                                                                                      true,
                                                                                      1.0e-6,
                                                                                      &scratch );
            require( actual.allFinite( ) && ( expected.array( ) == actual.array( ) ).all( ), "buffer resize/reuse mismatch" );
            for( auto& term : scratch ) term.setConstant( TUDAT_NAN );
            const auto poisoned = computeGeodesyNormalizedGravitationalAccelerationSum( position,
                                                                                        4.2828e13,
                                                                                        3.396e6,
                                                                                        cosine,
                                                                                        sine,
                                                                                        cache,
                                                                                        unused,
                                                                                        false,
                                                                                        Eigen::Matrix3d::Identity( ),
                                                                                        true,
                                                                                        1.0e-6,
                                                                                        &scratch );
            require( poisoned.allFinite( ) && ( expected.array( ) == poisoned.array( ) ).all( ), "stale term consumed" );
        }
    }
    require( experimentTermBufferResizes == 3, "unexpected buffer resizing" );
    auto untouched = scratch;
    experimentParallelAcceleration = false;
    Eigen::MatrixXd cosine = Eigen::MatrixXd::Identity( 4, 4 ), sine = Eigen::MatrixXd::Zero( 4, 4 );
    basic_mathematics::SphericalHarmonicsCache cache( 3, 3 );
    std::map< std::pair< int, int >, Eigen::Vector3d > unused;
    computeGeodesyNormalizedGravitationalAccelerationSum( Eigen::Vector3d( 4e6, 1e5, 2e5 ),
                                                          4.2828e13,
                                                          3.396e6,
                                                          cosine,
                                                          sine,
                                                          cache,
                                                          unused,
                                                          false,
                                                          Eigen::Matrix3d::Identity( ),
                                                          true,
                                                          1e-6,
                                                          &scratch );
    require( scratch.size( ) == untouched.size( ) && experimentTermBufferResizes == 3, "serial path changed scratch" );
    omp_set_num_threads( 1 );
}

int main( int argc, char** argv )
{
    try
    {
        int degree = argc > 1 ? std::stoi( argv[ 1 ] ) : 100;
        double duration = argc > 2 ? std::stod( argv[ 2 ] ) : 21600.;
        std::string mode = argc > 3 ? argv[ 3 ] : "variational";
        int repeats = argc > 4 ? std::stoi( argv[ 4 ] ) : 3;
        require( degree >= 2 && degree <= 120 && duration >= 300 && duration / 30 == std::floor( duration / 30 ) && repeats >= 1 &&
                         repeats <= 10,
                 "invalid workload" );
        require( mode == "dynamics" || mode == "variational", "invalid mode" );
        bool variational = mode == "variational";
        omp_set_dynamic( 0 );
        cacheChecks( );
        bufferChecks( );
        spice_interface::loadStandardSpiceKernels( );
        auto reference = makeRun( degree, duration, variational, true );
        integrate( reference );
        struct Choice {
            const char* name;
            int threads;
            bool acceleration;
            bool legendre;
            bool reuse;
        };
        std::vector< Choice > choices = { { "sequential", 1, false, false, false } };
        for( int threads : { 2, 4, 8 } )
        {
            choices.push_back( { "temporary", threads, true, true, false } );
            choices.push_back( { "reuse", threads, true, true, true } );
        }
        std::cout << std::setprecision( 12 );
        for( int rep = 0; rep < repeats; ++rep )
            for( int k = 0; k < static_cast< int >( choices.size( ) ); ++k )
            {
                auto choice = choices[ ( k + rep * 3 ) % choices.size( ) ];
                omp_set_num_threads( choice.threads );
                gravitation::experimentParallelAcceleration = choice.acceleration;
                basic_mathematics::experimentParallelLegendre = choice.legendre;
                auto start = Clock::now( );
                auto run = makeRun( degree, duration, variational, false, choice.reuse );
                auto setup = Clock::now( );
                gravitation::experimentAccelerationCalls = 0;
                gravitation::experimentTermBufferResizes = 0;
                gravitation::experimentTemporaryTermBuffers = 0;
                basic_mathematics::experimentLegendreUpdates = 0;
                integrate( run );
                auto end = Clock::now( );
                require( gravitation::experimentAccelerationCalls > 0, "experimental acceleration was not linked" );
                require( basic_mathematics::experimentLegendreUpdates > 0, "experimental Legendre update was not linked" );
                if( choice.acceleration )
                    require( gravitation::experimentAccelerationThreads == choice.threads, "wrong acceleration team size" );
                if( choice.legendre ) require( basic_mathematics::experimentLegendreThreads == choice.threads, "wrong Legendre team size" );
                if( choice.reuse )
                {
                    require( gravitation::experimentTermBufferResizes == 1, "expected one model-owned buffer resize" );
                    require( gravitation::experimentTemporaryTermBuffers == 0, "unexpected temporary buffer" );
                }
                else if( choice.acceleration )
                {
                    require( gravitation::experimentTemporaryTermBuffers == gravitation::experimentAccelerationCalls,
                             "expected temporary buffer for every acceleration call" );
                }
                maxAbs = 0;
                maxRatio = 0;
                auto actual = run.dynamics->getSingleArcPropagationResults( );
                auto expected = reference.dynamics->getSingleArcPropagationResults( );
                compare( actual->getEquationsOfMotionNumericalSolution( ), expected->getEquationsOfMotionNumericalSolution( ) );
                if( variational )
                {
                    auto a = run.solver->getSingleArcVariationalPropagationResults( );
                    auto b = reference.solver->getSingleArcVariationalPropagationResults( );
                    compare( a->getStateTransitionSolution( ), b->getStateTransitionSolution( ) );
                    compare( a->getSensitivitySolution( ), b->getSensitivitySolution( ) );
                }
                std::cout << "{\"degree\":" << degree << ",\"duration\":" << duration << ",\"mode\":\"" << mode << "\",\"variant\":\""
                          << choice.name << "\",\"repetition\":" << rep + 1 << ",\"threads\":" << choice.threads
                          << ",\"setup_seconds\":" << std::chrono::duration< double >( setup - start ).count( )
                          << ",\"integration_seconds\":" << std::chrono::duration< double >( end - setup ).count( )
                          << ",\"acceleration_calls\":" << gravitation::experimentAccelerationCalls
                          << ",\"legendre_updates\":" << basic_mathematics::experimentLegendreUpdates
                          << ",\"buffer_resizes\":" << gravitation::experimentTermBufferResizes
                          << ",\"temporary_buffers\":" << gravitation::experimentTemporaryTermBuffers
                          << ",\"epochs\":" << actual->getEquationsOfMotionNumericalSolution( ).size( )
                          << ",\"max_abs_difference\":" << maxAbs << ",\"max_tolerance_ratio\":" << maxRatio << ",\"validation\":\""
                          << ( maxRatio <= 1 ? "passed" : "failed" ) << "\"}" << std::endl;
                require( maxRatio <= 1, "numerical parity failed" );
            }
    }
    catch( const std::exception& e )
    {
        std::cerr << e.what( ) << std::endl;
        return 1;
    }
}
