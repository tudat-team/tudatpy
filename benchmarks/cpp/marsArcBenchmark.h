#ifndef TUDAT_MARS_ARC_BENCHMARK_H
#define TUDAT_MARS_ARC_BENCHMARK_H
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
inline double seconds( Clock::time_point start )
{
    return std::chrono::duration< double >( Clock::now( ) - start ).count( );
}
inline std::mutex spiceMutex;

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

inline SystemOfBodies makeBodies( int degree, bool locked )
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

inline std::shared_ptr< Settings > makeSettings( SystemOfBodies& bodies, int arc, int degree, double days, double step )
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

inline std::vector< Context > makeContexts( int arcs, int degree, double days, double step )
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
inline Eigen::MatrixXd combine( const Eigen::MatrixXd& stm, const Eigen::MatrixXd& sensitivity, int arc, int arcs )
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

inline void checkMatrix( const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected, const std::string& label )
{
    if( actual.rows( ) != expected.rows( ) || actual.cols( ) != expected.cols( ) || !actual.allFinite( ) || !expected.allFinite( ) )
    {
        throw std::runtime_error( label + ": invalid dimensions or values" );
    }
    if( !( ( actual - expected ).array( ).abs( ) <= ( 1.0e-14 + 1.0e-11 * expected.array( ).abs( ) ) ).all( ) )
    {
        std::ostringstream message;
        message << std::setprecision( 17 ) << label << ": numerical mismatch; max absolute difference "
                << ( actual - expected ).cwiseAbs( ).maxCoeff( ) << "; actual " << actual.transpose( ) << "; expected "
                << expected.transpose( );
        throw std::runtime_error( message.str( ) );
    }
}

inline void testExecutor( )
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

#endif
