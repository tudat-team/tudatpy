#ifndef TUDATPY_INDEPENDENT_ARC_MODELS_H
#define TUDATPY_INDEPENDENT_ARC_MODELS_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>
#include "scalarTypes.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{
namespace independent_arcs
{

namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tep = tudat::estimatable_parameters;
namespace te = tudat::ephemerides;
using Scalar = STATE_SCALAR_TYPE;
using Epoch = TIME_TYPE;
using Vector = Eigen::Matrix< Scalar, Eigen::Dynamic, 1 >;
using Settings = tp::SingleArcPropagatorSettings< Scalar, Epoch >;
using SingleSolver = tp::SingleArcVariationalEquationsSolver< Scalar, Epoch >;
using MultiSolver = tp::MultiArcVariationalEquationsSolver< Scalar, Epoch >;
using MultiSettings = tp::MultiArcPropagatorSettings< Scalar, Epoch >;
using Parameters = tep::EstimatableParameterSet< Scalar >;

// One mutex for this binding's supported SPICE delegates, not one per worker.
inline std::mutex& spiceMutex( )
{
    static std::mutex mutex;
    return mutex;
}

class LockedRotation final : public te::RotationalEphemeris
{
    std::shared_ptr< te::RotationalEphemeris > source_;

public:
    explicit LockedRotation( std::shared_ptr< te::RotationalEphemeris > source ):
        RotationalEphemeris( source->getBaseFrameOrientation( ), source->getTargetFrameOrientation( ) ), source_( std::move( source ) )
    {}
    Eigen::Quaterniond getRotationToBaseFrame( double time ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex( ) );
        return source_->getRotationToBaseFrame( time );
    }
    Eigen::Quaterniond getRotationToTargetFrame( double time ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex( ) );
        return source_->getRotationToTargetFrame( time );
    }
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( double time ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex( ) );
        return source_->getDerivativeOfRotationToBaseFrame( time );
    }
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( double time ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex( ) );
        return source_->getDerivativeOfRotationToTargetFrame( time );
    }
    void getFullRotationalQuantitiesToTargetFrame( Eigen::Quaterniond& rotation,
                                                   Eigen::Matrix3d& derivative,
                                                   Eigen::Vector3d& velocity,
                                                   double time ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex( ) );
        source_->getFullRotationalQuantitiesToTargetFrame( rotation, derivative, velocity, time );
    }
};

class LockedEphemeris final : public te::Ephemeris
{
    std::shared_ptr< te::Ephemeris > source_;

public:
    explicit LockedEphemeris( std::shared_ptr< te::Ephemeris > source ):
        Ephemeris( source->getReferenceFrameOrigin( ), source->getReferenceFrameOrientation( ) ), source_( std::move( source ) )
    {}
    Eigen::Vector6d getCartesianState( double time ) override
    {
        std::lock_guard< std::mutex > lock( spiceMutex( ) );
        return source_->getCartesianState( time );
    }
};

struct MarsGravityArcModel {
    const int degree;
    const int order;
    double referenceRadius;
    double gravitationalParameter;

    MarsGravityArcModel( int degreeValue = 100, int orderValue = -1 ):
        degree( degreeValue ), order( orderValue == -1 ? degreeValue : orderValue )
    {
        if( degree < 2 || degree > 120 || order < 0 || order > degree )
        {
            throw std::invalid_argument( "require 2 <= degree <= 120 and 0 <= order <= degree (or order=-1)" );
        }
        const auto bodies = makeBodies( false );
        auto gravity = std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >(
                bodies.at( "Mars" )->getGravityFieldModel( ) );
        referenceRadius = gravity->getReferenceRadius( );
        gravitationalParameter = gravity->getGravitationalParameter( );
    }

    tss::SystemOfBodies makeBodies( bool locked ) const
    {
        auto settings = tss::getDefaultBodySettings( { "Mars" }, "Mars", "J2000" );
        settings.at( "Mars" )->gravityFieldSettings =
                std::make_shared< tss::FromFileSphericalHarmonicsGravityFieldSettings >( tss::jgmro120d, degree );
        auto bodies = tss::createSystemOfBodies< Scalar, Epoch >( settings );
        if( locked )
        {
            auto mars = bodies.at( "Mars" );
            mars->setRotationalEphemeris( std::make_shared< LockedRotation >( mars->getRotationalEphemeris( ) ) );
            mars->setEphemeris( std::make_shared< LockedEphemeris >( mars->getEphemeris( ) ) );
        }
        bodies.createEmptyBody< Scalar, Epoch >( "Satellite" );
        return bodies;
    }
};

struct ArcConfiguration {
    MarsGravityArcModel model;
    std::vector< double > starts;
    double duration;
    double step;

    std::shared_ptr< Settings > makeSettings( tss::SystemOfBodies& bodies, std::size_t arc, const Vector& state ) const
    {
        tss::SelectedAccelerationMap selected;
        selected[ "Satellite" ][ "Mars" ].push_back(
                std::make_shared< tss::SphericalHarmonicAccelerationSettings >( model.degree, model.order ) );
        const auto accelerations = tss::createAccelerationModelsMap( bodies, selected, { "Satellite" }, { "Mars" } );
        return std::make_shared< tp::TranslationalStatePropagatorSettings< Scalar, Epoch > >(
                std::vector< std::string >{ "Mars" },
                accelerations,
                std::vector< std::string >{ "Satellite" },
                state,
                Epoch( starts.at( arc ) ),
                tudat::numerical_integrators::rungeKuttaFixedStepSettings< Epoch >(
                        Epoch( step ), tudat::numerical_integrators::rungeKutta87DormandPrince ),
                std::make_shared< tp::PropagationTimeTerminationSettings >( starts.at( arc ) + duration, true ) );
    }
};

inline std::vector< std::shared_ptr< tep::EstimatableParameterSettings > > singleParameters( const std::shared_ptr< Settings >& settings,
                                                                                             const tss::SystemOfBodies& bodies )
{
    auto result = tss::getInitialStateParameterSettings< Scalar, Epoch >( settings, bodies );
    result.push_back( std::make_shared< tep::EstimatableParameterSettings >( "Mars", tep::gravitational_parameter ) );
    return result;
}

// Deliberately not derived from/exposing the mutable Python simulator wrapper.
// No member or worker capture contains py::object, Python callbacks or borrowed arrays.
class IndependentArcVariationalSimulator
{
    ArcConfiguration configuration_;
    bool transfer_;
    std::size_t workers_ = 1;
    enum class ResultMode { none, dynamics, variational };
    ResultMode mode_ = ResultMode::none;
    tss::SystemOfBodies bodies_;
    std::shared_ptr< MultiSettings > settings_;
    std::shared_ptr< Parameters > parameters_;
    std::shared_ptr< MultiSolver > solver_;

    void validateStates( const Vector& states, bool allowTransferNaNs ) const
    {
        if( states.size( ) != static_cast< Eigen::Index >( 6 * configuration_.starts.size( ) ) )
        {
            throw std::invalid_argument( "initial_states must contain six components per arc" );
        }
        for( std::size_t arc = 0; arc < configuration_.starts.size( ); ++arc )
        {
            const auto state = states.segment( 6 * arc, 6 );
            if( allowTransferNaNs && arc > 0 && state.array( ).isNaN( ).all( ) )
            {
                continue;
            }
            if( !state.allFinite( ) || !std::isfinite( state.head( 3 ).stableNorm( ) ) ||
                state.head( 3 ).stableNorm( ) <= configuration_.model.referenceRadius )
            {
                throw std::invalid_argument( "initial states must be finite and above the Mars reference surface" );
            }
        }
    }

    void setInitialStates( const Vector& states )
    {
        validateStates( states, transfer_ );
        std::vector< Vector > split;
        for( std::size_t arc = 0; arc < configuration_.starts.size( ); ++arc )
        {
            split.push_back( states.segment( 6 * arc, 6 ) );
        }
        settings_->resetInitialStatesList( split );
    }

    void requireResults( bool variational ) const
    {
        if( mode_ == ResultMode::none || ( variational && mode_ != ResultMode::variational ) )
        {
            throw std::runtime_error( variational ? "no current successful variational results; integrate_full_equations first"
                                                  : "no current successful propagation results" );
        }
    }

    void updateParametersWithoutIntegration( const Vector& values )
    {
        // The native resetParameterEstimate(values, false) integrates dynamics
        // only: false does NOT mean 'defer all integration'. Keep the Python
        // facade's explicit reintegrate=False contract without calling it.
        parameters_->resetParameterValues( values );
        tss::setInitialStateVectorFromParameterSet< Scalar, Epoch >( parameters_, settings_ );
        const auto localParameters = solver_->getArcWiseParametersToEstimate( );
        for( std::size_t arc = 0; arc < configuration_.starts.size( ); ++arc )
        {
            Vector localValues = localParameters.at( arc )->template getFullParameterValues< Scalar >( );
            localValues.head( 6 ) = settings_->getSingleArcSettings( ).at( arc )->getInitialStates( );
            localParameters.at( arc )->resetParameterValues( localValues );
        }
    }

public:
    IndependentArcVariationalSimulator( const MarsGravityArcModel& model,
                                        std::vector< double > starts,
                                        const Eigen::MatrixXd& initialStates,
                                        double duration,
                                        double step = 30.0,
                                        std::int64_t workers = 1,
                                        bool transfer = false ):
        configuration_{ model, std::move( starts ), duration, step }, transfer_( transfer )
    {
        const auto count = configuration_.starts.size( );
        if( count == 0 || count > static_cast< std::size_t >( std::numeric_limits< int >::max( ) / 6 ) ||
            initialStates.rows( ) != static_cast< Eigen::Index >( count ) || initialStates.cols( ) != 6 )
        {
            throw std::invalid_argument( "initial_states must have shape (number_of_arcs, 6), with at least one arc" );
        }
        if( !std::isfinite( duration ) || !std::isfinite( step ) || duration <= 0 || step <= 0 || duration / step < 8 )
        {
            throw std::invalid_argument( "duration and step must be finite and positive, with at least eight steps per arc" );
        }
        // Time stores full hours in an int; reject values before conversion can overflow.
        const double timeLimit = ( static_cast< double >( std::numeric_limits< int >::max( ) ) - 1.0 ) * tudat::TIME_NORMALIZATION_TERM;
        if( duration >= timeLimit || step >= timeLimit )
        {
            throw std::invalid_argument( "duration/step exceeds the native Time representation" );
        }
        if( workers < 1 )
        {
            throw std::invalid_argument( "number_of_workers must be positive" );
        }
        for( std::size_t arc = 0; arc < count; ++arc )
        {
            const double epoch = configuration_.starts.at( arc );
            if( !std::isfinite( epoch ) || !std::isfinite( epoch + duration ) || std::abs( epoch ) >= timeLimit ||
                std::abs( epoch + duration ) >= timeLimit || epoch + duration <= epoch || epoch + step <= epoch ||
                epoch + duration - step >= epoch + duration || ( arc > 0 && epoch <= configuration_.starts.at( arc - 1 ) ) )
            {
                throw std::invalid_argument(
                        "arc_start_times must be finite and increasing; end epochs/steps must be representable and advance time" );
            }
        }
        Vector concatenated( 6 * count );
        for( std::size_t arc = 0; arc < count; ++arc )
        {
            concatenated.segment( 6 * arc, 6 ) = initialStates.row( arc ).transpose( ).cast< Scalar >( );
        }
        validateStates( concatenated, false );
        bodies_ = configuration_.model.makeBodies( false );
        bodies_.at( "Satellite" )
                ->setEphemeris( std::make_shared< te::MultiArcEphemeris >(
                        std::map< double, std::shared_ptr< te::Ephemeris > >{
                                { 0.0, std::make_shared< te::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "Mars", "J2000" ) } },
                        "Mars",
                        "J2000" ) );
        bodies_.processBodyFrameDefinitions< Scalar, Epoch >( );
        std::vector< std::shared_ptr< Settings > > arcs;
        for( std::size_t arc = 0; arc < count; ++arc )
        {
            arcs.push_back( configuration_.makeSettings( bodies_, arc, concatenated.segment( 6 * arc, 6 ) ) );
        }
        settings_ = std::make_shared< MultiSettings >( arcs, transfer );
        settings_->getOutputSettings( )->setIntegratedResult( true );
        auto parameters = tss::getInitialMultiArcParameterSettings< Scalar, Epoch >( settings_, bodies_, configuration_.starts );
        parameters.push_back( std::make_shared< tep::EstimatableParameterSettings >( "Mars", tep::gravitational_parameter ) );
        parameters_ = tss::createParametersToEstimate< Scalar, Epoch >( parameters, bodies_ );
        solver_ = std::make_shared< MultiSolver >( bodies_, settings_, parameters_, false );
        setIndependentArcExecution( workers );
    }

    void setIndependentArcExecution( std::int64_t workers )
    {
        if( workers < 1 )
        {
            throw std::invalid_argument( "number_of_workers must be positive" );
        }
        if( workers == 1 )
        {
            solver_->clearIndependentArcExecution( );
        }
        else
        {
            // Capture the configuration by value, never this/the Python owner.
            solver_->setIndependentArcExecution(
                    workers,
                    [ configuration = configuration_ ]( std::size_t arc, const std::shared_ptr< Settings >&, const Vector& values ) {
                        auto bodies = configuration.model.makeBodies( true );
                        auto settings = configuration.makeSettings( bodies, arc, values.head( 6 ) );
                        auto parameters = tss::createParametersToEstimate< Scalar, Epoch >( singleParameters( settings, bodies ), bodies );
                        parameters->resetParameterValues( values );
                        return std::make_shared< SingleSolver >( bodies, settings, parameters, true, false );
                    } );
        }
        workers_ = workers;
    }

    void clearIndependentArcExecution( )
    {
        setIndependentArcExecution( 1 );
    }
    std::size_t getNumberOfWorkers( ) const
    {
        return workers_;
    }
    bool lastIntegrationUsedIndependentArcs( ) const
    {
        return mode_ != ResultMode::none && solver_->getDynamicsSimulator( )->getLastIntegrationUsedIndependentArcs( );
    }
    Vector getInitialStates( ) const
    {
        return settings_->getInitialStates( );
    }
    Vector getParameterVector( ) const
    {
        return parameters_->getFullParameterValues< Scalar >( );
    }

    void integrateFullEquations( const std::optional< Vector >& states = std::nullopt )
    {
        if( states )
        {
            setInitialStates( *states );
        }
        mode_ = ResultMode::none;
        solver_->integrateVariationalAndDynamicalEquations( settings_->getInitialStateList( ), true );
        if( !solver_->getDynamicsSimulator( )->integrationCompletedSuccessfully( ) )
        {
            throw std::runtime_error( "native multi-arc integration did not complete successfully" );
        }
        mode_ = ResultMode::variational;
    }

    void integrateEquationsOfMotionOnly( const std::optional< Vector >& states = std::nullopt )
    {
        if( states )
        {
            setInitialStates( *states );
        }
        mode_ = ResultMode::none;
        solver_->integrateDynamicalEquationsOfMotionOnly( settings_->getInitialStates( ) );
        if( !solver_->getDynamicsSimulator( )->integrationCompletedSuccessfully( ) )
        {
            throw std::runtime_error( "native multi-arc integration did not complete successfully" );
        }
        Vector values = getParameterVector( );
        values.head( settings_->getInitialStates( ).size( ) ) = settings_->getInitialStates( );
        updateParametersWithoutIntegration( values );
        mode_ = ResultMode::dynamics;
    }

    void resetParameterEstimate( const Vector& values, bool reintegrate = true )
    {
        if( values.size( ) != static_cast< Eigen::Index >( 6 * configuration_.starts.size( ) + 1 ) || !values.allFinite( ) ||
            values( values.size( ) - 1 ) <= 0 )
        {
            throw std::invalid_argument( "parameter_vector must contain six finite state components per arc and one positive Mars mu" );
        }
        validateStates( values.head( values.size( ) - 1 ), false );
        mode_ = ResultMode::none;
        updateParametersWithoutIntegration( values );
        if( reintegrate )
        {
            integrateFullEquations( );
        }
    }

    std::shared_ptr< tp::MultiArcSimulationResults< tp::SingleArcSimulationResults, Scalar, Epoch > > getPropagationResults( ) const
    {
        requireResults( false );
        return solver_->getDynamicsSimulator( )->getMultiArcPropagationResults( );
    }
    std::shared_ptr< tp::MultiArcSimulationResults< tp::SingleArcVariationalSimulationResults, Scalar, Epoch > > getVariationalResults( )
            const
    {
        requireResults( true );
        return solver_->getMultiArcVariationalPropagationResults( );
    }
    std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionInterface( ) const
    {
        requireResults( true );
        return solver_->getStateTransitionMatrixInterface( );
    }
};

}  // namespace independent_arcs
}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
#endif
