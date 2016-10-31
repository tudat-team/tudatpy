#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/numericalObservationPartial.h"

namespace tudat
{

namespace observation_partials
{


void emptyFunction2( ){ }

boost::function< void( ) > getNumericalObservationPartialUpdateFunction(
        estimatable_parameters::EstimatebleParameterIdentifier parameter,
        std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > transmittingBody,
        std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > receivingBody,
        std::string transmittingStationName,
        std::string receivingStationName )
{

    using namespace tudat::estimatable_parameters;
    boost::function< void( ) > updateFunction = emptyFunction2;

    switch( parameter.first )
    {
    case gravitational_parameter:
        break;
    case constant_rotation_rate:
        break;
    case constant_drag_coefficient:
        break;
    case spherical_harmonics_cosine_coefficient_block:
        break;
    case spherical_harmonics_sine_coefficient_block:
        break;
    default:
        std::cerr<<"Parameter "<<parameter.first<<"of body: "<<parameter.second.first<<" unknown when generating update function."<<std::endl;
        break;
    }

    return boost::function< void( ) >( );
}

double receptionTimeFunction( double receptionTime,
                              double transmissionTime )
{
    return receptionTime;
}

double transmissionTimeFunction( double receptionTime,
                                 double transmissionTime )
{
    return transmissionTime;
}

void resetInitialDynamicalState(
        const boost::shared_ptr< propagators::PropagatorSettings< double > > & propagatorSettings,
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const double perturbedParameterValue, const int componentIndex )
{
    propagators::IntegratedStateType stateType;
    switch( parameter->getParameterName( ).first )
    {
    case estimatable_parameters::initial_body_state:
        stateType = propagators::transational_state;
        break;
    default:
        std::cerr<<"Error when resetting initial dynamical state, did not recognize type "<<parameter->getParameterName( ).first<<std::endl;

    }


    {
        Eigen::VectorXd initialDynamicalState = propagatorSettings->getInitialStates( );
        initialDynamicalState( componentIndex ) = perturbedParameterValue;
        propagatorSettings->resetInitialStates( initialDynamicalState );
    }
}


Eigen::Matrix< double, 1, 1 > calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< double( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction )
{

    double unperturbedParameterValue = parameter->getParameterValue( );

    parameter->setParameterValue( unperturbedParameterValue + parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, 1, 1 > upPerturbedValue;
    upPerturbedValue( 0, 0 ) = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue - parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, 1, 1 > downPerturbedValue;
    downPerturbedValue( 0, 0 ) = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return ( upPerturbedValue - downPerturbedValue ) / ( 2.0 * parameterPerturbation );
}

Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction )
{

    double unperturbedParameterValue = parameter->getParameterValue( );

    parameter->setParameterValue( unperturbedParameterValue + parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 >  upPerturbedValue = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue - parameterPerturbation );
    updateFunction( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 >  downPerturbedValue = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return ( upPerturbedValue - downPerturbedValue ) / ( 2.0 * parameterPerturbation );
}

Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateNumericalObservationParameterPartialWithSingleArcDynamicsUpdate(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings,
        const boost::shared_ptr< propagators::PropagatorSettings< double > > & propagatorSettings,
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction )
{

    double unperturbedParameterValue = parameter->getParameterValue( );

    parameter->setParameterValue( unperturbedParameterValue + parameterPerturbation );
    updateFunction( );
    propagators::SingleArcDynamicsSimulator< >( bodyMap, integratorSettings, propagatorSettings );
    Eigen::Matrix< double, Eigen::Dynamic, 1 >  upPerturbedValue = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue - parameterPerturbation );
    updateFunction( );
    propagators::SingleArcDynamicsSimulator< >( bodyMap, integratorSettings, propagatorSettings );
    Eigen::Matrix< double, Eigen::Dynamic, 1 >  downPerturbedValue = observationFunction( evaluationTime );

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return ( upPerturbedValue - downPerturbedValue ) / ( 2.0 * parameterPerturbation );
}

Eigen::MatrixXd calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const Eigen::VectorXd parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction )
{
    Eigen::MatrixXd parameterPartial = Eigen::MatrixXd::Zero( observationFunction( evaluationTime ).rows( ),
                                                              parameter->getParameterSize( ) );

    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );
    Eigen::VectorXd perturbedParameterValue;

    for( int i = 0; i < parameter->getParameterSize( ); i++ )
    {
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) += parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateFunction( );
        Eigen::VectorXd upPerturbedValue = observationFunction( evaluationTime );

        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) -= parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateFunction( );
        Eigen::VectorXd downPerturbedValue = observationFunction( evaluationTime );

        parameterPartial.block( 0, i, downPerturbedValue.rows( ), 1 ) = ( upPerturbedValue - downPerturbedValue ) /
                ( 2.0 * parameterPerturbation( i ) );
    }

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return parameterPartial;
}

Eigen::MatrixXd calculateNumericalObservationParameterPartialWithSingleArcDynamicsUpdate(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings,
        const boost::shared_ptr< propagators::PropagatorSettings< double > > & propagatorSettings,
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const Eigen::VectorXd parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction )
{
    Eigen::MatrixXd parameterPartial = Eigen::MatrixXd::Zero( observationFunction( evaluationTime ).rows( ),
                                                              parameter->getParameterSize( ) );

    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );
    Eigen::VectorXd perturbedParameterValue;

    for( int i = 0; i < parameter->getParameterSize( ); i++ )
    {
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) += parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        if( parameter->getParameterName( ).first < 0 )
        {
            resetInitialDynamicalState( propagatorSettings, parameter, perturbedParameterValue( i ), i );
        }

        propagators::SingleArcDynamicsSimulator< >( bodyMap, integratorSettings, propagatorSettings );
        updateFunction( );
        Eigen::VectorXd upPerturbedValue = observationFunction( evaluationTime );

        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) -= parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        if( parameter->getParameterName( ).first < 0 )
        {
            resetInitialDynamicalState( propagatorSettings, parameter, perturbedParameterValue( i ), i );
        }

        propagators::SingleArcDynamicsSimulator< >( bodyMap, integratorSettings, propagatorSettings );
        updateFunction( );
        Eigen::VectorXd downPerturbedValue = observationFunction( evaluationTime );

        if( parameter->getParameterName( ).first < 0 )
        {
            resetInitialDynamicalState( propagatorSettings, parameter, unperturbedParameterValue( i ), i );
        }

        parameterPartial.block( 0, i, downPerturbedValue.rows( ), 1 ) = ( upPerturbedValue - downPerturbedValue ) /
                ( 2.0 * parameterPerturbation( i ) );
    }

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return parameterPartial;
}

}

}
