#ifndef NUMERICALOBSERVATIONPARTIAL_H
#define NUMERICALOBSERVATIONPARTIAL_H

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace observation_partials
{

void emptyFunction2( );

boost::function< void( ) > getNumericalObservationPartialUpdateFunction(
        estimatable_parameters::EstimatebleParameterIdentifier parameter,
        std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > transmittingBody,
        std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > receivingBody,
        std::string transmittingStationName,
        std::string receivingStationName );

struct TimeEvaluator
{
public:
    TimeEvaluator( const int index ):index_( index ){ }

    double getTime( const std::vector< double > times ){ return times[ index_ ]; }

    int index_;
};

double receptionTimeFunction( double receptionTime,
                              double transmissionTime );

double transmissionTimeFunction( double receptionTime,
                                 double transmissionTime );

void resetInitialDynamicalState(
        const boost::shared_ptr< propagators::PropagatorSettings< double > > & propagatorSettings,
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const double parameterPerturbation, const int componentIndex );

Eigen::Matrix< double, 1, 1 > calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< double( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyFunction2 );

Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyFunction2 );


Eigen::MatrixXd calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const Eigen::VectorXd parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyFunction2 );


/*
Eigen::Matrix< double, 1, Eigen::Dynamic >
calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const Eigen::VectorXd parameterPerturbation,
        boost::function< double( double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyFunction2 )
{
    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );
    int parameterSize = unperturbedParameterValue.size( );

    Eigen::Matrix< double, 1, Eigen::Dynamic > numericalObservationPartials =
            Eigen::Matrix< double, 1, Eigen::Dynamic >::Zero( 1, parameterSize );

    Eigen::Matrix< double, 1, 1 > newObservation;
    Eigen::VectorXd perturbedParameterValue;
    for( int i = 0; i < parameterSize; i++ )
    {
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue[ i ] += parameterPerturbation[ i ];

        parameter->setParameterValue( perturbedParameterValue );
        updateFunction( );
        newObservation( 0, 0 ) = observationFunction( evaluationTime );
        numericalObservationPartials.block( 0, i, 1, 1 ) =
                numericalObservationPartials.block( 0, i, 1, 1 ) + newObservation;

        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue[ i ] -= parameterPerturbation[ i ];

        parameter->setParameterValue( perturbedParameterValue );
        updateFunction( );
        newObservation( 0, 0 ) = observationFunction( evaluationTime );
        numericalObservationPartials.block( 0, i, 1, 1 ) =
                numericalObservationPartials.block( 0, i, 1, 1 ) - newObservation;

        numericalObservationPartials.block( 0, i, 1, 1 ) /= ( 2.0 * parameterPerturbation[ i ] );
    }

    parameter->setParameterValue( unperturbedParameterValue );
    updateFunction( );

    return numericalObservationPartials;
}
*/

template< int ObservationSize = 1, typename ParameterType = double >
class NumericalObservationParameterPartial: public ObservationPartial< ObservationSize >
{
public:

    NumericalObservationParameterPartial(
            boost::shared_ptr< observation_models::ObservationModel< ObservationSize > > observationModel,
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameter,
            ParameterType relativeParameterPerturbation,
            ParameterType absoluteParameterPerturbation,
            boost::function< double( const std::vector< double > ) > evaluationTimeFunction,
            boost::function< void( ) > updateFunction = emptyFunction2 ):
        parameter_( parameter ),
        relativeParameterPerturbation_( relativeParameterPerturbation ),
        absoluteParameterPerturbation_( absoluteParameterPerturbation ),
        evaluationTimeFunction_( evaluationTimeFunction ),
        updateFunction_( updateFunction )
    {
        observationFunction_ = boost::bind( &observation_models::ObservationModel< ObservationSize >::computeObservations,
                                            observationModel, _1 );
    }

    ~NumericalObservationParameterPartial( ) { }

    std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > calculatePartial(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime )
    {
        ParameterType unperturbedParameterValue = parameter_->getParameterValue( );
        ParameterType parameterPerturbation = determineParameterPerturbation(
                    unperturbedParameterValue );

        return make_pair( calculateNumericalObservationParameterPartial(
                              parameter_, parameterPerturbation,
                              observationFunction_,
                              times ), evaluationTimeFunction_( times ) );
    }

private:
    boost::function< double( double ) > observationFunction_;

    boost::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType> > parameter_;

    ParameterType determineParameterPerturbation( const ParameterType parameterValue )
    {
        double perturbation = std::abs( parameterValue * relativeParameterPerturbation_ );
        assert( perturbation >= 0.0 );

        if( perturbation < absoluteParameterPerturbation_ )
        {
            perturbation = absoluteParameterPerturbation_;
        }

        return perturbation;
    }

    ParameterType relativeParameterPerturbation_;

    ParameterType absoluteParameterPerturbation_;

    boost::function< double( const std::vector< double > ) > evaluationTimeFunction_;

    boost::function< void( ) > updateFunction_;
};

}

}

#endif // NUMERICALOBSERVATIONPARTIAL_H
