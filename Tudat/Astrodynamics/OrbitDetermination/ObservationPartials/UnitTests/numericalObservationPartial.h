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

void emptyVoidFunction( );

//! Function to compute numerical partial derivative of double observable w.r.t. double parameter.
Eigen::Matrix< double, 1, 1 > calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< double( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyVoidFunction );

//! Function to compute numerical partial derivative of vector observable w.r.t. double parameter.
Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const double parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyVoidFunction );

//! Function to compute numerical partial derivative of vector observable w.r.t. vector parameter.
Eigen::MatrixXd calculateNumericalObservationParameterPartial(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const Eigen::VectorXd parameterPerturbation,
        boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double evaluationTime,
        boost::function< void( ) > updateFunction = &emptyVoidFunction );

}

}

#endif // NUMERICALOBSERVATIONPARTIAL_H
