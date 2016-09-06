#ifndef CREATEESTIMATABLEPARAMETERS_H
#define CREATEESTIMATABLEPARAMETERS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

template< typename InitialStateParameterType = double >
boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > >
createInitialDynamicalStateParameterToEstimate(
        const NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& parameterSettings )
{
    using namespace tudat::estimatable_parameters;

    boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > initialStateParameterToEstimate;

    if( !isParameterDynamicalPropertyInitialState( parameterSettings->parameterType_.first ) )
    {
        std::string errorMessage = "Error when requesting to make initial state parameter " +
                boost::lexical_cast< std::string >( parameterSettings->parameterType_.first ) + " of " +
                boost::lexical_cast< std::string >( parameterSettings->parameterType_.second.first ) +
                ", parameter is not an initial state parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        switch( parameterSettings->parameterType_.first )
        {
        case initial_body_state:
            if( boost::dynamic_pointer_cast< InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                        parameterSettings ) == NULL )
            {
                throw std::runtime_error( "Error when making body initial state parameter, settings type is incompatible" );
            }
            else
            {
                boost::shared_ptr< InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > > initialStateSettings =
                        boost::dynamic_pointer_cast< InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                if( ! ( initialStateSettings->initialTime_ == initialStateSettings->initialTime_  ) )
                {
                    initialStateParameterToEstimate = boost::make_shared< InitialTranslationalStateParameter< InitialStateParameterType > >(
                                initialStateSettings->parameterType_.second.first, initialStateSettings->initialStateValue_,
                                initialStateSettings->centralBody_,
                                initialStateSettings->frameOrientation_ );
                }
                else
                {
                    initialStateParameterToEstimate = boost::make_shared< InitialTranslationalStateParameter< InitialStateParameterType > >(
                                initialStateSettings->parameterType_.second.first, propagators::getInitialStateOfBody
                                < double, InitialStateParameterType >(
                                    initialStateSettings->parameterType_.second.first, initialStateSettings->centralBody_, bodyMap,
                                    initialStateSettings->initialTime_ ), initialStateSettings->centralBody_,
                                initialStateSettings->frameOrientation_ );
                }
            }
            break;
        default:
            std::string errorMessage = "Error, could not create parameter for initial state of type " +
                    boost::lexical_cast< std::string >( parameterSettings->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
    }

    return initialStateParameterToEstimate;
}


boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > createDoubleParameterToEstimate(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& doubleParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap = basic_astrodynamics::AccelerationMap( ) );

boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > createVectorParameterToEstimate(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& vectorParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap = basic_astrodynamics::AccelerationMap( ) );

template< typename InitialStateParameterType = double >
boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > createParametersToEstimate(
        const std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& parameterNames,
        const NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap = basic_astrodynamics::AccelerationMap( ) )

{
    using namespace tudat::estimatable_parameters;

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParametersToEstimate;
    std::vector< boost::shared_ptr< EstimatableParameter< double > > > doubleParametersToEstimate;
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParametersToEstimate;

    for( unsigned int i = 0; i < parameterNames.size( ); i++ )
    {
        if( isParameterDynamicalPropertyInitialState( parameterNames.at( i )->parameterType_.first ) )
        {
            initialDynamicalParametersToEstimate.push_back( createInitialDynamicalStateParameterToEstimate< InitialStateParameterType >(
                                                                bodyMap, parameterNames.at( i ) ) );
        }
        else if( isDoubleParameter( parameterNames[ i ]->parameterType_.first ) == true )
        {
            doubleParametersToEstimate.push_back( createDoubleParameterToEstimate(
                                                      parameterNames[ i ], bodyMap, accelerationModelMap ) );
        }
        else if( isDoubleParameter( parameterNames[ i ]->parameterType_.first ) == false )
        {
            vectorParametersToEstimate.push_back( createVectorParameterToEstimate(
                                                      parameterNames[ i ], bodyMap, accelerationModelMap ) );
        }
        else
        {
            std::string errorMessage = "Error, parameter type of  " +
                    boost::lexical_cast< std::string >( parameterNames[ i ]->parameterType_.second.first ) + "of " +
                    boost::lexical_cast< std::string >( parameterNames[ i ]->parameterType_.first ) +
                    "not recognized when making estimatable parameter set.";

            throw std::runtime_error( errorMessage );
        }
    }

    return boost::make_shared< EstimatableParameterSet< InitialStateParameterType > >(
                doubleParametersToEstimate, vectorParametersToEstimate, initialDynamicalParametersToEstimate );
}


}

}

#endif // CREATEESTIMATABLEPARAMETERS_H
