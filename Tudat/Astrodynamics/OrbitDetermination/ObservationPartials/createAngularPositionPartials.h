#ifndef CREATEANGULARPOSITIONPARTIALS_H
#define CREATEANGULARPOSITIONPARTIALS_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/createPositionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/angularPositionPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"


namespace tudat
{

namespace observation_partials
{

boost::shared_ptr< AngularPositionPartial > createAngularPositionPartialWrtInitialPosition(
        const observation_models::LinkEnds angularPositionLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const boost::shared_ptr< AngularPositionScaling > angularPositionScaler );


template< typename ParameterType >
boost::shared_ptr< AngularPositionPartial > createAngularPositionPartialWrtParameter(
        const observation_models::LinkEnds angularPositionLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const boost::shared_ptr< AngularPositionScaling > angularPositionScaler )
{
    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartials =
            createPositionPartialsWrtParameter( angularPositionLinkEnds, bodyMap, parameterToEstimate );
    boost::shared_ptr< AngularPositionPartial > angularPositionPartial;

    if( positionPartials.size( ) > 0 )
    {
        angularPositionPartial = boost::make_shared< AngularPositionPartial >( angularPositionScaler, positionPartials,
                                                                               parameterToEstimate->getParameterName( ) );
    }

    return angularPositionPartial;
}


template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 2 > > >, boost::shared_ptr< PositionPartialScaling > >
createAngularPositionPartials( const observation_models::LinkEnds angularPositionLinkEnds,
                               const simulation_setup::NamedBodyMap& bodyMap,
                               const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )

{
    boost::shared_ptr< AngularPositionScaling > angularPositionScaling =
            boost::make_shared< AngularPositionScaling >( );

    SingleLinkObservationTwoPartialList angularPositionPartials;

    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( boost::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                    initialDynamicalParameters.at( i ) ) == NULL )
        {
            std::cerr<<"Error when making one way range partials, could not identify parameter"<<std::endl;
        }

        std::string acceleratedBody = boost::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                    initialDynamicalParameters.at( i ) )->getParameterName( ).second.first;

        boost::shared_ptr< AngularPositionPartial > currentRangePartial = createAngularPositionPartialWrtInitialPosition(
                    angularPositionLinkEnds, bodyMap, acceleratedBody, angularPositionScaling );

        if( currentRangePartial != NULL )
        {
            currentPair = std::pair< int, int >( currentIndex, 3 );
            angularPositionPartials[ currentPair ] = currentRangePartial;
        }

        currentIndex += 6;
    }

    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator parameterIterator =
         doubleParametersToEstimate.begin( ); parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        boost::shared_ptr< AngularPositionPartial > currentRangePartial = createAngularPositionPartialWrtParameter(
                    angularPositionLinkEnds, bodyMap, parameterIterator->second, angularPositionScaling );

        if( currentRangePartial != NULL )
        {
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            angularPositionPartials[ currentPair ] = currentRangePartial;
        }
    }

    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParametersToEstimate =
            parametersToEstimate->getVectorParameters( );
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        boost::shared_ptr< AngularPositionPartial > currentRangePartial = createAngularPositionPartialWrtParameter(
                    angularPositionLinkEnds, bodyMap, parameterIterator->second, angularPositionScaling );

        if( currentRangePartial != NULL )
        {
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            angularPositionPartials[ currentPair ] = currentRangePartial;
        }

    }
    return std::make_pair( angularPositionPartials, angularPositionScaling );
}

template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationTwoPartialList, boost::shared_ptr< PositionPartialScaling > > >
createAngularPositionPartials( const std::vector< observation_models::LinkEnds > linkEnds,
                               const simulation_setup::NamedBodyMap& bodyMap,
                               const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
{
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationTwoPartialList , boost::shared_ptr< PositionPartialScaling > > > angularPositionPartials;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        assert( linkEnds[ i ].count( receiver ) > 0  );
        assert( linkEnds[ i ].count( transmitter ) > 0  );

        angularPositionPartials[ linkEnds[ i ] ] = createAngularPositionPartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate );
    }
    return angularPositionPartials;
}

}

}

#endif // CREATEANGULARPOSITIONPARTIALS_H
