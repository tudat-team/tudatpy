/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEEULERANGLEOBSERVATIONPARTIALS_H
#define TUDAT_CREATEEULERANGLEOBSERVATIONPARTIALS_H


#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialRotationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"


namespace tudat
{

namespace observation_partials
{


std::shared_ptr< ObservationPartial< 3 > > createEulerAngleObservablePartialWrtCurrentOrientation(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier );

template< typename ParameterType >
std::pair< SingleLinkObservationThreePartialList, std::shared_ptr< PositionPartialScaling > >
createEulerAngleObservablePartials(
        const observation_models::LinkEnds eulerAngleLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )

{
    // Create scaling object, to be used for each partial created here (i.e. same scaling for different parameters but same
    // observable).
    std::shared_ptr< PositionPartialScaling > eulerAngleScaling = nullptr;

    SingleLinkObservationThreePartialList eulerAnglePartials;

    // Define start index and size of current parameter under consideration
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( std::dynamic_pointer_cast< estimatable_parameters::InitialRotationalStateParameter< ParameterType > >(
                    initialDynamicalParameters.at( i ) ) != nullptr )
        {
            // Create partial (if needed)
            std::shared_ptr< ObservationPartial< 3 > > currentObservablePartial =
                    createEulerAngleObservablePartialWrtCurrentOrientation(
                       initialDynamicalParameters.at( i )->getParameterName( ) );

            // If partial exists, then dependency exists and parameter must be added.
            if( currentObservablePartial != nullptr )
            {
                currentPair = std::pair< int, int >( currentIndex, 7 );
                eulerAnglePartials[ currentPair ] = currentObservablePartial;
            }
        }
        currentIndex += initialDynamicalParameters.at( i )->getParameterSize( );

    }

    return std::make_pair( eulerAnglePartials, eulerAngleScaling );
}

template< typename ParameterType >
std::map< observation_models::LinkEnds,
std::pair< SingleLinkObservationThreePartialList, std::shared_ptr< PositionPartialScaling > > >
createEulerAngleObservablePartials(
        const std::vector<  observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
{
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationThreePartialList,
            std::shared_ptr< PositionPartialScaling > > > eulerAnlgeObservablePartials;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        if( linkEnds[ i ].count( observation_models::observed_body ) == 0 || linkEnds[ i ].size( ) != 1 )
        {
            throw std::runtime_error( "Error when making position observable partial, link ends are wrong" );
        }

        eulerAnlgeObservablePartials[ linkEnds[ i ] ] = createEulerAngleObservablePartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate );
    }
    return eulerAnlgeObservablePartials;
}


}

}

#endif // TUDAT_CREATEEULERANGLEOBSERVATIONPARTIALS_H

