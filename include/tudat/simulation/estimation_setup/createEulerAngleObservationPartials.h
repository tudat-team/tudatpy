/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialRotationalState.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"


namespace tudat
{

namespace observation_partials
{

//! Function to create object that comptes partial derivative of Euler angle observation w.r.t. current rotational state
/*!
 * Function to create object that comptes partial derivative of Euler angle observation w.r.t. current rotational state
 * \param parameterIdentifier Type and reference body of parameter to be estimated (initial rotational state)
 * \return Object that comptes partial derivative of Euler angle observation w.r.t. current rotational state
 */
std::shared_ptr< ObservationPartial< 3 > > createEulerAngleObservablePartialWrtCurrentRotationalState(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier );

//! Function to compute Euler angle observation partial objects for a single set of link ends
/*!
 *  Function to compute Euler angle observation partial objects for a single set of link ends
 *  \param eulerAngleLinkEnds Link ends (observed_body only) for which Euler angle partials are to be calculated
 *  \param bodies List of all bodies, for creating Euler angle partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary Euler angle partials of a single link end, and a nulptr position scaling pointer.
 */
template< typename ParameterType >
std::pair< SingleLinkObservationThreePartialList, std::shared_ptr< PositionPartialScaling > >
createEulerAngleObservablePartials(
        const observation_models::LinkEnds eulerAngleLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable = false,
        const bool isPartialForConcatenatedObservable = false )

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
                    createEulerAngleObservablePartialWrtCurrentRotationalState(
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

    // Iterate over all double parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator = doubleParametersToEstimate.begin( );
         parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameterIterator->second->getParameterName( ).first ) )
        {
            throw std::runtime_error(
                        "Error when making Euler angle partial, found kinematic rotation parameter: " +
                        std::to_string( estimatable_parameters::isParameterRotationMatrixProperty(
                                            parameterIterator->second->getParameterName( ).first ) ) +
                        ". Computation requires derivatives of Euler angles w.r.t. rotation matrix entries, which is not yet implemented." );
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate = parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator = vectorParametersToEstimate.begin( );
         parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameterIterator->second->getParameterName( ).first ) )
        {
            throw std::runtime_error(
                        "Error when making Euler angle partial, found kinematic rotation parameter: " +
                        std::to_string( estimatable_parameters::isParameterRotationMatrixProperty(
                                            parameterIterator->second->getParameterName( ).first ) ) +
                        ". Computation requires derivatives of Euler angles w.r.t. rotation matrix entries, which is not yet implemented." );
        }
    }

    return std::make_pair( eulerAnglePartials, eulerAngleScaling );
}

}

}

#endif // TUDAT_CREATEEULERANGLEOBSERVATIONPARTIALS_H

