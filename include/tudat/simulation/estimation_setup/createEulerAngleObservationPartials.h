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

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialRotationalState.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/eulerAngleObservablePartials.h"
#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

class SystemOfBodies;

}  // namespace simulation_setup

namespace observation_partials
{

//! Function to create object that comptes partial derivative of Euler angle observation w.r.t. current rotational state
/*!
 * Function to create object that comptes partial derivative of Euler angle observation w.r.t. current rotational state
 * \param parameterIdentifier Type and reference body of parameter to be estimated (initial rotational state)
 * \param rotationToBaseFrameFunction Function returning the rotation from body-fixed to base frame as a function of time
 * \return Object that comptes partial derivative of Euler angle observation w.r.t. current rotational state
 */
std::shared_ptr< ObservationPartial< 3 > > createEulerAngleObservablePartialWrtCurrentRotationalState(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
        const std::function< Eigen::Quaterniond( const double ) > rotationToBaseFrameFunction );

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
std::pair< SingleLinkObservationThreePartialList, std::shared_ptr< PositionPartialScaling > > createEulerAngleObservablePartials(
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

    // Rotation of the observed body, and partials of its rotation matrix w.r.t. the estimated parameters
    const std::string observedBody = eulerAngleLinkEnds.at( observation_models::observed_body ).bodyName_;
    const std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel = bodies.at( observedBody )->getRotationalEphemeris( );
    const std::function< Eigen::Quaterniond( const double ) > rotationToBaseFrameFunction =
            std::bind( &ephemerides::RotationalEphemeris::getRotationToBaseFrameTemplated< double >, rotationModel, std::placeholders::_1 );
    const RotationMatrixPartialNamedList rotationMatrixPartials =
            createRotationMatrixPartials( parametersToEstimate, observedBody, bodies );

    // Iterate over the estimated initial states
    int currentIndex = 0;
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParameters = parametersToEstimate->getEstimatedInitialStateParameters( );
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        const estimatable_parameters::EstimatebleParameterIdentifier parameterName =
                initialDynamicalParameters.at( i )->getParameterName( );
        if( parameterName.first == estimatable_parameters::initial_rotational_body_state && parameterName.second.first == observedBody )
        {
            eulerAnglePartials[ std::make_pair( currentIndex, 7 ) ] =
                    createEulerAngleObservablePartialWrtCurrentRotationalState( parameterName, rotationToBaseFrameFunction );
        }
        else if( ( parameterName.first == estimatable_parameters::initial_body_state ||
                   parameterName.first == estimatable_parameters::arc_wise_initial_body_state ) &&
                 rotationMatrixPartials.count( std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
        {
            // Rotation model that depends on the state of the observed body w.r.t. its central body
            const std::string centralBody =
                    std::dynamic_pointer_cast< ephemerides::SynchronousRotationalEphemeris >( rotationModel )->getCentralBodyName( );
            if( parameterName.second.first == observedBody || parameterName.second.first == centralBody )
            {
                eulerAnglePartials[ std::make_pair( currentIndex, 6 ) ] =
                        std::make_shared< EulerAngleObservationPartialWrtRotationParameter >(
                                parameterName,
                                rotationMatrixPartials.at( std::make_pair( estimatable_parameters::initial_body_state, "" ) ),
                                rotationToBaseFrameFunction,
                                parameterName.second.first == observedBody ? 1.0 : -1.0 );
            }
        }
        currentIndex += initialDynamicalParameters.at( i )->getParameterSize( );
    }

    // Iterate over the estimated rotation model parameters of the observed body
    for( auto parameterIterator : parametersToEstimate->getDoubleParameters( ) )
    {
        const estimatable_parameters::EstimatebleParameterIdentifier parameterName = parameterIterator.second->getParameterName( );
        if( parameterName.second.first == observedBody && estimatable_parameters::isParameterRotationMatrixProperty( parameterName.first ) )
        {
            eulerAnglePartials[ std::make_pair( parameterIterator.first, 1 ) ] =
                    std::make_shared< EulerAngleObservationPartialWrtRotationParameter >(
                            parameterName,
                            rotationMatrixPartials.at(
                                    std::make_pair( parameterName.first, parameterIterator.second->getSecondaryIdentifier( ) ) ),
                            rotationToBaseFrameFunction );
        }
    }
    for( auto parameterIterator : parametersToEstimate->getVectorParameters( ) )
    {
        const estimatable_parameters::EstimatebleParameterIdentifier parameterName = parameterIterator.second->getParameterName( );
        if( parameterName.second.first == observedBody && estimatable_parameters::isParameterRotationMatrixProperty( parameterName.first ) )
        {
            eulerAnglePartials[ std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) ] =
                    std::make_shared< EulerAngleObservationPartialWrtRotationParameter >(
                            parameterName,
                            rotationMatrixPartials.at(
                                    std::make_pair( parameterName.first, parameterIterator.second->getSecondaryIdentifier( ) ) ),
                            rotationToBaseFrameFunction );
        }
    }

    return std::make_pair( eulerAnglePartials, eulerAngleScaling );
}

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_CREATEEULERANGLEOBSERVATIONPARTIALS_H
