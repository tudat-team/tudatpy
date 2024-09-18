/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEPOSITIONPARTIALS_H
#define TUDAT_CREATEPOSITIONPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"


namespace tudat
{

namespace observation_partials
{

//! Function to return partial(s) of position of reference point w.r.t state of a single body.
/*!
 *  Function to return partial(s) of position of reference point w.r.t state of a single body. A set of link ends and the
 *  name  of the body wrt the position of which the partials are to be created. A map is returned, with the LinkEndType as
 *  key and  pointer to state partial as value. An entry for the map is created for each link end which corresponds to
 *  the body wrt the position of which the partial is to be taken.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether the body corresponds to the
 *  requested body and, if so, a partial object is created.
 *  \param bodies Map of body objects, used in the creation of the partials.
 *  \param bodyToEstimate Name of body wrt the position of which partials are to be created.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtBodyState(
        const observation_models::LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate );

//! Function to return partial(s) of position of ground station(s) w.r.t. rotational state of a single body.
/*!
 *  Function to return partial(s) of position of reference point w.r.t rotational state of a single body. A set of link ends
 *  and the name of the body wrt the position of which the partials are to be created. A map is returned, with the LinkEndType as
 *  key and pointer to state partial as value. An entry for the map is created for each link end which corresponds to
 *  the body wrt the position of which the partial is to be taken.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether the body corresponds to the
 *  requested body and, if so, a partial object is created.
 *  \param bodies Map of body objects, used in the creation of the partials.
 *  \param bodyToEstimate Name of body wrt the position of which partials are to be created.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtBodyRotationalState(
        const observation_models::LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToEstimate );

//! Function to return partial object(s) of position of reference point w.r.t. a (double) parameter.
/*!
 *  Function to return partial object(s) of position of reference point w.r.t. a (double) parameter. A set of link ends and
 *  parameter object  wrt which the partials are to be created. A map is returned, with the LinkEndType as key and pointer
 *  to position partial as value. An entry for the map is created for each entry of linkEnds for which there is a direct
 *  dependency between its position and  the parameter in question.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether a there is a direct dependency
 *  between its position and the parameter in question.
 *  \param bodies Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Map of position partial objects, one entry for each link end corresponding to the parameterToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate );

//! Function to return partial object(s) of position of reference point w.r.t. a (vector) parameter.
/*!
 *  Function to return partial object(s) of position of reference point w.r.t. a (vector) parameter. A set of link ends and
 *  parameter object  wrt which the partials are to be created. A map is returned, with the LinkEndType as key and pointer
 *  to position partial as value. An entry for the map is created for each entry of linkEnds for which there is a direct
 *  dependency between its position and  the parameter in question.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether a there is a direct dependency
 *  between its position and the parameter in question.
 *  \param bodies Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Map of position partial objects, one entry for each link end corresponding to the parameterToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate );

//! Function to create partial object(s) of rotation matrix wrt translational state
/*!
 *  Function to create partial object(s) of rotation matrix wrt a state parameter
 *  \param currentBody Body for which partial is to be created (must have a synchronous rotation model for output to be non-null)
 *  \return Rotation matrix partial object
 */
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtTranslationalState(
        const std::shared_ptr< simulation_setup::Body > currentBody );

//! Function to create partial object(s) of rotation matrix wrt a state parameter
/*!
 *  Function to create partial object(s) of rotation matrix wrt a state parameter
 *  \param bodies Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Rotation matrix partial object
 */
template< typename InitialStateParameterType >
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtStateParameter(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix<
        InitialStateParameterType, Eigen::Dynamic, 1 > > > parameterToEstimate )
{
    using namespace simulation_setup;
    using namespace ephemerides;

    // Declare return object.
    std::shared_ptr< RotationMatrixPartial >  rotationMatrixPartial;

    // Get body for rotation of which partial is to be created.
    std::shared_ptr< Body > currentBody = bodies.at( parameterToEstimate->getParameterName( ).second.first );

    // Check for which rotation model parameter the partial object is to be created.
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::initial_rotational_body_state:

        // Create rotation matrix partial object
        rotationMatrixPartial = std::make_shared< RotationMatrixPartialWrtQuaternion >(
                    std::bind( &Body::getCurrentRotationToGlobalFrame, currentBody ) );
        break;

    default:
        std::string errorMessage = "Warning, rotation matrix partial not implemented for state parameter " +
                std::to_string( parameterToEstimate->getParameterName( ).first );
        throw std::runtime_error( errorMessage );
        break;
    }

    return rotationMatrixPartial;

}

//! Function to create partial object(s) of rotation matrix wrt a (double) parameter.
/*!
 *  Function to create partial object(s) of rotation matrix from a body-fixed to inertial frame wrt a (double) parameter.
 *  \param bodies Map of body objects, used in the creation of the partial.
 *  \param parameterToEstimate Parameter wrt which rotation matrix partial object is to be created.
 *  \return Requested rotation matrix partial object.
 */
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate );

//! Function to create partial object(s) of rotation matrix wrt a (vector) parameter.
/*!
 *  Function to create partial object(s) of rotation matrix from a body-fixed to inertial frame wrt a (vector) parameter.
 *  \param bodies Map of body objects, used in the creation of the partial.
 *  \param parameterToEstimate Parameter wrt which rotation matrix partial object is to be created.
 *  \return Requested rotation matrix partial object.
 */
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate );

//! Function to create rotation matrix partial objects for all rotation model parameters of given body in given set of
//! parameters.
/*!
 *  Function to create rotation matrix partial objects for all rotation model parameters of given body in given set of
 *  parameters. All parameter objects of given EstimatableParameterSet are checked whether they are firstly a property of
 *  the requested body, and secondly if they represent a property of the rotation from the body-fixed to base frame.
 *  If both conditions are met, a partial  is created and added to the return list.
 *  \param parametersToEstimate Set of all parameters which are to be estimated.
 *  \param bodyName Name of body for which to create partials of rotation from body-fixed to base frame
 *  \param bodies Map of body objects, used in the creation of the partials.
 *  \return Rotation matrix partial objects for all rotation model parameters of bodyName in parametersToEstimate.
 */
template< typename ParameterType >
RotationMatrixPartialNamedList createRotationMatrixPartials(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::string& bodyName, const simulation_setup::SystemOfBodies& bodies )

{
    //std::vector< std::shared_ptr< EstimatableParameter< Eigen:: > > >
    //getEstimatedInitialStateParameters( )

    // Declare map to return.
    std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
            std::shared_ptr< RotationMatrixPartial > >
            rotationMatrixPartials;

    // Check if parametersToEstimate pointer isn't empty
    if ( !parametersToEstimate )
    {
        throw std::runtime_error(
                "Error when creating rotation matrix partials, the provided parameters to estimate are empty.");
    }
    // Retrieve double and vector parameters from total set of parameters.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > stateParameters =
            parametersToEstimate->getInitialStateParameters( );
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
            parametersToEstimate->getDoubleParameters( );
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getVectorParameters( );

    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialWrtTranslationalState =
            createRotationMatrixPartialsWrtTranslationalState( bodies.at( bodyName ) );
    if( rotationMatrixPartialWrtTranslationalState != nullptr )
    {
        rotationMatrixPartials[ std::make_pair(
                    estimatable_parameters::initial_body_state, "" ) ] = rotationMatrixPartialWrtTranslationalState;
    }

    for( auto parameterIterator = stateParameters.begin( ); parameterIterator != stateParameters.end( );
         parameterIterator++ )
    {
        // Check parameter is rotational property of requested body.
        if( ( parameterIterator->second->getParameterName( ).second.first == bodyName ) &&
                ( estimatable_parameters::isParameterRotationMatrixProperty(
                      parameterIterator->second->getParameterName( ).first ) ) )
        {
            // Create partial object.
            rotationMatrixPartials[ std::make_pair(
                        parameterIterator->second->getParameterName( ).first,
                        parameterIterator->second->getSecondaryIdentifier( )) ] =
                    createRotationMatrixPartialsWrtStateParameter( bodies, parameterIterator->second );
        }
    }

    // Iterate over double parameters.
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator = doubleParameters.begin( ); parameterIterator != doubleParameters.end( ); parameterIterator++ )
    {
        // Check parameter is rotational property of requested body.
        if( ( parameterIterator->second->getParameterName( ).second.first == bodyName ) &&
                ( estimatable_parameters::isParameterRotationMatrixProperty(
                      parameterIterator->second->getParameterName( ).first ) ) )
        {
            // Create partial object.
            rotationMatrixPartials[ std::make_pair(
                        parameterIterator->second->getParameterName( ).first,
                        parameterIterator->second->getSecondaryIdentifier( )) ] =
                    createRotationMatrixPartialsWrtParameter( bodies, parameterIterator->second );
        }
    }

    // Iterate over vector parameters.
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >::iterator
         parameterIterator = vectorParameters.begin( ); parameterIterator != vectorParameters.end( ); parameterIterator++ )
    {
        // Check parameter is rotational property of requested body.
        if( ( parameterIterator->second->getParameterName( ).second.first == bodyName ) &&
                ( estimatable_parameters::isParameterRotationMatrixProperty(
                      parameterIterator->second->getParameterName( ).first ) ) )
        {
            // Create partial object.
            rotationMatrixPartials[ std::make_pair(
                        parameterIterator->second->getParameterName( ).first,
                        parameterIterator->second->getSecondaryIdentifier( )) ] =
                    createRotationMatrixPartialsWrtParameter( bodies, parameterIterator->second );
        }
    }

    return rotationMatrixPartials;
}

//! Typedef of list of RotationMatrixPartial objects, ordered by parameter.
typedef std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
std::shared_ptr< RotationMatrixPartial > > RotationMatrixPartialNamedList;


}

}

#endif // TUDAT_CREATEPOSITIONPARTIALS_H
