/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"


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
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param bodyToEstimate Name of body wrt the position of which partials are to be created.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtBodyState(
        const observation_models::LinkEnds& linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate );

//! Function to return partial(s) of position of ground station(s) w.r.t. rotational state of a single body.
/*!
 *  Function to return partial(s) of position of reference point w.r.t rotational state of a single body. A set of link ends
 *  and the name of the body wrt the position of which the partials are to be created. A map is returned, with the LinkEndType as
 *  key and pointer to state partial as value. An entry for the map is created for each link end which corresponds to
 *  the body wrt the position of which the partial is to be taken.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether the body corresponds to the
 *  requested body and, if so, a partial object is created.
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param bodyToEstimate Name of body wrt the position of which partials are to be created.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtBodyRotationalState(
        const observation_models::LinkEnds& linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string& bodyToEstimate );

//! Function to return partial object(s) of position of reference point w.r.t. a (double) parameter.
/*!
 *  Function to return partial object(s) of position of reference point w.r.t. a (double) parameter. A set of link ends and
 *  parameter object  wrt which the partials are to be created. A map is returned, with the LinkEndType as key and pointer
 *  to position partial as value. An entry for the map is created for each entry of linkEnds for which there is a direct
 *  dependency between its position and  the parameter in question.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether a there is a direct dependency
 *  between its position and the parameter in question.
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Map of position partial objects, one entry for each link end corresponding to the parameterToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate );

//! Function to return partial object(s) of position of reference point w.r.t. a (vector) parameter.
/*!
 *  Function to return partial object(s) of position of reference point w.r.t. a (vector) parameter. A set of link ends and
 *  parameter object  wrt which the partials are to be created. A map is returned, with the LinkEndType as key and pointer
 *  to position partial as value. An entry for the map is created for each entry of linkEnds for which there is a direct
 *  dependency between its position and  the parameter in question.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether a there is a direct dependency
 *  between its position and the parameter in question.
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Map of position partial objects, one entry for each link end corresponding to the parameterToEstimate.
 */
std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > createCartesianStatePartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate );



//! Function to create partial object(s) of rotation matrix wrt a state parameter
/*!
 *  Function to create partial object(s) of rotation matrix wrt a state parameter
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Rotation matrix partial object
 */
template< typename InitialStateParameterType >
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtStateParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix<
        InitialStateParameterType, Eigen::Dynamic, 1 > > > parameterToEstimate )
{
    using namespace simulation_setup;
    using namespace ephemerides;

    // Declare return object.
    std::shared_ptr< RotationMatrixPartial >  rotationMatrixPartial;

    // Get body for rotation of which partial is to be created.
    std::shared_ptr< Body > currentBody = bodyMap.at( parameterToEstimate->getParameterName( ).second.first );

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
 *  \param bodyMap Map of body objects, used in the creation of the partial.
 *  \param parameterToEstimate Parameter wrt which rotation matrix partial object is to be created.
 *  \return Requested rotation matrix partial object.
 */
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate );

//! Function to create partial object(s) of rotation matrix wrt a (vector) parameter.
/*!
 *  Function to create partial object(s) of rotation matrix from a body-fixed to inertial frame wrt a (vector) parameter.
 *  \param bodyMap Map of body objects, used in the creation of the partial.
 *  \param parameterToEstimate Parameter wrt which rotation matrix partial object is to be created.
 *  \return Requested rotation matrix partial object.
 */
std::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
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
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \return Rotation matrix partial objects for all rotation model parameters of bodyName in parametersToEstimate.
 */
template< typename ParameterType >
RotationMatrixPartialNamedList createRotationMatrixPartials(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::string& bodyName, const simulation_setup::NamedBodyMap& bodyMap )

{
    //std::vector< std::shared_ptr< EstimatableParameter< Eigen:: > > >
    //getEstimatedInitialStateParameters( )

    // Declare map to return.
    std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
            std::shared_ptr< RotationMatrixPartial > >
            rotationMatrixPartials;

    // Retrieve double and vector parameters from total set of parameters.
   std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > stateParameters =
            parametersToEstimate->getInitialStateParameters( );
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
            parametersToEstimate->getDoubleParameters( );
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getVectorParameters( );

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
                        parameterIterator->second->getSecondaryIdentifier( ) ) ] =
                    createRotationMatrixPartialsWrtStateParameter( bodyMap, parameterIterator->second );
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
                        parameterIterator->second->getSecondaryIdentifier( ) ) ] =
                    createRotationMatrixPartialsWrtParameter( bodyMap, parameterIterator->second );
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
                        parameterIterator->second->getSecondaryIdentifier( ) ) ] =
                    createRotationMatrixPartialsWrtParameter( bodyMap, parameterIterator->second );
        }
    }

    return rotationMatrixPartials;
}

//! Typedef of list of RotationMatrixPartial objects, ordered by parameter.
typedef std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
std::shared_ptr< RotationMatrixPartial > > RotationMatrixPartialNamedList;

//! Function to create an objects that computes the partial derivatives of a three-dimensional position observable w.r.t.
//! the position of a body.
/*!
 *  Function to create an objects that computes the partial derivatives of a three-dimensional position observable w.r.t.
 *  the position of a body.
 *  \param linkEnds Set of link ends used for observation model of three-dimensional position
 *  \param bodyMap List of bodies that comprise the environment
 *  \param bodyToEstimate Name of body w.r.t. the position of which a partial is to be compured
 *  \param positionObservableScaler Object that scales position partial to observable partial.
 *  \return Single object that computes partial of given observable w.r.t. given parameter.
 */
std::shared_ptr< PositionObervationPartial > createPositionObservablePartialWrtPosition(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< PositionObservationScaling > positionObservableScaler );

//! Function to create a list of objects that compute the partial derivatives of a three-dimensional position observable.
/*!
 *  Function to create a list of objects that compute the partial derivatives of a three-dimensional position observable
 *  A single object is created for each parameter w.r.t. whih a partial derivative is to be taken. Note that the
 *  three-dimensional position observable is only sensitive to the position of the body under observation.
 *  \param positionObservableLinkEnds Set of link ends used for observation model of three-dimensional position
 *  \param bodyMap List of bodies that comprise the environment
 *  \param parametersToEstimate List of parameters that is to be estimated.
 *  \return Pair, first entry is map (key is start index and size of parameter; value is partial object), secod entry is
 *  scaling object to be used for all partials.
 */
template< typename ParameterType >
std::pair< SingleLinkObservationThreePartialList, std::shared_ptr< PositionPartialScaling > >
createPositionObservablePartials(
        const observation_models::LinkEnds positionObservableLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )

{
    // Create scaling object, to be used for each partial created here (i.e. same scaling for different parameters but same
    // observable).
    std::shared_ptr< PositionObservationScaling > positionObservableScaling =
            std::make_shared< PositionObservationScaling >( );

    SingleLinkObservationThreePartialList positionObservablePartials;

    // Define start index and size of current parameter under consideration
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( std::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                    initialDynamicalParameters.at( i ) ) != nullptr )
        {

            // Retrieve name of body w.r.t. position of which partial is to be taken.
            std::string acceleratedBody = std::dynamic_pointer_cast<
                    estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                        initialDynamicalParameters.at( i ) )->getParameterName( ).second.first;

            // Create partial (if needed)
            std::shared_ptr< PositionObervationPartial > currentObservablePartial =
                    createPositionObservablePartialWrtPosition(
                        positionObservableLinkEnds, bodyMap, acceleratedBody, positionObservableScaling );

            // If partial exists, then dependency exists and parameter must be added.
            if( currentObservablePartial != nullptr )
            {
                currentPair = std::pair< int, int >( currentIndex, 6 );
                positionObservablePartials[ currentPair ] = currentObservablePartial;
            }

            currentIndex += 6;
        }
    }

    return std::make_pair( positionObservablePartials, positionObservableScaling );
}

//! Function to create a list of objects that compute the partial derivatives of a list of 3-dimensional position observable.
/*!
 *  Function to create a list of objects that compute the partial derivatives of a list of 3-dimensional position observable
 *  A single object is created for each parameter w.r.t. whih a partial derivative is to be taken, separately for each set of
 *  link ends. Note that the three-dimensional position observable is only sensitive to the position of the body under
 * observation.
 *  \param linkEnds List of sets of link ends used for observation models of three-dimensional position
 *  \param bodyMap List of bodies that comprise the environment
 *  \param parametersToEstimate List of parameters that is to be estimated.
 *  \return For each set of link ends a single pair, containing:
 *  First entry is map (key is start index and size of parameter; value is partial object), secod entry is
 *  scaling object to be used for all partials.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds,
std::pair< SingleLinkObservationThreePartialList, std::shared_ptr< PositionPartialScaling > > >
createPositionObservablePartials(
        const std::vector<  observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
{
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationThreePartialList,
            std::shared_ptr< PositionPartialScaling > > > positionObservablePartials;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        if( linkEnds[ i ].count( observation_models::observed_body ) == 0 || linkEnds[ i ].size( ) != 1 )
        {
            throw std::runtime_error( "Error when making position observable partial, link ends are wrong" );
        }

        positionObservablePartials[ linkEnds[ i ] ] = createPositionObservablePartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate );
    }
    return positionObservablePartials;
}

}

}

#endif // TUDAT_CREATEPOSITIONPARTIALS_H
