#ifndef TUDAT_CREATEPOSITIONPARTIALS_H
#define TUDAT_CREATEPOSITIONPARTIALS_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

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

//! Typedef of list of RotationMatrixPartial objects, ordered by parameter.
typedef std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
boost::shared_ptr< RotationMatrixPartial > > RotationMatrixPartialNamedList;

//! Function to return partial(s) of position of ground station(s) w.r.t. state of a single body.
/*!
 *  Function to return partial(s) of position of ground station(s) w.r.t. state of a single body. A set of link ends and the name
 *  of the body wrt the position of which the partials are to be created. A map is returned, with the LinkEndType as key and
 *  pointer to position partial as value. An entry for the map is created for each link end which corresponds to the body wrt the position of
 *  which the partial is to be taken.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether the body corresponds to the requested body and, if so,
 *  a partial object is created.
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param bodyToEstimate Name of body wrt the position of which partials are to be created.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */
std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > createPositionPartialsWrtBodyPosition(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate );

//! Function to return partial object(s) of position of ground station(s) w.r.t. a (double) parameter.
/*!
 *  Function to return partial object(s) of position of ground station(s) w.r.t. a (double) parameter. A set of link ends and parameter object
 *  wrt which the partials are to be created. A map is returned, with the LinkEndType as key and pointer to position partial as value.
 *  An entry for the map is created for each entry of linkEnds for which there is a direct dependency between its position and
 *  the parameter in question.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether a there is a direct dependency between its position and
 *  the parameter in question.
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */

std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > createPositionPartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate );

//! Function to return partial object(s) of position of ground station(s) w.r.t. a (vector) parameter.
/*!
 *  Function to return partial object(s) of position of ground station(s) w.r.t. a (vector) parameter. A set of link ends and parameter object
 *  wrt which the partials are to be created. A map is returned, with the LinkEndType as key and pointer to position partial as value.
 *  An entry for the map is created for each entry of linkEnds for which there is a direct dependency between its position and
 *  the parameter in question.
 *  \param linkEnds Set of link ends, for each entry of this map, it is checked whether a there is a direct dependency between its position and
 *  the parameter in question.
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \param parameterToEstimate Parameter object wrt which partials are to be calculated.
 *  \return Map of position partial objects, one entry for each link end corresponding to the bodyToEstimate.
 */
std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > createPositionPartialsWrtParameter(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate );

//! Function to create partial object(s) of rotation matrix wrt a (double) parameter.
/*!
 *  Function to create partial object(s) of rotation matrix from a body-fixed to inertial frame wrt a (double) parameter.
 *  \param bodyMap Map of body objects, used in the creation of the partial.
 *  \param parameterToEstimate Parameter wrt which rotation matrix partial object is to be created.
 *  \return Requested rotation matrix partial object.
 */
boost::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate );

//! Function to create partial object(s) of rotation matrix wrt a (vector) parameter.
/*!
 *  Function to create partial object(s) of rotation matrix from a body-fixed to inertial frame wrt a (vector) parameter.
 *  \param bodyMap Map of body objects, used in the creation of the partial.
 *  \param parameterToEstimate Parameter wrt which rotation matrix partial object is to be created.
 *  \return Requested rotation matrix partial object.
 */
boost::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate );

//! Function to create rotation matrix partial objects for all rotation model parameters of given body in given set of
//! parameters.
/*!
 *  Function to create rotation matrix partial objects for all rotation model parameters of given body in given set of
 *  parameters. All parameter objects of given EstimatableParameterSet are checked whether they are firstly a property of
 *  the requested body, and secondly if they represent a property of the rotation from the body-fixed to base frame.
 *  If both conditions are met, a partial  is created and added to teh return list.
 *  \param parametersToEstimate Set of all parameters which are to be estimated.
 *  \param bodyName Name of body for which to create partials of rotation from body-fixed to base frame
 *  \param bodyMap Map of body objects, used in the creation of the partials.
 *  \return Rotation matrix partial objects for all rotation model parameters of bodyName in parametersToEstimate.
 */
template< typename ParameterType >
RotationMatrixPartialNamedList createRotationMatrixPartials(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::string& bodyName, const simulation_setup::NamedBodyMap& bodyMap )

{
    // Declare map to return.
    std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
            boost::shared_ptr< RotationMatrixPartial > >
            rotationMatrixPartials;

    // Retrieve double and vector parameters from total set of parameters.
    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
            parametersToEstimate->getDoubleParameters( );
    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getVectorParameters( );

    // Iterate over double parameters.
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
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
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >::iterator
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

boost::shared_ptr< PositionObervationPartial > createPositionObservablePartialWrtPosition(
        const observation_models::LinkEnds linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const boost::shared_ptr< PositionObservationScaling > positionObservableScaler );

template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 3 > > >, boost::shared_ptr< PositionPartialScaling > >
createPositionObservablePartials( const observation_models::LinkEnds positionObservableLinkEnds,
                                  const simulation_setup::NamedBodyMap& bodyMap,
                                  const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )

{
    boost::shared_ptr< PositionObservationScaling > positionObservableScaling =
            boost::make_shared< PositionObservationScaling >( );

    SingleLinkObservationThreePartialList positionObservablePartials;

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
            std::cerr<<"Error when making position observable partials, could not identify parameter"<<std::endl;
        }

        std::string acceleratedBody = boost::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                    initialDynamicalParameters.at( i ) )->getParameterName( ).second.first;

        boost::shared_ptr< PositionObervationPartial > currentObservablePartial = createPositionObservablePartialWrtPosition(
                    positionObservableLinkEnds, bodyMap, acceleratedBody, positionObservableScaling );

        if( currentObservablePartial != NULL )
        {
            currentPair = std::pair< int, int >( currentIndex, 3 );
            positionObservablePartials[ currentPair ] = currentObservablePartial;
        }

        currentIndex += 6;
    }

    return std::make_pair( positionObservablePartials, positionObservableScaling );
}

template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationThreePartialList, boost::shared_ptr< PositionPartialScaling > > >
createPositionObservablePartials( const std::vector<  observation_models::LinkEnds > linkEnds,
                                  const simulation_setup::NamedBodyMap& bodyMap,
                                  const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
{
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationThreePartialList , boost::shared_ptr< PositionPartialScaling > > > positionObservablePartials;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        if( linkEnds[ i ].count( observation_models::observed_body ) == 0 || linkEnds[ i ].size( ) != 1 )
        {
            std::cerr<<"Error when making position observable partial, link ends are wrong"<<std::endl;
        }

        positionObservablePartials[ linkEnds[ i ] ] = createPositionObservablePartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate );
    }
    return positionObservablePartials;
}

}

}

#endif // CREATEPOSITIONPARTIALS_H
