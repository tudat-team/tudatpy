/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"


namespace tudat
{

namespace observation_partials
{


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


} // namespace observation_partials

} // namespace tudat

#endif // TUDAT_CREATEPOSITIONPARTIALS_H
