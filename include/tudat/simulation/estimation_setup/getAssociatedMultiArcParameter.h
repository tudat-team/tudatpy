/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GETASSOCIATEDMULTIARCPARAMETER_H
#define TUDAT_GETASSOCIATEDMULTIARCPARAMETER_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to get the multi-arc parameter equivalent of a single-arc initial state parameter
/*!
 *  Function to get the multi-arc parameter equivalent of a single-arc initial state parameter. The initial state arcs are
 *  provided as input to this function.
 *  \param singleArcParameter Single-arc parameter object for which the multi-arc equivalent is to be created.
 *  \param arcStartTimes Vector of start times for separate arcs.
 *  \return Multi-arc parameter equivalent of single-arc initial state parameter input
 */
template< typename StateScalarType >
std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
getAssociatedMultiArcParameter(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
                singleArcParameter,
        const std::vector< double >& arcStartTimes )
{
    std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
            multiArcParameter;

    // Check state type
    switch( singleArcParameter->getParameterName( ).first )
    {
        case estimatable_parameters::initial_body_state: {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::InitialTranslationalStateParameter< StateScalarType > >
                    singleArcTranslationalStateParameter =
                            std::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< StateScalarType > >(
                                    singleArcParameter );
            if( singleArcTranslationalStateParameter == nullptr )
            {
                throw std::runtime_error(
                        "Error when getting multi-arc parameter from single-arc equivalent, single-arc translational state is "
                        "inconsistent " );
            }

            // Retrieve single-arc initial state
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > singleArcInitialState =
                    singleArcTranslationalStateParameter->getParameterValue( );

            // Create multi-arc initial states. First arc initial state is taken from single-arc, other initial states set to zero.
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > multiArcInitialStates =
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 6 * arcStartTimes.size( ) );
            multiArcInitialStates.segment( 0, 6 ) = singleArcInitialState;

            // Creater multi-arc parameter
            std::vector< std::string > centralBodyList;
            for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
            {
                centralBodyList.push_back( singleArcTranslationalStateParameter->getCentralBody( ) );
            }
            multiArcParameter = std::make_shared< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< StateScalarType > >(
                    singleArcTranslationalStateParameter->getParameterName( ).second.first,
                    arcStartTimes,
                    multiArcInitialStates,
                    centralBodyList,
                    singleArcTranslationalStateParameter->getFrameOrientation( ) );
            break;
        }
        default:
            throw std::runtime_error( "Error when getting multi-arc parameter from single-arc equivalent, parameter type " +
                                      getParameterTypeString( singleArcParameter->getParameterName( ).first ) + " not recognized." );
    }
    return multiArcParameter;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_GETASSOCIATEDMULTIARCPARAMETER_H
