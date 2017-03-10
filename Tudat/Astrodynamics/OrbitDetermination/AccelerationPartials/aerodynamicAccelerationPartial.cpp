/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
CentralGravitationPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter )
    {
        // If parameter is gravitational parameter, check and create dependency function .
        partialFunctionPair = getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
    }
    else
    {
        partialFunctionPair = std::make_pair( boost::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}

//! Function to create a function returning the current partial w.r.t. a gravitational parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
CentralGravitationPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    // Check if parameter is gravitational parameter.
    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check if parameter body is central body.
        if( parameterId.second.first == acceleratingBody_ )
        {
            partialFunction = boost::bind( &CentralGravitationPartial::wrtGravitationalParameterOfCentralBody,
                                           this, _1 );
            numberOfColumns = 1;

        }

        // Check if parameter body is accelerated body, and if the mutual acceleration is used.
        if( parameterId.second.first == acceleratedBody_ )
        {
            if( accelerationUsesMutualAttraction_ )
            {
                partialFunction = boost::bind( &CentralGravitationPartial::wrtGravitationalParameterOfCentralBody,
                                               this, _1 );
                numberOfColumns = 1;
            }
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void CentralGravitationPartial::wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial = computePartialOfCentralGravityWrtGravitationalParameter(
                currentAcceleratedBodyState_, currentCentralBodyState_ );
}



}

}
