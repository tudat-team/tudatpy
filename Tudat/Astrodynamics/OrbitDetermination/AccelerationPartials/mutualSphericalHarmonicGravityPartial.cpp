/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/mutualSphericalHarmonicGravityPartial.h"

namespace tudat
{

namespace acceleration_partials
{

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > MutualSphericalHarmonicsGravityPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > parameterPartial;

    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        if( parameter->getParameterName( ).second.first == acceleratingBody_ ||
                ( parameter->getParameterName( ).second.first == acceleratedBody_ && accelerationUsesMutualAttraction_ ) )
        {
            parameterPartial = std::make_pair(
                        boost::bind( &MutualSphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody, this, _1 ), 1 );
        }
        else
        {
            parameterPartial = std::make_pair( boost::function< void( Eigen::MatrixXd& ) >( ), 0 );
        }
    }
    else
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyExertingAcceleration =
                accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyUndergoingAcceleration =
                accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );

        parameterPartial = orbit_determination::createMergedParameterPartialFunction(
                    partialFunctionFromBodyExertingAcceleration, partialFunctionFromBodyUndergoingAcceleration );
    }
    return parameterPartial;
}

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > MutualSphericalHarmonicsGravityPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyExertingAcceleration =
            accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromBodyUndergoingAcceleration =
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );

    return  orbit_determination::createMergedParameterPartialFunction(
                partialFunctionFromBodyExertingAcceleration, partialFunctionFromBodyUndergoingAcceleration );
}


void MutualSphericalHarmonicsGravityPartial::update( const double currentTime )
{
    accelerationPartialOfShExpansionOfBodyExertingAcceleration_->update( currentTime );
    accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->update( currentTime );

    currentTime_ = currentTime;

}

}

}
