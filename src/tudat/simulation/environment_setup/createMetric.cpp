/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/relativity/solarSystemMetric.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"
#include <stdexcept>

namespace tudat
{


namespace simulation_setup
{

std::shared_ptr< relativity::Metric > createSpaceTimeMetric(
        const std::shared_ptr< SpaceTimeMetricSettings >& spaceTimeMetricSettings,
        const SystemOfBodies& bodies )
{
    if( spaceTimeMetricSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating space-time metric: input settings are nullptr." );
    }

    std::shared_ptr< relativity::Metric > spaceTimeMetric;
    std::shared_ptr< relativity::PPNParameterSet > ppnSet = bodies.getSpaceTimeProperties( )->getPpnParameterSet( );

    switch( spaceTimeMetricSettings->getMetricType( ) )
    {
    case schwardschild_metric:
    {
        std::shared_ptr< SchwardschildSpaceTimeMetricSettings > schwarzschildSettings =
                std::dynamic_pointer_cast< SchwardschildSpaceTimeMetricSettings >( spaceTimeMetricSettings );

        if (schwarzschildSettings == nullptr)
        {
            throw std::runtime_error( "Error: Expected Schwarzschild metric settings." );
        }

        const std::string& bodyName = schwarzschildSettings->getBodyName( );
        if( !bodies.doesBodyExist( bodyName ) )
        {
            throw std::runtime_error( "Error: Body " + bodyName + " not found in body map." );
        }

        std::shared_ptr< Body > centralBody = bodies.getBody( bodyName );
        if (centralBody->getGravityFieldModel( ) == nullptr)
        {
            throw std::runtime_error( "Error: Body " + bodyName + " has no gravity field model." );
        }

        spaceTimeMetric = std::make_shared< relativity::HarmonicSchwarzschildMetric >(
            std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                       centralBody->getGravityFieldModel( ) ),
            ppnSet,
            bodyName,
            schwarzschildSettings->getIncludeSecondPostNewtonianOrder( ) );

        break;
    }

    case solar_system_metric:
    {
        std::shared_ptr< SolarSystemSpaceTimeMetricSettings > solarSettings =
                std::dynamic_pointer_cast< SolarSystemSpaceTimeMetricSettings >( spaceTimeMetricSettings );

        if (solarSettings == nullptr)
        {
            throw std::runtime_error( "Error: Expected solar system metric settings." );
        }

        std::vector< std::string > bodyList;
        std::vector< std::function< double( ) > > bodyGravitationalParameterFunctions;
        std::vector< std::function< Eigen::Vector6d( ) > > bodyStateFunctions;
        std::vector< std::function< Eigen::Vector3d( const double ) > > bodyAccelerationFunctions;
        std::vector< int > secondOrderBodyList;
        std::map< int, std::function< double( ) > > bodyAngularMomentumFunctions;
        std::map< int, std::shared_ptr< SphericalHarmonicWrapper > > bodySphericalHarmonicGravityWrappers;
        std::map< int, std::function< void( const double ) > > rotationUpdateFunctions;

        const std::vector< std::string >& firstOrderBodies = solarSettings->getBodiesWithFirstOrderExpansion( );
        const std::vector< std::string >& secondOrderBodies = solarSettings->getBodiesWithSecondOrderExpansion( );
        const std::map< std::string, std::pair< int, int > >& harmonicMap = solarSettings->getBodySphericalHarmonicExpansions( );

        for (unsigned int i = 0; i < firstOrderBodies.size( ) + secondOrderBodies.size( ); ++i)
        {
            std::string currentBody;
            if (i < firstOrderBodies.size( ))
            {
                currentBody = firstOrderBodies.at( i );
            }
            else
            {
                currentBody = secondOrderBodies.at( i - firstOrderBodies.size( ) );
                secondOrderBodyList.push_back( i );
            }

            bodyList.push_back( currentBody );

            if (!bodies.doesBodyExist( currentBody ) )
            {
                throw std::runtime_error( "Error: Body " + currentBody + " not found." );
            }

            std::shared_ptr< Body > body = bodies.getBody( currentBody );

            if (body->getGravityFieldModel( ) == nullptr)
            {
                throw std::runtime_error( "Error: Body " + currentBody + " has no gravity field." );
            }

            bodyStateFunctions.push_back( std::bind( &Body::getState, body ) );

            if( body->getEphemeris( ) != nullptr )
            {
                const double accelerationStep = 10.0;
                std::shared_ptr< ephemerides::Ephemeris > ephemeris = body->getEphemeris( );
                bodyAccelerationFunctions.push_back(
                    [ ephemeris, accelerationStep ]( const double t )
                    {
                        return ephemeris->getCartesianAcceleration( t, accelerationStep );
                    } );
            }
            else
            {
                bodyAccelerationFunctions.push_back(
                    []( const double ){ return Eigen::Vector3d::Zero( ); } );
            }

            bodyGravitationalParameterFunctions.push_back(
                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                           body->getGravityFieldModel( ) ) );

            if (harmonicMap.count( currentBody ) > 0)
            {
                std::shared_ptr< gravitation::SphericalHarmonicsGravityField > harmonicField =
                        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                            body->getGravityFieldModel( ) );

                if (harmonicField == nullptr)
                {
                    throw std::runtime_error( "Error: Body " + currentBody + " does not have SH gravity model." );
                }
                if( body->getRotationalEphemeris( ) == nullptr )
                {
                    throw std::runtime_error( "Error: Body " + currentBody +
                                              " has SH metric settings but no rotational ephemeris." );
                }

                std::function< Eigen::Matrix3d( ) > rotationDerivativeFunction =
                        std::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, body );
                std::function< Eigen::Quaterniond( ) > rotationFunction =
                        std::bind( &Body::getCurrentRotationToLocalFrame, body );

                rotationUpdateFunctions[ i ] =
                        [ body, currentBody ]( const double t )
                        {
                            body->setCurrentRotationalStateToLocalFrameFromEphemeris( t );
                        };

                bodySphericalHarmonicGravityWrappers[ i ] = std::make_shared< SphericalHarmonicWrapper >(
                    harmonicField,
                    std::bind( &Body::getPosition, body ),
                    harmonicMap.at( currentBody ),
                    std::make_shared< basic_mathematics::SphericalHarmonicsCache >(
                        harmonicMap.at( currentBody ).first + 1,
                        harmonicMap.at( currentBody ).second + 1 ),
                    rotationFunction,
                    rotationDerivativeFunction );
            }
        }

        if( !solarSettings->getUseBodyAccelerations( ) )
        {
            bodyAccelerationFunctions.clear( );
            for( unsigned int i = 0; i < bodyList.size( ); ++i )
            {
                bodyAccelerationFunctions.push_back(
                    []( const double ){ return Eigen::Vector3d::Zero( ); } );
            }
        }

        spaceTimeMetric = std::make_shared< relativity::SolarSystemMetric >(
            bodyList,
            bodyGravitationalParameterFunctions,
            bodyStateFunctions,
            secondOrderBodyList,
            ppnSet,
            bodyAccelerationFunctions,
            bodyAngularMomentumFunctions,
            bodySphericalHarmonicGravityWrappers,
            rotationUpdateFunctions );

        break;
    }
    default:
        throw std::runtime_error(
            "Error when creating space-time metric: metric type " +
            std::to_string( static_cast< int >( spaceTimeMetricSettings->getMetricType( ) ) ) +
            " is not recognized." );
    }

    return spaceTimeMetric;
}

} // namespace simulation_setup

} // namespace tudat
