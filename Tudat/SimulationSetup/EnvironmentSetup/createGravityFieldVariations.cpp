/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/tabulatedGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{


//! Function to create a set of gravity field variations, stored in the associated interface class
std::shared_ptr< gravitation::GravityFieldVariationsSet > createGravityFieldModelVariationsSet(
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< std::shared_ptr< GravityFieldVariationSettings > >&
            gravityFieldVariationSettings )
{

    using namespace tudat::gravitation;

    // Declare lists for input to GravityFieldVariationsSet
    std::vector< std::shared_ptr< GravityFieldVariations > > variationObjects;
    std::vector< BodyDeformationTypes > variationTypes;
    std::vector< std::string > variationIdentifiers;
    std::map< int, std::shared_ptr< interpolators::InterpolatorSettings > > createInterpolators;
    std::map< int, double > initialTimes;
    std::map< int, double > finalTimes;
    std::map< int, double > timeSteps;

    // Iterate over all variations to create.
    for( unsigned int i = 0; i < gravityFieldVariationSettings.size( ); i++ )
    {
        // Get current type of deformation
        variationTypes.push_back( gravityFieldVariationSettings.at( i )->getBodyDeformationType( ) );

        // Set current variation object in list.
        variationObjects.push_back( createGravityFieldVariationsModel(
                                        gravityFieldVariationSettings.at( i ), body, bodyMap  ) );

        if( gravityFieldVariationSettings.at( i )->getBodyDeformationType( ) == basic_solid_body )
        {
           std::vector< std::string > deformingBodies =
                   std::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >(
                       gravityFieldVariationSettings.at( i ) )->getDeformingBodies( );
           std::string currentIdentifier = deformingBodies.at( 0 );
           for( unsigned int j = 1; j < deformingBodies.size( ); j++ )
           {
               currentIdentifier += "_" + deformingBodies.at( j );
           }
           variationIdentifiers.push_back( currentIdentifier );
        }
        else
        {
            variationIdentifiers.push_back( "" );
        }

        // Check if current variation is interpolated, and set settings if necessary.
        if( gravityFieldVariationSettings.at( i )->getInterpolatorSettings( ) != nullptr )
        {
            createInterpolators[ i ]
                    = gravityFieldVariationSettings.at( i )->getInterpolatorSettings( )
                          ->interpolatorSettings_;
            initialTimes[ i ]
                    = gravityFieldVariationSettings.at( i )->getInterpolatorSettings( )
                          ->initialTime_;
            finalTimes[ i ]
                    = gravityFieldVariationSettings.at( i )->getInterpolatorSettings( )->finalTime_;
            timeSteps[ i ]
                    = gravityFieldVariationSettings.at( i )->getInterpolatorSettings( )->timeStep_;
        }
    }

    // Create object with settings for updating variations from new parameter values.
    std::shared_ptr< GravityFieldVariationsSet > fieldVariationsSet =
            std::make_shared< GravityFieldVariationsSet >(
                variationObjects, variationTypes, variationIdentifiers,
                createInterpolators, initialTimes, finalTimes, timeSteps );

    if( std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                bodyMap.at( body )->getGravityFieldModel( ) ) == nullptr )
    {
        throw std::runtime_error( "Error when making gravity field variations of body " + body +
                                  ", base type is not time dependent" );
    }

    std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                bodyMap.at( body )->getGravityFieldModel( ) )->setFieldVariationSettings( fieldVariationsSet, false );

    return fieldVariationsSet;
}

//! Function to create a single gravity field variation object.
std::shared_ptr< gravitation::GravityFieldVariations > createGravityFieldVariationsModel(
        const std::shared_ptr< GravityFieldVariationSettings > gravityFieldVariationSettings,
        const std::string body,
        const NamedBodyMap& bodyMap )
{
    using namespace tudat::gravitation;
    
    std::shared_ptr< GravityFieldVariations > gravityFieldVariationModel;
    
    // Check type of variation
    switch( gravityFieldVariationSettings->getBodyDeformationType( ) )
    {
    case basic_solid_body:
    {
        // Check consistency
        std::shared_ptr< BasicSolidBodyGravityFieldVariationSettings >
                basicSolidBodyGravityVariationSettings =
                std::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >(
                    gravityFieldVariationSettings );
        if( basicSolidBodyGravityVariationSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected basic solid body gravity field settings for " + body );
        }
        else
        {
            // Define list of required input.
            std::vector< std::string > deformingBodies
                    = basicSolidBodyGravityVariationSettings->getDeformingBodies( );
            std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction;
            std::function< Eigen::Quaterniond( const double ) > deformedBodyOrientationFunction;
            std::vector< std::function< Eigen::Vector6d( const double ) > >
                    deformingBodyStateFunctions;
            std::vector< std::function< double( ) > > gravitionalParametersOfDeformingBodies;

            // Iterate over all bodies causing tidal perturbation.
            for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
            {
                // Check if perturbing body exists.
                if( bodyMap.count( deformingBodies[ i ] ) == 0 )
                {
                    throw std::runtime_error( "Error when making basic solid body gravity field variation, " +
                                               deformingBodies[ i ] + " deforming body not found." );
                }
                
                // Create body state functions (depending on whether the variation is calculated
                // directly during propagation, or a priori by an interpolator
                if( gravityFieldVariationSettings->getInterpolatorSettings( ) != nullptr )
                {
                    deformingBodyStateFunctions.push_back(
                                std::bind(
                                    &Body::getStateInBaseFrameFromEphemeris< double, double >,
                                    bodyMap.at( deformingBodies[ i ] ), std::placeholders::_1 ) );
                }
                else
                {
                    deformingBodyStateFunctions.push_back(
                                std::bind( &Body::getState, bodyMap.at( deformingBodies[ i ] ) ) );
                }

                // Get gravitational parameter of perturbing bodies.
                if( bodyMap.at( deformingBodies[ i ] )->getGravityFieldModel( ) == nullptr )
                {
                    throw std::runtime_error(
                                "Error, could not find gravity field model in body " + deformingBodies[ i ] +
                                " when making basic sh variation for body " + body );
                }
                else
                {
                    gravitionalParametersOfDeformingBodies.push_back(
                                std::bind( &GravityFieldModel::getGravitationalParameter,
                                             bodyMap.at( deformingBodies[ i ] )
                                             ->getGravityFieldModel( ) ) );
                }
            }

            // Set state and orientation functions of perturbed body.
            if( gravityFieldVariationSettings->getInterpolatorSettings( ) != nullptr )
            {
                deformedBodyStateFunction = std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >,
                                                         bodyMap.at( body ), std::placeholders::_1 );
                deformedBodyOrientationFunction = std::bind(
                            &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                            bodyMap.at( body )->getRotationalEphemeris( ), std::placeholders::_1 );
            }
            else
            {
                deformedBodyStateFunction = std::bind( &Body::getState, bodyMap.at( body ) );
                deformedBodyOrientationFunction = std::bind( &Body::getCurrentRotationToLocalFrame,
                                                               bodyMap.at( body ) );

                
            }

            std::function< double( ) > gravitionalParameterOfDeformedBody =
                    std::bind( &GravityFieldModel::getGravitationalParameter,
                                 bodyMap.at( body )->getGravityFieldModel( ) );
            
            // Create basic tidal variation object.
            gravityFieldVariationModel
                    = std::make_shared< BasicSolidBodyTideGravityFieldVariations >(
                        deformedBodyStateFunction,
                        deformedBodyOrientationFunction,
                        deformingBodyStateFunctions,
                        basicSolidBodyGravityVariationSettings->getBodyReferenceRadius( ),
                        gravitionalParameterOfDeformedBody,
                        gravitionalParametersOfDeformingBodies,
                        basicSolidBodyGravityVariationSettings->getLoveNumbers( ),
                        deformingBodies );
        }
        break;
    }    
    case tabulated_variation:
    {
        // Check input consistency
        std::shared_ptr< TabulatedGravityFieldVariationSettings > tabulatedGravityFieldVariationSettings
                = std::dynamic_pointer_cast< TabulatedGravityFieldVariationSettings >(
                    gravityFieldVariationSettings );
        if( tabulatedGravityFieldVariationSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected tabulated gravity field variation settings for " + body );
        }
        else
        {
            // Create variation.
            gravityFieldVariationModel = std::make_shared< TabulatedGravityFieldVariations >
                    (  tabulatedGravityFieldVariationSettings->getCosineCoefficientCorrections( ),
                       tabulatedGravityFieldVariationSettings->getSineCoefficientCorrections( ),
                       tabulatedGravityFieldVariationSettings->getMinimumDegree( ),
                       tabulatedGravityFieldVariationSettings->getMinimumOrder( ),
                       tabulatedGravityFieldVariationSettings->getInterpolatorSettings( )->interpolatorSettings_ );
        }
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, case " + std::to_string(
                                       gravityFieldVariationSettings->getBodyDeformationType( ) ) +
                                   " not implemented for gravity field variations." );
    }
        
    }
    
    if( gravityFieldVariationModel == nullptr )
    {
        throw std::runtime_error( "Gravity variation model IS nullptr after creation." );
    }
    
    return gravityFieldVariationModel;
    
}

} // namespace simulation_setup

} // namespace tudat
