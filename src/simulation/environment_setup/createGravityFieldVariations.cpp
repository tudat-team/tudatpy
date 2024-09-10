/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/interpolators/interpolator.h"

#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"
#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/gravitation/tabulatedGravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{


//! Function to create a set of gravity field variations, stored in the associated interface class
std::shared_ptr< gravitation::GravityFieldVariationsSet > createGravityFieldModelVariationsSet(
        const std::string& body,
        const SystemOfBodies& bodies,
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


    if( bodies.count( body ) == 0 )
    {
        throw std::runtime_error( "Error when making gravity field variations of body " + body +
                                  ", body not found" );
    }
    // Iterate over all variations to create.
    for( unsigned int i = 0; i < gravityFieldVariationSettings.size( ); i++ )
    {
        // Get current type of deformation
        variationTypes.push_back( gravityFieldVariationSettings.at( i )->getBodyDeformationType( ) );

        // Set current variation object in list.
        variationObjects.push_back( createGravityFieldVariationsModel(
                                        gravityFieldVariationSettings.at( i ), body, bodies  ) );

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
                bodies.at( body )->getGravityFieldModel( ) ) == nullptr )
    {
        throw std::runtime_error( "Error when making gravity field variations of body " + body +
                                  ", base type is not time dependent" );
    }

    std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                bodies.at( body )->getGravityFieldModel( ) )->setFieldVariationSettings( fieldVariationsSet, false );

    return fieldVariationsSet;
}

//! Function to create a single gravity field variation object.
std::shared_ptr< gravitation::GravityFieldVariations > createGravityFieldVariationsModel(
        const std::shared_ptr< GravityFieldVariationSettings > gravityFieldVariationSettings,
        const std::string body,
        const SystemOfBodies& bodies )
{
    using namespace tudat::gravitation;
    
    std::shared_ptr< GravityFieldVariations > gravityFieldVariationModel;
    
    // Check type of variation
    switch( gravityFieldVariationSettings->getBodyDeformationType( ) )
    {
    case basic_solid_body:
    case mode_coupled_solid_body:
    {
        std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > gravityField =
                std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                            bodies.at( body )->getGravityFieldModel( ) );
        if( gravityField == nullptr )
        {
            throw std::runtime_error( "Error when making basic solid-body gravity field variation of body " + body +
                                      ", base type is not time dependent" );
        }

        // Check consistency
        std::vector< std::string > deformingBodies;
        if( gravityFieldVariationSettings->getBodyDeformationType( ) == basic_solid_body )
        {
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
                deformingBodies  = basicSolidBodyGravityVariationSettings->getDeformingBodies( );
            }
        }
        else if( gravityFieldVariationSettings->getBodyDeformationType( ) == mode_coupled_solid_body )
        {
            std::shared_ptr< ModeCoupledSolidBodyGravityFieldVariationSettings >
                modeCoupledSolidBodyGravityVariationSettings =
                std::dynamic_pointer_cast< ModeCoupledSolidBodyGravityFieldVariationSettings >(
                    gravityFieldVariationSettings );
            if( modeCoupledSolidBodyGravityVariationSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected mode-coupled solid body gravity field settings for " + body );
            }
            else
            {
                deformingBodies  = modeCoupledSolidBodyGravityVariationSettings->getDeformingBodies( );
            }
        }
        else
        {
            throw std::runtime_error( "Error, incorrect variation type when making tidal gravity field variation" );
        }

        // Define list of required input.
        std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction;
        std::function< Eigen::Quaterniond( const double ) > deformedBodyOrientationFunction;
        std::vector< std::function< Eigen::Vector6d( const double ) > >
                deformingBodyStateFunctions;
        std::vector< std::function< double( ) > > gravitionalParametersOfDeformingBodies;

        // Iterate over all bodies causing tidal perturbation.
        for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
        {
            // Check if perturbing body exists.
            if( bodies.count( deformingBodies[ i ] ) == 0 )
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
                                bodies.at( deformingBodies[ i ] ), std::placeholders::_1 ) );
            }
            else
            {
                deformingBodyStateFunctions.push_back(
                            std::bind( &Body::getState, bodies.at( deformingBodies[ i ] ) ) );
            }

            // Get gravitational parameter of perturbing bodies.
            if( bodies.at( deformingBodies[ i ] )->getGravityFieldModel( ) == nullptr )
            {
                throw std::runtime_error(
                            "Error, could not find gravity field model in body " + deformingBodies[ i ] +
                            " when making basic sh variation for body " + body );
            }
            else
            {
                gravitionalParametersOfDeformingBodies.push_back(
                            std::bind( &GravityFieldModel::getGravitationalParameter,
                                         bodies.at( deformingBodies[ i ] )
                                         ->getGravityFieldModel( ) ) );
            }
        }

        // Set state and orientation functions of perturbed body.
        if( gravityFieldVariationSettings->getInterpolatorSettings( ) != nullptr )
        {
            deformedBodyStateFunction = std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >,
                                                     bodies.at( body ), std::placeholders::_1 );
            deformedBodyOrientationFunction = std::bind(
                        &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                        bodies.at( body )->getRotationalEphemeris( ), std::placeholders::_1 );
        }
        else
        {
            deformedBodyStateFunction = std::bind( &Body::getState, bodies.at( body ) );
            deformedBodyOrientationFunction = std::bind( &Body::getCurrentRotationToLocalFrame,
                                                           bodies.at( body ) );


        }

        std::function< double( ) > gravitionalParameterOfDeformedBody =
                std::bind( &GravityFieldModel::getGravitationalParameter,
                             bodies.at( body )->getGravityFieldModel( ) );


        // Create basic tidal variation object.
        if( gravityFieldVariationSettings->getBodyDeformationType( ) == basic_solid_body )
        {
            std::shared_ptr< BasicSolidBodyGravityFieldVariationSettings >
                basicSolidBodyGravityVariationSettings =
                std::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >(
                    gravityFieldVariationSettings );
            gravityFieldVariationModel
                = std::make_shared<BasicSolidBodyTideGravityFieldVariations>(
                deformedBodyStateFunction,
                deformedBodyOrientationFunction,
                deformingBodyStateFunctions,
                gravityField->getReferenceRadius( ),
                gravitionalParameterOfDeformedBody,
                gravitionalParametersOfDeformingBodies,
                basicSolidBodyGravityVariationSettings->getLoveNumbers( ),
                deformingBodies );
        }
        else if( gravityFieldVariationSettings->getBodyDeformationType( ) == mode_coupled_solid_body )
        {
            std::shared_ptr< ModeCoupledSolidBodyGravityFieldVariationSettings >
                modeCoupledSolidBodyGravityVariationSettings =
                std::dynamic_pointer_cast< ModeCoupledSolidBodyGravityFieldVariationSettings >(
                    gravityFieldVariationSettings );
            gravityFieldVariationModel
                = std::make_shared<ModeCoupledSolidBodyTideGravityFieldVariations>(
                deformedBodyStateFunction,
                deformedBodyOrientationFunction,
                deformingBodyStateFunctions,
                gravityField->getReferenceRadius( ),
                gravitionalParameterOfDeformedBody,
                gravitionalParametersOfDeformingBodies,
                modeCoupledSolidBodyGravityVariationSettings->getLoveNumbers( ),
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
    case periodic_variation:
    {
        // Check input consistency
        std::shared_ptr< PeriodicGravityFieldVariationsSettings > periodicGravityFieldVariationSettings
                = std::dynamic_pointer_cast< PeriodicGravityFieldVariationsSettings >(
                    gravityFieldVariationSettings );
        if( periodicGravityFieldVariationSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected periodic gravity field variation settings for " + body );
        }
        else
        {
            if( periodicGravityFieldVariationSettings->getCosineShAmplitudesCosineTime( ).size( ) > 0 ||
                periodicGravityFieldVariationSettings->getCosineShAmplitudesSineTime( ).size( ) > 0 ||
                periodicGravityFieldVariationSettings->getSineShAmplitudesCosineTime( ).size( ) > 0 ||
                periodicGravityFieldVariationSettings->getSineShAmplitudesSineTime( ).size( ) > 0)
            {
                // Create variation.
                gravityFieldVariationModel = std::make_shared< PeriodicGravityFieldVariations >
                        (  periodicGravityFieldVariationSettings->getCosineShAmplitudesCosineTime(),
                           periodicGravityFieldVariationSettings->getCosineShAmplitudesSineTime( ),
                           periodicGravityFieldVariationSettings->getSineShAmplitudesCosineTime( ),
                           periodicGravityFieldVariationSettings->getSineShAmplitudesSineTime( ),
                           periodicGravityFieldVariationSettings->getFrequencies( ),
                           periodicGravityFieldVariationSettings->getReferenceEpoch( ),
                           periodicGravityFieldVariationSettings->getMinimumDegree( ),
                           periodicGravityFieldVariationSettings->getMinimumOrder( ) );
            }
            else
            {
                throw std::runtime_error( "Error when creating periodic gravity field variations; no amplitudes found" );
            }
        }

        break;
    }
    case polynomial_variation:
    {
        // Check input consistency
        std::shared_ptr< PolynomialGravityFieldVariationsSettings > polynomialGravityFieldVariationSettings
            = std::dynamic_pointer_cast< PolynomialGravityFieldVariationsSettings >(
                gravityFieldVariationSettings );
        if( polynomialGravityFieldVariationSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected polynomial gravity field variation settings for " + body );
        }
        else
        {
                // Create variation.
                gravityFieldVariationModel = std::make_shared< PolynomialGravityFieldVariations >
                    (  polynomialGravityFieldVariationSettings->getCosineAmplitudes( ),
                       polynomialGravityFieldVariationSettings->getSineAmplitudes( ),
                       polynomialGravityFieldVariationSettings->getReferenceEpoch( ),
                       polynomialGravityFieldVariationSettings->getMinimumDegree( ),
                       polynomialGravityFieldVariationSettings->getMinimumOrder( ) );
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


//! Function to create constant complex Love number list for a range of degrees and orders.
std::map< int, std::vector< std::complex< double > > > getFullLoveNumbersVector(
        const std::complex< double > constantLoveNumber, const int maximumDegree, const int maximumOrder )
{
    // Define list of Love numbers
    std::map< int, std::vector< std::complex< double > > > fullLoveNumbersVector;

    // Iterate over all degrees and orders.
    for( unsigned int i = 2; i <= static_cast< unsigned int >( maximumDegree );i++ )
    {
        for( unsigned int j = 0; j <= ( i + 2 ); j++ )
        {
            // Set current Love numbers.
            if( static_cast< int >( j ) <= maximumOrder  )
            {
                fullLoveNumbersVector[ i ].push_back( constantLoveNumber );
            }
            else
            {
                fullLoveNumbersVector[ i ].push_back( std::complex< double >( 0.0, 0.0 ) );
            }
        }
    }

    return fullLoveNumbersVector;
}

//! Function to create constant real Love number list for a range of degrees and orders.
std::map< int, std::vector< std::complex< double > > > getFullLoveNumbersVector(
        const double constantLoveNumber, const int maximumDegree, const int maximumOrder )
{
    return getFullLoveNumbersVector( std::complex< double >( constantLoveNumber, 0.0 ), maximumDegree, maximumOrder );
}


std::vector< std::complex< double > > getLoveNumberPerDegree(
        const std::complex< double > loveNumber,
        const int degree )
{
    std::vector< std::complex< double > > loveNumbers;
    for( int i = 0; i <= degree; i++ )
    {
        loveNumbers.push_back( loveNumber );
    }
    return loveNumbers;
}

std::vector< std::complex< double > > getLoveNumberPerDegree(
        const double loveNumber,
        const int degree )
{
    return getLoveNumberPerDegree( std::complex< double >( loveNumber, 0.0 ), degree );
}


} // namespace simulation_setup

} // namespace tudat
