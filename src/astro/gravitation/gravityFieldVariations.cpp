/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/basics/utilities.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/math/interpolators/linearInterpolator.h"

namespace tudat
{

namespace gravitation
{

//! Function to add sine and cosine corrections at given time to coefficient matrices.
void PairInterpolationInterface::getCosineSinePair(
        const double time, Eigen::MatrixXd& sineCoefficients, Eigen::MatrixXd& cosineCoefficients )
{
    // Interpolate corrections
    Eigen::MatrixXd cosineSinePair = cosineSineInterpolator_->interpolate( time );

    // Split combined interpolated cosine/sine correction block and add to existing values
    cosineCoefficients.block( startDegree_, startOrder_, numberOfDegrees_, numberOfOrders_ ) +=
            cosineSinePair.block( 0, 0, cosineSinePair.rows( ), cosineSinePair.cols( ) / 2 );
    sineCoefficients.block( startDegree_, startOrder_, numberOfDegrees_, numberOfOrders_ ) +=
            cosineSinePair.block( 0, cosineSinePair.cols( ) / 2,
                                  cosineSinePair.rows( ), cosineSinePair.cols( ) / 2 );
}


//! Function to add sine and cosine corrections at given time to coefficient matrices.
void GravityFieldVariations::addSphericalHarmonicsCorrections(
        const double time, Eigen::MatrixXd& sineCoefficients, Eigen::MatrixXd& cosineCoefficients )
{
    // Calculate corrections.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > correctionPair =
            calculateSphericalHarmonicsCorrections( time );

    // Add corrections to existing values
    sineCoefficients.block( minimumDegree_, minimumOrder_, numberOfDegrees_, numberOfOrders_ )
            += correctionPair.second;
    lastSineCorrection_.block( minimumDegree_, minimumOrder_, numberOfDegrees_, numberOfOrders_ )
            = correctionPair.second;

    cosineCoefficients.block( minimumDegree_, minimumOrder_, numberOfDegrees_, numberOfOrders_ )
            += correctionPair.first;
    lastCosineCorrection_.block( minimumDegree_, minimumOrder_, numberOfDegrees_, numberOfOrders_ )
            = correctionPair.first;
}

//! Function to retrieve a variation object of given type (and name if necessary).
std::pair< bool, std::shared_ptr< gravitation::GravityFieldVariations > >
GravityFieldVariationsSet::getGravityFieldVariation(
        const BodyDeformationTypes deformationType,
        const std::string identifier )
{
    // Declare return pointer.
    std::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVariation;

    // Check how many variation objects of request type are in list.
    int numberOfEntries =
            std::count( variationType_.begin( ), variationType_.end( ), deformationType );

    // Check if number of variation objects of requested type is not zero.
    bool isEntryFound = 1;
    if( numberOfEntries == 0 )
    {
        isEntryFound = 0;
    }
    // If number of objects of requested type is 1, no further search is required
    else if( numberOfEntries == 1 )
    {
        // Retrieve index of object and set as return pointer.
        std::vector< BodyDeformationTypes >::iterator findIterator =
                std::find( variationType_.begin( ), variationType_.end( ), deformationType );
        int vectorIndex = std::distance( variationType_.begin( ), findIterator );
        gravityFieldVariation = variationObjects_[ vectorIndex ];
    }
    // If number of objects of requested type is > 1, search for object with required identifier
    else
    {
        bool isCorrectIdentifierFound = 0;

        // Loop over all entries
        for( unsigned int i = 0; i < variationType_.size( ); i++ )
        {
            // Check if type and identifer match
            if( ( variationType_[ i ] == deformationType ) &&
                    ( ( variationIdentifier_[ i ] == identifier ) || ( identifierPerType_.at( variationType_[ i ] ).size( ) == 1 ) ) )
            {
                // Set return pointer and exit loop.
                gravityFieldVariation = variationObjects_[ i ];
                isCorrectIdentifierFound = 1;
                break;
            }
        }

        // Provide warning if no matches are found
        if( isCorrectIdentifierFound == 0 )
        {
            std::string errorMessage = "Error when retrieving gravity field variation of type " +
                    std::to_string( deformationType ) +
                    ", none of " +
                    std::to_string( numberOfEntries ) +
                    " potential entries match identifier.";
            throw std::runtime_error( errorMessage );
        }
    }

    return std::make_pair( isEntryFound, gravityFieldVariation );
}

//! Function to create a function linearly interpolating the sine and cosine correction coefficients
//! produced by an object of GravityFieldVariations type.
std::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) >
createInterpolatedSphericalHarmonicCorrectionFunctions(
        std::shared_ptr< GravityFieldVariations > variationObject,
        const double initialTime,
        const double finalTime,
        const double timeStep,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    // Declare map of combined cosine and since corrections, to be filled and passed to interpolator
    std::map< double, Eigen::MatrixXd > cosineSineCorrectionsMap;

    // Perform once for matrices size determination.
    double currentTime = initialTime;
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > singleCorrections;
    singleCorrections = variationObject->calculateSphericalHarmonicsCorrections( currentTime );
    int correctionDegrees = singleCorrections.first.rows( );
    int correctionOrders = singleCorrections.first.cols( );

    // Loop over all times at which corrections are to be calculated.
    Eigen::MatrixXd cosineSineCorrections;
    while( currentTime < finalTime )
    {
        // Calculate current corrections.
        singleCorrections = variationObject->calculateSphericalHarmonicsCorrections( currentTime );

        // Set current corrections in single block.
        cosineSineCorrections = Eigen::MatrixXd::Zero( correctionDegrees, 2 * correctionOrders );
        cosineSineCorrections.block( 0, 0, correctionDegrees, correctionOrders ) +=
                singleCorrections.first;
        cosineSineCorrections.block( 0, correctionOrders, correctionDegrees, correctionOrders ) +=
                singleCorrections.second;
        cosineSineCorrectionsMap[ currentTime ] = cosineSineCorrections;

        // Increment time.
        currentTime += timeStep;
    }


    // Create interpolator
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            cosineSineCorrectionInterpolator =
            interpolators::createOneDimensionalInterpolator< double, Eigen::MatrixXd >(
                cosineSineCorrectionsMap, interpolatorSettings );

    // Create pair interpolation interface for TimeDependentSphericalHarmonicsGravityField
    std::shared_ptr< PairInterpolationInterface > interpolationInterface =
            std::make_shared< PairInterpolationInterface >(
                cosineSineCorrectionInterpolator, variationObject->getMinimumDegree( ),
                variationObject->getMinimumOrder( ),
                variationObject->getNumberOfDegrees( ), variationObject->getNumberOfOrders( ) );

    // Create update function.
    return std::bind( &PairInterpolationInterface::getCosineSinePair,
                        interpolationInterface, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );
}

//! Class constructor.
GravityFieldVariationsSet::GravityFieldVariationsSet(
        const std::vector< std::shared_ptr< GravityFieldVariations > > variationObjects,
        const std::vector< BodyDeformationTypes > variationType,
        const std::vector< std::string > variationIdentifier,
        const std::map< int, std::shared_ptr< interpolators::InterpolatorSettings > >
        createInterpolator,
        const std::map< int, double > initialTimes,
        const std::map< int, double > finalTimes,
        const std::map< int, double > timeSteps ):
    variationObjects_( variationObjects ), variationType_( variationType ),
    variationIdentifier_( variationIdentifier ),
    createInterpolator_( createInterpolator ),
    initialTimes_( initialTimes ), finalTimes_( finalTimes ), timeSteps_( timeSteps )
{
    // Check consistency of input data vector sizes.
    if( variationObjects_.size( ) != variationType_.size( ) )
    {
        throw std::runtime_error( "Error when making GravityFieldVariationsSet, inconsistent input, type 1" );
    }
    if( variationObjects_.size( ) != variationIdentifier_.size( ) )
    {
        throw std::runtime_error(  "Error when making GravityFieldVariationsSet, inconsistent input, type 2" );
    }

    for( unsigned int i = 0; i < variationType_.size( ); i++ )
    {
        identifierPerType_[ variationType_.at( i ) ].push_back( variationIdentifier_.at( i ) );
    }

    // Check if interpolation information is provided where required.
    for( std::map< int, std::shared_ptr< interpolators::InterpolatorSettings > >::iterator
         interpolatorSettingsIterator =
         createInterpolator_.begin( ); interpolatorSettingsIterator != createInterpolator_.end( );
         interpolatorSettingsIterator++ )
    {
        if( initialTimes_.count( interpolatorSettingsIterator->first ) == 0 )
        {
            std::string errorMessage = "Error when making GravityFieldVariationsSet, inconsistent input, type 4, " +
                    std::to_string( interpolatorSettingsIterator->first );
            throw std::runtime_error( errorMessage );
        }

        if( finalTimes_.count( interpolatorSettingsIterator->first ) == 0 )
        {
            std::string errorMessage = "Error when making GravityFieldVariationsSet, inconsistent input, type 5, " +
                    std::to_string( interpolatorSettingsIterator->first );
            throw std::runtime_error( errorMessage );
        }

        if( timeSteps_.count( interpolatorSettingsIterator->first ) == 0 )
        {
            std::string errorMessage = "Error when making GravityFieldVariationsSet, inconsistent input, type 6, " +
                    std::to_string( interpolatorSettingsIterator->first );
            throw std::runtime_error( errorMessage );
        }
    }
}

//! Function to retrieve list of variation functions.
std::vector< std::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
GravityFieldVariationsSet::getVariationFunctions( )
{
    // Declare list
    std::vector< std::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
            variationFunctions;

    // Iterate over all corrections and add correction function.
    for( unsigned int i = 0; i < variationObjects_.size( ); i++ )
    {
        // If current variation is to be interpolated, create interpolation function and add to list
        if( createInterpolator_.count( i ) > 0 )
        {
            variationFunctions.push_back(
                        createInterpolatedSphericalHarmonicCorrectionFunctions(
                            variationObjects_[ i ], initialTimes_[ i ], finalTimes_[ i ],
                            timeSteps_[ i ], createInterpolator_[ i ] ) );
        }
        // If current variation is to not be interpolated, create function directly by function
        // pointer by to current GravityFieldVariations object.
        else
        {
            variationFunctions.push_back(
                        std::bind(
                            &GravityFieldVariations::addSphericalHarmonicsCorrections,
                            variationObjects_[ i ], std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 ) );
        }
    }

    // Return list.
    return variationFunctions;
}

//! Function to retrieve the tidal gravity field variation with the specified bodies causing deformation
std::shared_ptr< GravityFieldVariations > GravityFieldVariationsSet::getDirectTidalGravityFieldVariation(
        const std::vector< std::string >& namesOfBodiesCausingDeformation,
        const BodyDeformationTypes tideType )
{
    std::shared_ptr< GravityFieldVariations > gravityFieldVariation;

    // Check number of variation objects of correct type
    int numberOfBasicModels = std::count( variationType_.begin( ), variationType_.end( ), tideType );

    // If one model of correct type is found, check if it is consistent with input
    if( ( numberOfBasicModels ) == 1 )
    {
        // Retrieve variation object. It is the return object if no namesOfBodiesCausingDeformation are given
        gravityFieldVariation = variationObjects_.at(
                    ( std::distance( variationType_.begin( ), std::find( variationType_.begin( ), variationType_.end( ),
                                                                         tideType ) ) ) );
        std::shared_ptr< SolidBodyTideGravityFieldVariations > tidalGravityFieldVariation =
                std::dynamic_pointer_cast< SolidBodyTideGravityFieldVariations >( gravityFieldVariation );

        // Check consistency
        if( tidalGravityFieldVariation == nullptr )
        {
            throw std::runtime_error( "Error when getting direct tidal gravity field variation, one model found, but type does not match" );
        }
        // Check if consistency with input needs to be determined
        else if( namesOfBodiesCausingDeformation.size( ) != 0 )
        {
            bool doBodiesMatch = utilities::doStlVectorContentsMatch(
                        tidalGravityFieldVariation->getDeformingBodies( ), namesOfBodiesCausingDeformation );

            if( !doBodiesMatch )
            {
                throw std::runtime_error(
                            "Error when getting direct tidal gravity field variation, one model found, but deforming bodies do not match" );
            }
        }
    }
    // If no models found, throw exception
    else if( numberOfBasicModels == 0 )
    {
        throw std::runtime_error( "Error when getting direct tidal gravity field variation, no such model found" );
    }
    else
    {
        // If no input bodies are given, no object can be returned: no information is available to choose from list
        if( namesOfBodiesCausingDeformation.size( ) == 0 )
        {
            throw std::runtime_error(
                        "Error when getting direct tidal gravity field variation, found multiple models, but not id is provided" );
        }
        bool isVariationFound = 0;

        // Iterate over all objects and check consistency with input
        for( unsigned int i = 0; i < variationType_.size( ); i++ )
        {
            if( std::dynamic_pointer_cast< gravitation::SolidBodyTideGravityFieldVariations >(
                        variationObjects_.at( i ) ) != nullptr )
            {
                bool doBodiesMatch = utilities::doStlVectorContentsMatch(
                            std::dynamic_pointer_cast< gravitation::SolidBodyTideGravityFieldVariations >(
                                variationObjects_.at( i ) )->getDeformingBodies( ), namesOfBodiesCausingDeformation );

                if( doBodiesMatch )
                {
                    gravityFieldVariation = variationObjects_.at( i );
                    isVariationFound = 1;
                }
            }
        }

        if( isVariationFound == 0 )
        {
            throw std::runtime_error(
                        "Error when getting direct tidal gravity field variation, found multiple models of correct type, but none coincide with ids in list" );
        }
    }

    return gravityFieldVariation;
}

//! Function to retrieve the tidal gravity field variations
std::vector< std::shared_ptr< GravityFieldVariations > > GravityFieldVariationsSet::getDirectTidalGravityFieldVariations( )
{
    std::vector< std::shared_ptr< GravityFieldVariations > > directTidalVariations;

    // Iterate over variation objects and check type
    for( unsigned int i = 0; i < variationType_.size( ); i++ )
    {
        if( variationType_[ i ] == basic_solid_body || variationType_[ i ] == mode_coupled_solid_body )
        {
            directTidalVariations.push_back( variationObjects_[ i ] );
        }
    }

    return directTidalVariations;
}

} // namespace gravitation

} // namespace tudat
