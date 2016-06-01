/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"

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
    cosineCoefficients.block( minimumDegree_, minimumOrder_, numberOfDegrees_, numberOfOrders_ )
            += correctionPair.first;
}

//! Function to retrieve a variation object of given type (and name if necessary).
std::pair< bool, boost::shared_ptr< gravitation::GravityFieldVariations > >
GravityFieldVariationsSet::getGravityFieldVariation(
        const BodyDeformationTypes deformationType,
        const std::string identifier )
{
    // Declare return pointer.
    boost::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVariation;

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
                    ( variationIdentifier_[ i ] == identifier ) )
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
            std::cerr << "Error when retrieving gravity field variation of type " <<
                       deformationType << ", none of " <<
                       numberOfEntries << " potential entries match identifier." << std::endl;
        }
    }

    return std::make_pair( isEntryFound, gravityFieldVariation );
}

//! Function to create a function linearly interpolating the sine and cosine correction coefficients
//! produced by an object of GravityFieldVariations type.
boost::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) >
createInterpolatedSphericalHarmonicCorrectionFunctions(
        boost::shared_ptr< GravityFieldVariations > variationObject,
        const double initialTime,
        const double finalTime,
        const double timeStep,
        const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
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
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
            cosineSineCorrectionInterpolator =
            interpolators::createOneDimensionalInterpolator< double, Eigen::MatrixXd >(
                cosineSineCorrectionsMap, interpolatorSettings );

    // Create pair interpolation interface for TimeDependentSphericalHarmonicsGravityField
    boost::shared_ptr< PairInterpolationInterface > interpolationInterface =
            boost::make_shared< PairInterpolationInterface >(
                cosineSineCorrectionInterpolator, variationObject->getMinimumDegree( ),
                variationObject->getMinimumOrder( ),
                variationObject->getNumberOfDegrees( ), variationObject->getNumberOfOrders( ) );

    // Create update function.
    return boost::bind( &PairInterpolationInterface::getCosineSinePair,
                        interpolationInterface, _1, _2, _3 );
}

//! Class constructor.
GravityFieldVariationsSet::GravityFieldVariationsSet(
        const std::vector< boost::shared_ptr< GravityFieldVariations > > variationObjects,
        const std::vector< BodyDeformationTypes > variationType,
        const std::vector< std::string > variationIdentifier,
        const std::map< int, boost::shared_ptr< interpolators::InterpolatorSettings > >
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
        std::cerr << "Error when making GravityFieldVariationsSet, inconsistent input, type 1" << std::endl;
    }
    if( variationObjects_.size( ) != variationIdentifier_.size( ) )
    {
        std::cerr << "Error when making GravityFieldVariationsSet, inconsistent input, type 2" << std::endl;
    }

    // Check if interpolation information is provided where required.
    for( std::map< int, boost::shared_ptr< interpolators::InterpolatorSettings > >::iterator
         interpolatorSettingsIterator =
         createInterpolator_.begin( ); interpolatorSettingsIterator != createInterpolator_.end( );
         interpolatorSettingsIterator++ )
    {
        if( initialTimes_.count( interpolatorSettingsIterator->first ) == 0 )
        {
            std::cerr << "Error when making GravityFieldVariationsSet, inconsistent input, type 4, " <<
                       interpolatorSettingsIterator->first << std::endl;
        }

        if( finalTimes_.count( interpolatorSettingsIterator->first ) == 0 )
        {
            std::cerr << "Error when making GravityFieldVariationsSet, inconsistent input, type 5, " <<
                       interpolatorSettingsIterator->first << std::endl;
        }

        if( timeSteps_.count( interpolatorSettingsIterator->first ) == 0 )
        {
            std::cerr << "Error when making GravityFieldVariationsSet, inconsistent input, type 6, " <<
                       interpolatorSettingsIterator->first << std::endl;
        }
    }
}

//! Function to retrieve list of variation functions.
std::vector< boost::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
GravityFieldVariationsSet::getVariationFunctions( )
{
    // Declare list
    std::vector< boost::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
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
                        boost::bind(
                            &GravityFieldVariations::addSphericalHarmonicsCorrections,
                            variationObjects_[ i ], _1, _2, _3 ) );
        }
    }

    // Return list.
    return variationFunctions;
}

} // namespace gravitation

} // namespace tudat
