#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Function to read BLQ file.
std::vector<std::map<std::string, Eigen::Matrix<double, 3, 11> > > readBlqFile( const std::string &fileName )
{
    std::cout << "Reading BLQ file " << fileName << std::endl;
    using namespace tudat::unit_conversions;

    // Create stream from given file name.
    std::fstream stream( fileName.c_str( ), std::ios::in );

    // Check if stream was successfully created.
    if ( stream.fail( ))
    {
        throw std::runtime_error( "Blq data file " + fileName + " could not be opened." );
    }

    std::vector<std::string> vectorOfIndividualStrings;
    std::vector<std::vector<std::string> > unparsedSiteBlock;
    unparsedSiteBlock.resize( 11 );

//    // Initialized number of parsed sites.
//    int numberOfSitesParsed = 0;

    // Declare current site characteristics
    std::string currentSiteName;
    Eigen::Matrix<double, 6, 11> currentDataMatrix;

    // Declare maps containin return values.
    std::map<std::string, Eigen::Matrix<double, 3, 11> > phaseMap;
    std::map<std::string, Eigen::Matrix<double, 3, 11> > amplitudeMap;

    // Declare variable to determine if file headerhas been passed.
    bool isFileHeaderParsed = 0;

    std::string line;

    // While stream is successful and not at end, continue.
    while ( !stream.fail( ) && !stream.eof( ))
    {
        // Continue this loop until file header has been passed.
        while ( !isFileHeaderParsed )
        {
            // Get line from file.
            std::getline( stream, line );

            // Trim input string (removes all leading and trailing whitespaces).
            boost::algorithm::trim( line );

            // Split string into multiple strings, each containing one element from a line from the
            // data file.
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( " " ),
                                     boost::algorithm::token_compress_on );

            // Check whether initial part of file header is finished.
            if ( vectorOfIndividualStrings[ 0 ] != "$$" )
            {
                isFileHeaderParsed = 1;
            }
        }

        unparsedSiteBlock[ 0 ] = vectorOfIndividualStrings;

        for ( int i = 1; i < 11; i++ )
        {
            std::getline( stream, line );
            boost::algorithm::trim( line );
            boost::algorithm::split( unparsedSiteBlock[ i ],
                                     line,
                                     boost::algorithm::is_any_of( " " ),
                                     boost::algorithm::token_compress_on );
        }

        for ( int i = 0; i < 11; i++ )
        {
            assert( unparsedSiteBlock[ 0 ].size( ) == 1 );
            currentSiteName = unparsedSiteBlock[ 0 ][ 0 ];
            for ( int i = 1; i < 4; i++ )
            {
                assert( unparsedSiteBlock[ i ][ 0 ] == "$$" );
            }
            for ( int i = 4; i < 10; i++ )
            {
                for ( int j = 0; j < 11; j++ )
                {
                    currentDataMatrix( i - 4, j ) = boost::lexical_cast<double>( unparsedSiteBlock[ i ][ j ] );
                }
            }
            assert( unparsedSiteBlock[ 10 ][ 0 ] == "$$" );
            amplitudeMap[ currentSiteName ] = currentDataMatrix.block( 0, 0, 3, 11 );
            phaseMap[ currentSiteName ] =
                convertDegreesToRadians<Eigen::MatrixXd>( currentDataMatrix.block( 3, 0, 3, 11 ));
//            numberOfSitesParsed++;
        }

        if ( !stream.fail( ) && !stream.eof( ))
        {
            std::getline( stream, line );

            // Trim input string (removes all leading and trailing whitespaces).
            boost::algorithm::trim( line );

            // Split string into multiple strings, each containing one element from a line from the
            // data file.
            boost::algorithm::split( vectorOfIndividualStrings,
                                     line,
                                     boost::algorithm::is_any_of( " " ),
                                     boost::algorithm::token_compress_on );
        }
    }

    std::vector<std::map<std::string, Eigen::Matrix<double, 3, 11> > > outputVector;
    outputVector.resize( 2 );
    outputVector[ 0 ] = amplitudeMap;
    outputVector[ 1 ] = phaseMap;

    return outputVector;
}

Eigen::Vector3d convertLocalOceanTideDisplacementToPlentocentricDisplacement(
    const Eigen::Vector3d &localDisplacement,
    const std::vector<Eigen::Vector3d> &earthNorthradialUnitVectors )
{
    Eigen::Vector3d planetocentricDisplacement = Eigen::Vector3d::Zero( );

    for ( int i = 0; i < 3; i++ )
    {
        planetocentricDisplacement += -localDisplacement( 2 - i ) * earthNorthradialUnitVectors[ i ];
    }
    return planetocentricDisplacement;
}

//! Constructor, reads and initializes tide displacement amplitudes and phases.
OceanTideEarthDeformation::OceanTideEarthDeformation(
    const std::vector< std::string >& blqFiles,
    const std::function<Eigen::Vector6d( const double )> doodsonArgumentFunction ) :
    doodsonArgumentFunction_( doodsonArgumentFunction )
{
    for ( unsigned int i = 0; i < blqFiles.size( ); i++ )
    {
        addBlqFile( blqFiles.at( i ) );
    }

    // Initialize characteristics of 11 main tides, from Moyer(2000), Table 5-1 and Eq. 5-88
    doodsonMultipliers_ << 2.0, 0.0, 0.0, 0.0,
        2.0, 2.0, -2.0, 0.0,
        2.0, -1.0, 0.0, 1.0,
        2.0, 2.0, 0.0, 0.0,

        1.0, 1.0, 0.0, 0.0,
        1.0, -1.0, 0.0, 0.0,
        1.0, 1.0, -2.0, 0.0,
        1.0, -2.0, 0.0, 1.0,

        0.0, 2.0, 0.0, 0.0,
        0.0, 1.0, 0.0, -1.0,
        0.0, 0.0, 2.0, 0.0;

    schwiderskiFactors_ <<
                        0.0, 0.0, 0.0, 0.0, 0.5, -.50, -.50, -.50, 0.0, 0.0, 0.0;
    schwiderskiFactors_ *= mathematical_constants::PI;

}

void OceanTideEarthDeformation::addBlqFile( const std::string& blqFile )
{
        // Read BLQ file.
    std::vector<std::map<std::string, Eigen::Matrix<double, 3, 11> > > blqFileData = readBlqFile( blqFile );
    std::map<std::string, Eigen::Matrix<double, 3, 11> > currentAmplitudes = blqFileData.at( 0 );
    std::map<std::string, Eigen::Matrix<double, 3, 11> > currentPhases = blqFileData.at( 1 );

    // Set data retriebed from BLQ files as member variables.
    siteOceanTideAmplitudes_.insert( currentAmplitudes.begin( ), currentAmplitudes.end( ) );
    siteOceanTidePhases_.insert( currentPhases.begin( ), currentPhases.end( ) );
}


Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacementInEnuFrame(
    const double time,
    const std::string& siteIdentifier )
{
    // Check to see if requested site is available.
    if( siteOceanTideAmplitudes_.count( siteIdentifier ) == 0 )
    {
        throw std::runtime_error( "Error, requested site: " + siteIdentifier + ", which is not available for ocean tide displacement calculation (blq file missing)" );
    }

    // Calculate doodson arguments at given time.
    Eigen::Vector6d doodsonArguments = doodsonArgumentFunction_( time );

    // Initialize displacement to zero.
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );

    Eigen::Matrix< double, 11, 1 > arguments;

    // Iterate over all tidal constituents and calculate current phase and resulting displacement.
    for( int i = 0; i < 11; i++ )
    {
        arguments[ i ] = 0.0;

        // Combine influence of doodson arguments (argument 5 and 6 are zero for all 11 constituent tides).
        for( int j = 0; j < 4; j++ )
        {
            arguments[ i ] += doodsonArguments[ j ] * doodsonMultipliers_( i, j );
        }

        // Add schwiderski phase factor for applicable tide.s
        arguments[ i ] += schwiderskiFactors_( i );

        // Calculate site displacement due to tide in three dimensions.
        for( int j = 0; j < 3; j++ )
        {
            displacement( j ) += siteOceanTideAmplitudes_[ siteIdentifier ]( j, i ) *
                                 std::cos ( arguments[ i ] - siteOceanTidePhases_[ siteIdentifier ]( j, i ) );
        }
    }

    return ( Eigen::Vector3d( ) << -displacement( 1 ), -displacement( 2 ), displacement( 0 ) ).finished( );
}


//! Function to calculate site displacement at given time.
Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacement(
    const double time,
    const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    return stationState->getRotationFromBodyFixedToTopocentricFrame( time ).inverse( ) * calculateDisplacementInEnuFrame( time, stationState->getSiteId( ) );
}

}

}
