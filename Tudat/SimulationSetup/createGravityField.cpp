/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#endif

#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"

#include "Tudat/SimulationSetup/createGravityField.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to read a gravity field file
std::pair< double, double  > readGravityFieldFile(
        const std::string& fileName, const int maximumDegree, const int maximumOrder,
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd >& coefficients,
        const int gravitationalParameterIndex, const int referenceRadiusIndex )
{
    // Attempt to open gravity file.
    std::fstream stream( fileName.c_str( ), std::ios::in );
    if( stream.fail( ) )
    {
        boost::throw_exception(
                    std::runtime_error( "Pds gravity field data file could not be opened." ) );
    }

    // Declare variables for reading file.
    std::vector< std::string > vectorOfIndividualStrings;
    vectorOfIndividualStrings.resize( 4 );
    std::string line;


    double gravitationalParameter = TUDAT_NAN;
    double referenceRadius = TUDAT_NAN;

    if( ( gravitationalParameterIndex >= 0 ) &&
            ( referenceRadiusIndex >= 0 ) )
    {
        // Get first line of file.
        std::getline( stream, line );

        // Get reference radius and gravitational parameter from first line of file.
        boost::algorithm::trim( line );
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( "\t, " ),
                                 boost::algorithm::token_compress_on );
        if( gravitationalParameterIndex >= static_cast< int >( vectorOfIndividualStrings.size( ) ) ||
                referenceRadiusIndex >= static_cast< int >( vectorOfIndividualStrings.size( ) ) )
        {
            throw std::runtime_error( "Error when reading gravity field file, requested header index exceeds file contents" );
        }

        gravitationalParameter = boost::lexical_cast< double >( vectorOfIndividualStrings[ gravitationalParameterIndex ] );
        referenceRadius = boost::lexical_cast< double >( vectorOfIndividualStrings[ referenceRadiusIndex ] );
    }
    else if( ( !( gravitationalParameterIndex >= 0 ) &&
              ( referenceRadiusIndex >= 0 ) ) ||
             ( ( gravitationalParameterIndex >= 0 ) &&
                           !( referenceRadiusIndex >= 0 ) ) )
    {
        throw std::runtime_error( "Error when reading gravity field file, must retrieve either both or neither of Re and mu" );
    }


    // Declare variables for reading in cosine and sine coefficients.
    int currentDegree = 0, currentOrder = 0;
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd( maximumDegree + 1, maximumOrder + 1 );
    cosineCoefficients.setZero( );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd( maximumDegree + 1, maximumOrder + 1 );
    sineCoefficients.setZero( );

    // Read coefficients up to required maximum degree and order.
    while ( !stream.fail( ) && !stream.eof( ) &&
            ( currentDegree <= maximumDegree || currentOrder <= maximumOrder )  )
    {
        // Read current line
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Split string into multiple strings, each containing one element from a line from the
        // data file.
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( ", " ),
                                 boost::algorithm::token_compress_on );

        // Check current line for consistency
        if( vectorOfIndividualStrings.size( ) < 4 )
        {
            std::cerr<<"Error when reading pds gravity field file, number of fields is "
                    <<vectorOfIndividualStrings.size( )<<std::endl;
        }
        else
        {
            // Read current degree and orde from line.
            currentDegree = boost::lexical_cast< int >( vectorOfIndividualStrings[ 0 ] );
            currentOrder = boost::lexical_cast< int >( vectorOfIndividualStrings[ 1 ] );

            // Set cosine and sine coefficients for current degree and order.
            if( currentDegree <= maximumDegree && currentOrder <= maximumOrder )
            {
                cosineCoefficients( currentDegree, currentOrder ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings[ 2 ] );
                sineCoefficients( currentDegree, currentOrder ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] );
            }
        }
    }

    // Set cosine coefficient at (0,0) to 1.
    cosineCoefficients( 0, 0 ) = 1.0;
    coefficients = std::make_pair( cosineCoefficients, sineCoefficients );

    return std::make_pair( gravitationalParameter, referenceRadius );
}

//! Function to create a gravity field model.
boost::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const boost::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings )
{
    using namespace tudat::gravitation;

    // Declare return object.
    boost::shared_ptr< GravityFieldModel > gravityFieldModel;

    // Check which type of gravity field model is to be created.
    switch( gravityFieldSettings->getGravityFieldType( ) )
    {
    case central:
    {
        // Check whether settings for point mass gravity field model are consistent with its type.
        boost::shared_ptr< CentralGravityFieldSettings > centralFieldSettings =
                boost::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        if( centralFieldSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected central field settings when making gravity field model for body " +
                        body);
        }
        else if( gravityFieldVariationSettings.size( ) != 0 )
        {
            std::cerr<<"Error, requested central gravity field, but field variations settings are not empty."<<std::endl;
        }
        else
        {
            // Create and initialize point mass gravity field model.
            gravityFieldModel = boost::make_shared< GravityFieldModel >(
                        centralFieldSettings->getGravitationalParameter( ) );
        }
        break;
    }
    #if USE_CSPICE
    case central_spice:
    {
        if( gravityFieldVariationSettings.size( ) != 0 )
        {
            std::cerr<<"Error, requested central gravity field, but field variations settings are not empty."<<std::endl;
        }
        else
        {
            // Create and initialize point mass gravity field model from Spice.
            gravityFieldModel = boost::make_shared< GravityFieldModel >(
                        spice_interface::getBodyGravitationalParameter( body ) );
        }

        break;
    }
    #endif
    case spherical_harmonic:
    {
        // Check whether settings for spherical harmonic gravity field model are consistent with
        // its type.
        boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicFieldSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                    gravityFieldSettings );

        if( sphericalHarmonicFieldSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected spherical harmonic field settings when making gravity field model of "
                        + body );
        }
        else
        {
            // Check consistency of cosine and sine coefficients.
            if( ( sphericalHarmonicFieldSettings->getCosineCoefficients( ).rows( ) !=
                  sphericalHarmonicFieldSettings->getSineCoefficients( ).rows( ) ) ||
                    ( sphericalHarmonicFieldSettings->getCosineCoefficients( ).cols( ) !=
                      sphericalHarmonicFieldSettings->getSineCoefficients( ).cols( ) ) )
            {
                throw std::runtime_error(
                            std::string( "Error when making spherical harmonic field, sine and " ) +
                            std::string( "cosine matrix  sizes are not equal for body " ) + body );
            }
            else
            {

                if( gravityFieldVariationSettings.size( ) == 0 &&
                        sphericalHarmonicFieldSettings->getCreateTimeDependentField( ) == 0 )
                {
                    // Create and initialize spherical harmonic gravity field model.
                    gravityFieldModel = boost::make_shared< SphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( ) );
                }
                else
                {
                    if( bodyMap.at( body )->getGravityFieldModel( ) != NULL )
                    {
                        std::cerr<<"Warning when making time-dependent gravity field model for body "<<body<<" existing gravity field "
                                <<" is not empty but overwritten in Body! "<<std::endl;
                    }

                    // Create preliminary TimeDependentSphericalHarmonicsGravityField, without actual variation settings.
                    gravityFieldModel = boost::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( ) );
                }


            }
        }
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, did not recognize gravity field model settings type " +
                    boost::lexical_cast< std::string >(
                        gravityFieldSettings->getGravityFieldType( ) ) );
    }

    return gravityFieldModel;
}

} // namespace simulation_setup

} // namespace tudat
