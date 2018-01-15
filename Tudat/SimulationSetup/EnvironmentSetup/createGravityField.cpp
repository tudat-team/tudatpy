/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "Tudat/Astrodynamics/Gravitation/triAxialEllipsoidGravity.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace simulation_setup
{

//! Get the path of the SH file for a SH model.
std::string getPathForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel )
{
    switch ( sphericalHarmonicsModel )
    {
    case egm96:
        return input_output::getGravityModelsPath( ) + "Earth/egm96.txt";
    case ggm02c:
        return input_output::getGravityModelsPath( ) + "Earth/ggm02c.txt";
    case ggm02s:
        return input_output::getGravityModelsPath( ) + "Earth/ggm02s.txt";
    case glgm3150:
        return input_output::getGravityModelsPath( ) + "Moon/glgm3150.txt";
    case lpe200:
        return input_output::getGravityModelsPath( ) + "Moon/lpe200.txt";
    case jgmro120d:
        return input_output::getGravityModelsPath( ) + "Mars/jgmro120d.txt";
    default:
        std::cerr << "No path known for Spherical Harmonics Model " << sphericalHarmonicsModel << std::endl;
        throw;
    }
}

//! Get the associated reference frame for a SH model.
std::string getReferenceFrameForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel )
{
    switch ( sphericalHarmonicsModel )
    {
    case egm96:
    case ggm02c:
    case ggm02s:
        return "IAU_Earth";
    case glgm3150:
    case lpe200:
        return "IAU_Moon";
    case jgmro120d:
        return "IAU_Mars";
    default:
        std::cerr << "No reference frame known for Spherical Harmonics Model " << sphericalHarmonicsModel << std::endl;
        throw;
    }
}

//! Constructor with custom model.
FromFileSphericalHarmonicsGravityFieldSettings::FromFileSphericalHarmonicsGravityFieldSettings(
        const std::string& filePath, const std::string& associatedReferenceFrame,
        const int maximumDegree, const int maximumOrder,
        const int gravitationalParameterIndex, const int referenceRadiusIndex,
        const double gravitationalParameter, const double referenceRadius ) :
    SphericalHarmonicsGravityFieldSettings( gravitationalParameter, referenceRadius, Eigen::MatrixXd( ),
                                            Eigen::MatrixXd( ), associatedReferenceFrame ),
    filePath_( filePath ),
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    gravitationalParameterIndex_( gravitationalParameterIndex ),
    referenceRadiusIndex_( referenceRadiusIndex )
{
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;
    std::pair< double, double > referenceData =
            readGravityFieldFile( filePath, maximumDegree, maximumOrder, coefficients,
                                  gravitationalParameterIndex, referenceRadiusIndex );
    gravitationalParameter_ = gravitationalParameterIndex >= 0 ? referenceData.first : gravitationalParameter;
    referenceRadius_ = referenceRadiusIndex >= 0 ? referenceData.second : referenceRadius;
    cosineCoefficients_ = coefficients.first;
    sineCoefficients_ = coefficients.second;
}

//! Constructor with model included in Tudat.
FromFileSphericalHarmonicsGravityFieldSettings::FromFileSphericalHarmonicsGravityFieldSettings(
        const SphericalHarmonicsModel sphericalHarmonicsModel ) :
    FromFileSphericalHarmonicsGravityFieldSettings( getPathForSphericalHarmonicsModel( sphericalHarmonicsModel ),
                                                 getReferenceFrameForSphericalHarmonicsModel( sphericalHarmonicsModel ),
                                                 50, 50, 0, 1 )
{
    sphericalHarmonicsModel_ = sphericalHarmonicsModel;
}


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
        throw std::runtime_error( "Pds gravity field data file could not be opened: " + fileName );
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

        gravitationalParameter = std::stod( vectorOfIndividualStrings[ gravitationalParameterIndex ] );
        referenceRadius = std::stod( vectorOfIndividualStrings[ referenceRadiusIndex ] );
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
            std::string errorMessage = "Error when reading pds gravity field file, number of fields is " +
                    std::to_string( vectorOfIndividualStrings.size( ) );
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Read current degree and orde from line.
            currentDegree = std::stoi( vectorOfIndividualStrings[ 0 ] );
            currentOrder = std::stoi( vectorOfIndividualStrings[ 1 ] );

            // Set cosine and sine coefficients for current degree and order.
            if( currentDegree <= maximumDegree && currentOrder <= maximumOrder )
            {
                cosineCoefficients( currentDegree, currentOrder ) =
                        std::stod( vectorOfIndividualStrings[ 2 ] );
                sineCoefficients( currentDegree, currentOrder ) =
                        std::stod( vectorOfIndividualStrings[ 3 ] );
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
            throw std::runtime_error( "Error, requested central gravity field, but field variations settings are not empty." );
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
            throw std::runtime_error( "Error, requested central gravity field, but field variations settings are not empty." );
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
                        std::string errorMessage = "Warning when making time-dependent gravity field model for body " + body +
                                " existing gravity field is not empty but overwritten in Body! ";
                        throw std::runtime_error( errorMessage );
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
                    std::to_string(
                        gravityFieldSettings->getGravityFieldType( ) ) );
    }

    return gravityFieldModel;
}

//! Function to create gravity field settings for a homogeneous triaxial ellipsoid
boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettings(
        const double axisA, const double axisB, const double axisC, const double ellipsoidDensity,
        const int maximumDegree, const int maximumOrder,
        const std::string& associatedReferenceFrame  )
{
    // Compute reference quantities
    double ellipsoidGravitationalParameter = gravitation::calculateTriAxialEllipsoidVolume(
                axisA, axisB, axisC ) * ellipsoidDensity * physical_constants::GRAVITATIONAL_CONSTANT;
    double ellipsoidReferenceRadius = gravitation::calculateTriAxialEllipsoidVolume(
                axisA, axisB, axisC );

    // Compute coefficients
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients =
            gravitation::createTriAxialEllipsoidNormalizedSphericalHarmonicCoefficients(
                axisA, axisB, axisC, maximumDegree, maximumOrder );

    return boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                ellipsoidGravitationalParameter, ellipsoidReferenceRadius, coefficients.first,
                coefficients.second, associatedReferenceFrame );
}

} // namespace simulation_setup

} // namespace tudat
