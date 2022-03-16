/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/triAxialEllipsoidGravity.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/io/basicInputOutput.h"

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
        return paths::getGravityModelsPath( ) + "/Earth/egm96.txt";
    case ggm02c:
        return paths::getGravityModelsPath( ) + "/Earth/ggm02c.txt";
    case ggm02s:
        return paths::getGravityModelsPath( ) + "/Earth/ggm02s.txt";
    case goco05c:
        return paths::getGravityModelsPath( ) + "/Earth/GOCO05c.txt";
    case glgm3150:
        return paths::getGravityModelsPath( ) + "/Moon/glgm3150.txt";
    case gggrx1200:
        return paths::getGravityModelsPath( ) + "/Moon/gggrx_1200l_sha.tab";
    case lpe200:
        return paths::getGravityModelsPath( ) + "/Moon/lpe200.txt";
    case jgmro120d:
        return paths::getGravityModelsPath( ) + "/Mars/jgmro120d.txt";
    default:
        std::cerr << "No path known for Spherical Harmonics Model " << sphericalHarmonicsModel << std::endl;
        throw;
    }
}

int getMaximumGravityFieldDegreeOrder( const SphericalHarmonicsModel sphericalHarmonicsModel )
{
    int maximumDegreeOrder = 0;
    switch ( sphericalHarmonicsModel )
    {
    case egm96:
        maximumDegreeOrder = 360;
        break;
    case ggm02c:
        maximumDegreeOrder = 200;
        break;
    case ggm02s:
        maximumDegreeOrder = 160;
        break;
    case goco05c:
        maximumDegreeOrder = 719;
        break;
    case glgm3150:
        maximumDegreeOrder = 150;
        break;
    case gggrx1200:
        maximumDegreeOrder = 1199;
        break;
    case lpe200:
        maximumDegreeOrder = 200;
        break;
    case jgmro120d:
        maximumDegreeOrder = 120;
        break;
    default:
        throw std::runtime_error( "No maximum degree known for Spherical Harmonics Model " + std::to_string(
                                      static_cast< int >( sphericalHarmonicsModel ) ) );
    }
    return maximumDegreeOrder;
}

//! Get the associated reference frame for a SH model.
std::string getReferenceFrameForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel )
{
    switch ( sphericalHarmonicsModel )
    {
    case egm96:
    case ggm02c:
    case ggm02s:
    case goco05c:
        return "IAU_Earth";
    case glgm3150:
    case gggrx1200:
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
                                                    getMaximumGravityFieldDegreeOrder( sphericalHarmonicsModel ),
                                                    getMaximumGravityFieldDegreeOrder( sphericalHarmonicsModel ), 0, 1 )
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
            ( currentDegree <= maximumDegree && currentOrder < maximumOrder )  )
    {
        // Read current line
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Split string into multiple strings, each containing one element from a line from the
        // data file.
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( ", \t" ),
                                 boost::algorithm::token_compress_on );

        // Check current line for consistency
        if( vectorOfIndividualStrings.size( ) != 0 )
        {
            if( vectorOfIndividualStrings.size( ) < 4 )
            {
                std::string errorMessage = "Error when reading pds gravity field file, number of fields is " +
                        std::to_string( vectorOfIndividualStrings.size( ) );
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Read current degree and orde from line.
                currentDegree = static_cast< int >( std::round( std::stod( vectorOfIndividualStrings[ 0 ] ) ) );
                currentOrder = static_cast< int >( std::round( std::stod( vectorOfIndividualStrings[ 1 ] ) ) );
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
    }

    // Set cosine coefficient at (0,0) to 1.
    cosineCoefficients( 0, 0 ) = 1.0;
    coefficients = std::make_pair( cosineCoefficients, sineCoefficients );

    return std::make_pair( gravitationalParameter, referenceRadius );
}

//! Function to create a gravity field model.
std::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const std::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings )
{
    using namespace tudat::gravitation;

    // Declare return object.
    std::shared_ptr< GravityFieldModel > gravityFieldModel;

    // Check which type of gravity field model is to be created.
    switch( gravityFieldSettings->getGravityFieldType( ) )
    {
    case central:
    {
        // Check whether settings for point mass gravity field model are consistent with its type.
        std::shared_ptr< CentralGravityFieldSettings > centralFieldSettings =
                std::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        if( centralFieldSettings == nullptr )
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
            gravityFieldModel = std::make_shared< GravityFieldModel >(
                        centralFieldSettings->getGravitationalParameter( ) );
        }
        break;
    }
    case central_spice:
    {
        if( gravityFieldVariationSettings.size( ) != 0 )
        {
            throw std::runtime_error( "Error, requested central gravity field, but field variations settings are not empty." );
        }
        else
        {
            // Create and initialize point mass gravity field model from Spice.
            gravityFieldModel = std::make_shared< GravityFieldModel >(
                        spice_interface::getBodyGravitationalParameter( body ) );
        }

        break;
    }
    case spherical_harmonic:
    {
        // Check whether settings for spherical harmonic gravity field model are consistent with
        // its type.
        std::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicFieldSettings =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                    gravityFieldSettings );

        if( sphericalHarmonicFieldSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected spherical harmonic field settings when making gravity field model of "
                        + body );
        }
        else
        {
            std::function< void( ) > inertiaTensorUpdateFunction;
            if( bodies.count( body ) == 0 )
            {
                inertiaTensorUpdateFunction = std::function< void( ) >( );
            }
            else
            {
                inertiaTensorUpdateFunction =
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( body ), true );
                if( sphericalHarmonicFieldSettings->getScaledMeanMomentOfInertia( ) == sphericalHarmonicFieldSettings->getScaledMeanMomentOfInertia( ) )
                {
                    bodies.at( body )->setBodyInertiaTensor(
                                sphericalHarmonicFieldSettings->getInertiaTensor( ),
                                sphericalHarmonicFieldSettings->getScaledMeanMomentOfInertia( )) ;

                }
            }

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
                std::string associatedReferenceFrame = sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( );
                if( associatedReferenceFrame == "" )
                {
                    std::shared_ptr< ephemerides::RotationalEphemeris> rotationalEphemeris =
                            bodies.at( body )->getRotationalEphemeris( );
                    if( rotationalEphemeris == nullptr )
                    {
                        throw std::runtime_error( "Error when creating spherical harmonic gravity field for body " + body +
                                                  ", neither a frame ID nor a rotational model for the body have been defined" );
                    }
                    else
                    {
                        associatedReferenceFrame = rotationalEphemeris->getTargetFrameOrientation( );
                    }
                }

                if( gravityFieldVariationSettings.size( ) == 0 &&
                        sphericalHarmonicFieldSettings->getCreateTimeDependentField( ) == 0 )
                {
                    // Create and initialize spherical harmonic gravity field model.
                    gravityFieldModel = std::make_shared< SphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                associatedReferenceFrame,
                                inertiaTensorUpdateFunction );
                }
                else
                {
                    if( bodies.at( body )->getGravityFieldModel( ) != nullptr )
                    {
                        std::string errorMessage = "Warning when making time-dependent gravity field model for body " + body +
                                " existing gravity field is not empty but overwritten in Body! ";
                        throw std::runtime_error( errorMessage );
                    }

                    // Create preliminary TimeDependentSphericalHarmonicsGravityField, without actual variation settings.
                    gravityFieldModel = std::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                associatedReferenceFrame,
                                inertiaTensorUpdateFunction );
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
std::shared_ptr< SphericalHarmonicsGravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettings(
        const double axisA, const double axisB, const double axisC, const double ellipsoidDensity,
        const int maximumDegree, const int maximumOrder,
        const std::string& associatedReferenceFrame  )
{
    // Compute reference quantities
    double ellipsoidGravitationalParameter = gravitation::calculateTriAxialEllipsoidVolume(
                axisA, axisB, axisC ) * ellipsoidDensity * physical_constants::GRAVITATIONAL_CONSTANT;
    double ellipsoidReferenceRadius = gravitation::calculateTriAxialEllipsoidReferenceRadius(
                axisA, axisB, axisC );

    // Compute coefficients
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients =
            gravitation::createTriAxialEllipsoidNormalizedSphericalHarmonicCoefficients(
                axisA, axisB, axisC, maximumDegree, maximumOrder );

    return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                ellipsoidGravitationalParameter, ellipsoidReferenceRadius, coefficients.first,
                coefficients.second, associatedReferenceFrame );
}

} // namespace simulation_setup

} // namespace tudat
