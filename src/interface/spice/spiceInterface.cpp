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

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/spice/spiceException.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/paths.hpp"

#include <math.h>

namespace tudat
{

namespace spice_interface
{

std::string getCorrectedTargetBodyName( const std::string &targetBodyName )
{
    std::string correctedTargetBodyName;
    if( targetBodyName == "Mercury" || targetBodyName == "Venus" || targetBodyName == "MERCURY" || targetBodyName == "VENUS" ||
        targetBodyName == "mercury" || targetBodyName == "venus" )
    {
        correctedTargetBodyName = targetBodyName + "_Barycenter";
    }
    else
    {
        correctedTargetBodyName = targetBodyName;
    }

    return correctedTargetBodyName;
}

//! Convert a Julian date to ephemeris time (equivalent to TDB in Spice).
double convertJulianDateToEphemerisTime( const double julianDate )
{
    return ( julianDate - j2000_c( ) ) * spd_c( );
}

double getApproximateUtcFromTdb( const double ephemerisTime )
{
    double timeOffset = TUDAT_NAN;
    deltet_c( ephemerisTime, "ET", &timeOffset );
    return ephemerisTime - timeOffset;
}

//! Convert ephemeris time (equivalent to TDB) to a Julian date.
double convertEphemerisTimeToJulianDate( const double ephemerisTime )
{
    return j2000_c( ) + ( ephemerisTime ) / spd_c( );
}

//! Converts a date string to ephemeris time.
double convertDateStringToEphemerisTime( const std::string &dateString )
{
    double ephemerisTime = 0.0;
    str2et_c( dateString.c_str( ), &ephemerisTime );
    return ephemerisTime;
}

//! Get Cartesian state of a body, as observed from another body.
Eigen::Vector6d getBodyCartesianStateAtEpoch( const std::string &targetBodyName,
                                              const std::string &observerBodyName,
                                              const std::string &referenceFrameName,
                                              const std::string &aberrationCorrections,
                                              const double ephemerisTime )
{
    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving Cartesian state from Spice, input time is " + std::to_string( ephemerisTime ) );
    }
    // Declare variables for cartesian state and light-time to be determined by Spice.
    double stateAtEpoch[ 6 ];
    double lightTime;

    // Call Spice function to calculate state and light-time.
    spkezr_c( getCorrectedTargetBodyName( targetBodyName ).c_str( ),
              ephemerisTime,
              referenceFrameName.c_str( ),
              aberrationCorrections.c_str( ),
              getCorrectedTargetBodyName( observerBodyName ).c_str( ),
              stateAtEpoch,
              &lightTime );

    // Put result in Eigen Vector.
    Eigen::Vector6d cartesianStateVector;
    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 6; i++ )
        {
            cartesianStateVector( i ) = stateAtEpoch[ i ];
        }
    }
    else
    {
        handleSpiceException( );
        // cartesianStateVector.setConstant( 1.0E12 );
    }

    // Convert from km(/s) to m(/s).
    return unit_conversions::convertKilometersToMeters< Eigen::Vector6d >( cartesianStateVector );
}

//! Get Cartesian position of a body, as observed from another body.
Eigen::Vector3d getBodyCartesianPositionAtEpoch( const std::string &targetBodyName,
                                                 const std::string &observerBodyName,
                                                 const std::string &referenceFrameName,
                                                 const std::string &aberrationCorrections,
                                                 const double ephemerisTime )
{
    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving Cartesian position from Spice, input time is " +
                                     std::to_string( ephemerisTime ) );
    }
    // Declare variables for cartesian position and light-time to be determined by Spice.
    double positionAtEpoch[ 3 ];
    double lightTime;

    // Call Spice function to calculate position and light-time.
    spkpos_c( getCorrectedTargetBodyName( targetBodyName ).c_str( ),
              ephemerisTime,
              referenceFrameName.c_str( ),
              aberrationCorrections.c_str( ),
              getCorrectedTargetBodyName( observerBodyName ).c_str( ),
              positionAtEpoch,
              &lightTime );

    // Put result in Eigen Vector.
    Eigen::Vector3d cartesianPositionVector;

    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 3; i++ )
        {
            cartesianPositionVector( i ) = positionAtEpoch[ i ];
        }
    }
    else
    {
        handleSpiceException( );
        cartesianPositionVector.setConstant( 1.0E12 );
    }

    // Convert from km to m.
    return unit_conversions::convertKilometersToMeters< Eigen::Vector3d >( cartesianPositionVector );
}

//! Get Cartesian state of a satellite from its two-line element set at a specified epoch.
Eigen::Vector6d getCartesianStateFromTleAtEpoch( double epoch, std::shared_ptr< ephemerides::Tle > tle )
{
    if( !( epoch == epoch ) )
    {
        throw std::invalid_argument( "Error when retrieving TLE from Spice, input time is " + std::to_string( epoch ) );
    }

    // Physical constants used by CSpice's implementation of SGP4.
    double physicalConstants[ 8 ] = { 1.082616E-3, -2.53881E-6, -1.65597E-6, 7.43669161e-2, 120.0, 78.0, 6378.135, 1.0 };

    // Declare variable that will hold the state as returned by Spice.
    double stateAtEpoch[ 6 ];

    // TODO: convert elements to units required by CSpice (?)
    double elements[ 10 ];
    elements[ 0 ] = 0.0;  // This element is mandatory as input to ev2lin_ but not used internally (used to be accessed in SGP).
    elements[ 1 ] = 0.0;  // Idem dito.
    elements[ 2 ] = tle->getBStar( );
    elements[ 3 ] = tle->getInclination( );
    elements[ 4 ] = tle->getRightAscension( );
    elements[ 5 ] = tle->getEccentricity( );
    elements[ 6 ] = tle->getArgOfPerigee( );
    elements[ 7 ] = tle->getMeanAnomaly( );
    elements[ 8 ] = tle->getMeanMotion( );
    elements[ 9 ] = tle->getEpoch( );  // TLE ephemeris epoch in seconds since J2000

    // Call Spice function. Return value is always 0, so no need to save it.
    ev2lin_( &epoch, physicalConstants, elements, stateAtEpoch );

    // Put result in Eigen Vector.
    Eigen::Vector6d cartesianStateVector;
    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 6; i++ )
        {
            cartesianStateVector( i ) = stateAtEpoch[ i ];
        }
    }
    else
    {
        handleSpiceException( );
        cartesianStateVector.setConstant( 1.0E12 );
    }

    // Convert from km to m.
    return unit_conversions::convertKilometersToMeters< Eigen::Vector6d >( cartesianStateVector );
}

//! Compute quaternion of rotation between two frames.
Eigen::Quaterniond computeRotationQuaternionBetweenFrames( const std::string &originalFrame,
                                                           const std::string &newFrame,
                                                           const double ephemerisTime )
{
    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving rotation quaternion from Spice, input time is " +
                                     std::to_string( ephemerisTime ) );
    }

    // Declare rotation matrix.
    double rotationArray[ 3 ][ 3 ];

    // Calculate rotation matrix.
    pxform_c( originalFrame.c_str( ), newFrame.c_str( ), ephemerisTime, rotationArray );

    // Put rotation matrix in Eigen Matrix3d.
    Eigen::Matrix3d rotationMatrix;
    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                rotationMatrix( i, j ) = rotationArray[ i ][ j ];
            }
        }
    }
    else
    {
        handleSpiceException( );
        rotationMatrix.setIdentity( );
    }

    // Convert matrix3d to Quaternion.
    return Eigen::Quaterniond( rotationMatrix );
}

Eigen::Matrix3d computeRotationMatrixBetweenFrames( const std::string &originalFrame,
                                                    const std::string &newFrame,
                                                    const double ephemerisTime )
{
    return Eigen::Matrix3d( computeRotationQuaternionBetweenFrames( originalFrame, newFrame, ephemerisTime ) );
}

//! Compute rotation matrix for state vector between two frames.
Eigen::Matrix6d computeStateRotationMatrixBetweenFrames( const std::string &originalFrame,
                                                         const std::string &newFrame,
                                                         const double ephemerisTime )
{
    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving state rotation matrix from Spice, input time is " +
                                     std::to_string( ephemerisTime ) );
    }

    double stateTransition[ 6 ][ 6 ];

    // Calculate state transition matrix.
    sxform_c( originalFrame.c_str( ), newFrame.c_str( ), ephemerisTime, stateTransition );

    // Put rotation matrix in Eigen Matrix6d
    Eigen::Matrix6d stateTransitionMatrix = Eigen::Matrix6d::Zero( );
    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; j < 6; j++ )
            {
                stateTransitionMatrix( i, j ) = stateTransition[ i ][ j ];
            }
        }
    }
    else
    {
        handleSpiceException( );
        stateTransitionMatrix.setIdentity( );
    }

    return stateTransitionMatrix;
}

//! Computes time derivative of rotation matrix between two frames.
Eigen::Matrix3d computeRotationMatrixDerivativeBetweenFrames( const std::string &originalFrame,
                                                              const std::string &newFrame,
                                                              const double ephemerisTime )
{
    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving rotation matrix derivative from Spice, input time is " +
                                     std::to_string( ephemerisTime ) );
    }

    double stateTransition[ 6 ][ 6 ];

    // Calculate state transition matrix.
    sxform_c( originalFrame.c_str( ), newFrame.c_str( ), ephemerisTime, stateTransition );

    // Put rotation matrix derivative in Eigen Matrix3d
    Eigen::Matrix3d matrixDerivative = Eigen::Matrix3d::Zero( );
    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                matrixDerivative( i, j ) = stateTransition[ i + 3 ][ j ];
            }
        }
    }
    else
    {
        handleSpiceException( );
        matrixDerivative.setZero( );
    }

    return matrixDerivative;
}

//! Computes the angular velocity of one frame w.r.t. to another frame.
Eigen::Vector3d getAngularVelocityVectorOfFrameInOriginalFrame( const std::string &originalFrame,
                                                                const std::string &newFrame,
                                                                const double ephemerisTime )
{
    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving angular velocity from Spice, input time is " +
                                     std::to_string( ephemerisTime ) );
    }

    double stateTransition[ 6 ][ 6 ];

    // Calculate state transition matrix.
    sxform_c( originalFrame.c_str( ), newFrame.c_str( ), ephemerisTime, stateTransition );

    double rotation[ 3 ][ 3 ];
    double angularVelocity[ 3 ];

    // Calculate angular velocity vector.
    xf2rav_c( stateTransition, rotation, angularVelocity );

    if( !failed_c( ) )
    {
        return ( Eigen::Vector3d( ) << angularVelocity[ 0 ], angularVelocity[ 1 ], angularVelocity[ 2 ] ).finished( );
    }
    else
    {
        handleSpiceException( );
        return Eigen::Vector3d::Zero( );
    }
}

std::pair< Eigen::Quaterniond, Eigen::Matrix3d > computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames(
        const std::string &originalFrame,
        const std::string &newFrame,
        const double ephemerisTime )
{
    double stateTransition[ 6 ][ 6 ];

    if( !( ephemerisTime == ephemerisTime ) )
    {
        throw std::invalid_argument( "Error when retrieving rotational state from Spice, input time is " +
                                     std::to_string( ephemerisTime ) );
    }

    sxform_c( originalFrame.c_str( ), newFrame.c_str( ), ephemerisTime, stateTransition );

    Eigen::Matrix3d matrixDerivative;
    Eigen::Matrix3d rotationMatrix;
    if( !failed_c( ) )
    {
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                rotationMatrix( i, j ) = stateTransition[ i ][ j ];
                matrixDerivative( i, j ) = stateTransition[ i + 3 ][ j ];
            }
        }
    }
    else
    {
        handleSpiceException( );
        rotationMatrix.setIdentity( );
        matrixDerivative.setZero( );
    }

    return std::make_pair( Eigen::Quaterniond( rotationMatrix ), matrixDerivative );
}

//! Get property of a body from Spice.
std::vector< double > getBodyProperties( const std::string &body, const std::string &property, const int maximumNumberOfValues )
{
    // Delcare variable in which raw result is to be put by Spice function.
    double *propertyArray = new double[ maximumNumberOfValues ];

    // Call Spice function to retrieve property.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str( ), property.c_str( ), maximumNumberOfValues, &numberOfReturnedParameters, propertyArray );

    // Put result in STL vector.
    std::vector< double > bodyProperties;
    bodyProperties.resize( numberOfReturnedParameters );
    for( int i = 0; i < numberOfReturnedParameters; i++ )
    {
        bodyProperties.at( i ) = propertyArray[ i ];
    }
    delete[] propertyArray;
    return bodyProperties;
}

//! Get gravitational parameter of a body.
double getBodyGravitationalParameter( const std::string &body )
{
    // Delcare variable in which raw result is to be put by Spice function.
    double gravitationalParameter[ 1 ];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str( ), "GM", 1, &numberOfReturnedParameters, gravitationalParameter );

    // Convert from km^3/s^2 to m^3/s^2
    return unit_conversions::convertKilometersToMeters< double >( unit_conversions::convertKilometersToMeters< double >(
            unit_conversions::convertKilometersToMeters< double >( gravitationalParameter[ 0 ] ) ) );
}

//! Get the (arithmetic) mean of the three principal axes of the tri-axial ellipsoid shape.
double getAverageRadius( const std::string &body )
{
    // Delcare variable in which raw result is to be put by Spice function.
    double radii[ 3 ];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str( ), "RADII", 3, &numberOfReturnedParameters, radii );

    // Compute average and convert from km to m.
    return unit_conversions::convertKilometersToMeters< double >( radii[ 0 ] + radii[ 1 ] + radii[ 2 ] ) / 3.0;
}

//! Get the (arithmetic) mean of the two equatorial axes of the tri-axial ellipsoid shape.
double getAverageEquatorialRadius( const std::string &body )
{
    // Declare variable in which raw result is to be put by Spice function.
    double radii[ 3 ];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str( ), "RADII", 3, &numberOfReturnedParameters, radii );

    // Compute average and convert from km to m.
    return unit_conversions::convertKilometersToMeters< double >( radii[ 0 ] + radii[ 1 ] ) / 2.0;
}

//! Get the polar radius of the tri-axial ellipsoid shape.
double getPolarRadius( const std::string &body )
{
    // Declare variable in which raw result is to be put by Spice function.
    double radii[ 3 ];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str( ), "RADII", 3, &numberOfReturnedParameters, radii );

    // Compute average and convert from km to m.
    return unit_conversions::convertKilometersToMeters< double >( radii[ 2 ] );
}

//! Convert a body name to its NAIF identification number.
int convertBodyNameToNaifId( const std::string &bodyName )
{
    // Convert body name to NAIF ID number.
    SpiceInt bodyNaifId;
    SpiceBoolean isIdFound;
    bods2c_c( bodyName.c_str( ), &bodyNaifId, &isIdFound );

    // Convert SpiceInt (typedef for long) to int and return.
    return static_cast< int >( bodyNaifId );
}

//! Convert a NAIF identification number to its body name.
std::string convertNaifIdToBodyName( int bodyNaifId )
{
    // Maximum SPICE name length is 32. Therefore, a name length of 33 is used (+1 for null terminator)
    SpiceChar bodyName[ 33 ];

    bodc2s_c( bodyNaifId, 33, bodyName );

    // Convert SpiceChar to std::string
    return static_cast< std::string >( bodyName );
}

//! Check if a certain property of a body is in the kernel pool.
bool checkBodyPropertyInKernelPool( const std::string &bodyName, const std::string &bodyProperty )
{
    // Convert body name to NAIF ID.
    const int naifId = convertBodyNameToNaifId( bodyName );

    // Determine if property is in pool.
    SpiceBoolean isPropertyInPool = bodfnd_c( naifId, bodyProperty.c_str( ) );
    return static_cast< bool >( isPropertyInPool );
}

//! Load a Spice kernel.
void loadSpiceKernelInTudat( const std::string &fileName )
{
    setSpiceErrorHandling( );

    furnsh_c( fileName.c_str( ) );
}

//! Get the amount of loaded Spice kernels.
int getTotalCountOfKernelsLoaded( )
{
    SpiceInt count;
    ktotal_c( "ALL", &count );
    return count;
}

//! Clear all Spice kernels.
void clearSpiceKernels( )
{
    kclear_c( );
}

//! Get all standard Spice kernels used in tudat.
std::vector< std::string > getStandardSpiceKernels( const std::vector< std::string > alternativeEphemerisKernels )
{
    std::vector< std::string > standardSpiceKernels;

    //    std::string kernelPath = paths::getSpiceKernelPath();
    //    standardSpiceKernels.push_back(kernelPath + "/pck00010.tpc");
    //    standardSpiceKernels.push_back(kernelPath + "/gm_de431.tpc");

    //    if (alternativeEphemerisKernels.size() == 0) {
    //        standardSpiceKernels.push_back(kernelPath + "/tudat_merged_spk_kernel.bsp");
    //    } else {
    //        for (unsigned int i = 0; i < alternativeEphemerisKernels.size(); i++) {
    //            standardSpiceKernels.push_back(alternativeEphemerisKernels.at(i));
    //        }
    //    }
    //    standardSpiceKernels.push_back(kernelPath + "/naif0012.tls");
    return standardSpiceKernels;
}

void loadStandardSpiceKernels( const std::vector< std::string > alternativeEphemerisKernels )
{
    std::string kernelPath = paths::getSpiceKernelPath( );
    loadSpiceKernelInTudat( kernelPath + "/pck00010.tpc" );
    //    loadSpiceKernelInTudat(kernelPath + "/gm_de431.tpc");
    loadSpiceKernelInTudat( kernelPath + "/inpop19a_TDB_m100_p100_spice.tpc" );
    loadSpiceKernelInTudat( kernelPath + "/NOE-4-2020.tpc" );
    loadSpiceKernelInTudat( kernelPath + "/NOE-5-2021.tpc" );
    loadSpiceKernelInTudat( kernelPath + "/NOE-6-2018-MAIN-v2.tpc" );

    if( alternativeEphemerisKernels.size( ) == 0 )
    {
        loadSpiceKernelInTudat( kernelPath + "/codes_300ast_20100725.bsp" );
        loadSpiceKernelInTudat( kernelPath + "/codes_300ast_20100725.tf" );
        loadSpiceKernelInTudat( kernelPath + "/inpop19a_TDB_m100_p100_spice.bsp" );
        loadSpiceKernelInTudat( kernelPath + "/NOE-4-2020.bsp" );
        loadSpiceKernelInTudat( kernelPath + "/NOE-5-2021.bsp" );
        loadSpiceKernelInTudat( kernelPath + "/NOE-6-2018-MAIN-v2.bsp" );
        loadSpiceKernelInTudat( kernelPath + "/juice_mat_crema_4_0_20220601_20330626_v01.bsp" );
    }
    else
    {
        for( unsigned int i = 0; i < alternativeEphemerisKernels.size( ); i++ )
        {
            loadSpiceKernelInTudat( alternativeEphemerisKernels.at( i ) );
        }
    }
    loadSpiceKernelInTudat( kernelPath + "/naif0012.tls" );
}

Eigen::Matrix3d getRotationFromJ2000ToEclipJ2000( )
{
    return spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "ECLIPJ2000", 0.0 ).toRotationMatrix( );
}

Eigen::Matrix3d getRotationFromEclipJ2000ToJ2000( )
{
    return spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "J2000", 0.0 ).toRotationMatrix( );
}

void toggleErrorReturn( )
{
    erract_c( "SET", 0, "RETURN" );
}

void toggleErrorAbort( )
{
    errdev_c( "SET", 0, "ABORT" );
}

void suppressErrorOutput( )
{
    errdev_c( "SET", 0, "NULL" );
}

std::string getErrorMessage( )
{
    if( failed_c( ) )
    {
        SpiceChar message[ 1841 ];
        getmsg_c( "LONG", 1841, message );
        return static_cast< std::string >( message );
    }
    else
    {
        return "";
    }
}

bool checkFailure( )
{
    if( failed_c( ) )
    {
        // reset_c( );
        return true;
    }
    else
    {
        return false;
    }
}

void setSpiceErrorHandling( )
{
    erract_c( "SET", 0, "RETURN" );
    errdev_c( "SET", 0, "NULL" );
}

void handleSpiceException( )
{
    SpiceChar shortMessage[ SPICE_ERROR_SMSGLN ];
    SpiceChar explanation[ SPICE_ERROR_XMSGLN ];
    SpiceChar longMessage[ SPICE_ERROR_LMSGLN ];
    SpiceChar traceback[ SPICE_ERROR_TRCLEN ];

    // Get the error messages from SPICE
    getmsg_c( "SHORT", SPICE_ERROR_SMSGLN, shortMessage );
    getmsg_c( "EXPLAIN", SPICE_ERROR_XMSGLN, explanation );
    getmsg_c( "LONG", SPICE_ERROR_LMSGLN, longMessage );

    // Get the traceback
    qcktrc_c( SPICE_ERROR_TRCLEN, traceback );

    auto spiceException = getExceptionFromShortMessage( shortMessage, explanation, longMessage, traceback );

    throw *spiceException;
    // throwSpiceException< spiceException >( shortMessage, explanation, longMessage, traceback );
}

// template< typename T >
// void throwSpiceException( const std::string &shortMessage,
//                           const std::string &explanation,
//                           const std::string &longMessage,
//                           const std::string &traceback )
// {
//     throw T( shortMessage, explanation, longMessage, traceback );
// }

std::unique_ptr< tudat::exceptions::SpiceError > getExceptionFromShortMessage( const std::string &shortMessage,
                                                                const std::string &explanation,
                                                                const std::string &longMessage,
                                                                const std::string &traceback )
{
    // clang-format off
    /* 
    Assuming you have the cspice source code in a directory called cspice,
    SPICE errors can be retrieved using the following commands:
    cd cspice
    grep -rhPo 'sigerr_\(\K"SPICE\([A-Za-z]+\)"' | sort | uniq
    */
    std::map< std::string, std::function<std::unique_ptr<tudat::exceptions::SpiceError>(const std::string&, const std::string&, const std::string&, const std::string&)> > spiceExceptions = {
    {"SPICE(ADDRESSOUTOFBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceADDRESSOUTOFBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AGENTLISTOVERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceAGENTLISTOVERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ALLGONE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceALLGONE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AMBIGTEMPL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceAMBIGTEMPL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ARRAYSIZEMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceARRAYSIZEMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ARRAYTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceARRAYTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AVALOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceAVALOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AXISUNDERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceAXISUNDERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADACTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADACTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADADDRESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADADDRESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADANGLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGLEUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADANGLEUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGRATEERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADANGRATEERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGULARRATE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADANGULARRATE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGULARRATEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADANGULARRATEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADARCHITECTURE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADARCHITECTURE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADARRAYSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADARRAYSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADATTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADATTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADATTRIBUTE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADATTRIBUTE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADATTRIBUTES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADATTRIBUTES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAUVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADAUVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAVFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADAVFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAVFRAMEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADAVFRAMEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAXIS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADAXIS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAXISLENGTH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADAXISLENGTH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAXISNUMBERS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADAXISNUMBERS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBLOCKSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADBLOCKSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBODYID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADBODYID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBORESIGHTSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADBORESIGHTSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBOUNDARY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADBOUNDARY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCATALOGFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCATALOGFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCENTERNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCENTERNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCHECKFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCHECKFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCKTYPESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCKTYPESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOARSEVOXSCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOARSEVOXSCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOLUMDECL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOLUMDECL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOLUMNCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOLUMNCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOLUMNDECL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOLUMNDECL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOMMENTAREA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOMMENTAREA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOMPNUMBER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOMPNUMBER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOORDBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOORDBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOORDSYS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCOORDSYS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCURVETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADCURVETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDAFTRANSFERFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDAFTRANSFERFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASCOMMENTAREA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDASCOMMENTAREA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASDIRECTORY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDASDIRECTORY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDASFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASTRANSFERFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDASTRANSFERFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATALINE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDATALINE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATAORDERTOKEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDATAORDERTOKEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATATYPEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDATATYPEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDEFAULTVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDEFAULTVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDESCRTIMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDESCRTIMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDIMENSION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDIMENSION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDIMENSIONS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDIMENSIONS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDIRECTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDIRECTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDOUBLEPRECISION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDOUBLEPRECISION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDOWNSAMPLINGTOL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADDOWNSAMPLINGTOL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADECCENTRICITY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADECCENTRICITY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADENDPOINTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADENDPOINTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADEULERANGLEUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADEULERANGLEUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFILEFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFILEFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFINEVOXELSCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFINEVOXELSCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFORMATSPECIFIER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFORMATSPECIFIER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAMECLASS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFRAMECLASS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAMECOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFRAMECOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAMESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFRAMESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFROMTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFROMTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFROMTIMESYSTEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFROMTIMESYSTEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFROMTIMETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADFROMTIMETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADGEOMETRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADGEOMETRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADGM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADGM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADHARDSPACE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADHARDSPACE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADHERMITDEGREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADHERMITDEGREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINDEX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINDEX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINITSTATE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINITSTATE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTDATALINE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINPUTDATALINE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTETTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINPUTETTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINPUTTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTUTCTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINPUTUTCTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINSTRUMENTID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINSTRUMENTID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINTEGER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADINTEGER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADKERNELTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADKERNELTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADKERNELVARTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADKERNELVARTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLAGRANGEDEGREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLAGRANGEDEGREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLATITUDEBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLATITUDEBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLATITUDERANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLATITUDERANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLATUSRECTUM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLATUSRECTUM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLEAPSECONDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLEAPSECONDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLIMBLOCUSMIX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLIMBLOCUSMIX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLINEPERRECCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLINEPERRECCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLISTFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLISTFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLONGITUDERANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADLONGITUDERANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMATRIX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMATRIX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMEANMOTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMEANMOTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMECCENTRICITY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMECCENTRICITY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMETHODSYNTAX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMETHODSYNTAX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMIDNIGHTTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMIDNIGHTTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMSEMIMAJOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMSEMIMAJOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMSOPQUATERNION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADMSOPQUATERNION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADNOFDIGITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADNOFDIGITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADNOFSTATES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADNOFSTATES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADNUMBEROFPOINTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADNUMBEROFPOINTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOBJECTID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOBJECTID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOBJECTNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOBJECTNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETANGLES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOFFSETANGLES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETANGUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOFFSETANGUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETAXESFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOFFSETAXESFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETAXISXYZ)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOFFSETAXISXYZ>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADORBITALPERIOD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADORBITALPERIOD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOUTPUTSPKTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOUTPUTSPKTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOUTPUTTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADOUTPUTTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPARTNUMBER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPARTNUMBER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPCKVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPCKVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPECCENTRICITY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPECCENTRICITY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPERIAPSEVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPERIAPSEVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPICTURE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPICTURE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPLATECOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPLATECOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPODLOCATION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPODLOCATION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPRECVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPRECVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPRIORITYSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADPRIORITYSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADQUATSIGN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADQUATSIGN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADQUATTHRESHOLD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADQUATTHRESHOLD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRADIUS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADRADIUS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRADIUSCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADRADIUSCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRATEFRAMEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADRATEFRAMEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRATETHRESHOLD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADRATETHRESHOLD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRECORDCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADRECORDCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADREFVECTORSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADREFVECTORSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTATIONAXISXYZ)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADROTATIONAXISXYZ>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTATIONSORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADROTATIONSORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTATIONTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADROTATIONTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTAXESFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADROTAXESFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROWCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADROWCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSCID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSCID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSEMIAXIS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSEMIAXIS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSEMILATUS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSEMILATUS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSHAPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSHAPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOLDAY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSOLDAY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOLINDEX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSOLINDEX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOLTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSOLTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOURCERADIUS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSOURCERADIUS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSPICEQUATERNION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSPICEQUATERNION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTARINDEX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSTARINDEX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTARTTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSTARTTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTDIONAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSTDIONAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTOPTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSTOPTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSUBSTR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSUBSTR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSUBSTRINGBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSUBSTRINGBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSURFACEMAP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADSURFACEMAP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTABLEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTABLEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTERMLOCUSMIX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTERMLOCUSMIX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMEBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMECASE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMECASE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMECOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMECOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMEFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEITEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMEITEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEOFFSET)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMEOFFSET>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMESTRING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMESTRING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMETYPEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTIMETYPEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTLECOVERAGEPAD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTLECOVERAGEPAD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTLEPADS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTLEPADS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTOTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTOTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTOTIMESYSTEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTOTIMESYSTEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTOTIMETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTOTIMETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTYPESHAPECOMBO)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADTYPESHAPECOMBO>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARASSIGN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVARASSIGN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARIABLESIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVARIABLESIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARIABLETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVARIABLETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVARNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVECTOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVECTOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVERTEXCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVERTEXCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVERTEXINDEX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADVERTEXINDEX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADWINDOWSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBADWINDOWSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BARRAYTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBARRAYTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BARYCENTEREPHEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBARYCENTEREPHEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BARYCENTERIDCODE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBARYCENTERIDCODE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BEFOREBEGSTR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBEFOREBEGSTR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKCOMMANDLINE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKCOMMANDLINE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKFILETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKFILETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKINPUTFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKINPUTFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKINPUTTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKINPUTTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKNAMEASSIGNED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKNAMEASSIGNED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKOUTPTFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKOUTPTFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKSCLKSTRING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKSCLKSTRING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKTIMEFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLANKTIMEFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLOCKSNOTEVEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBLOCKSNOTEVEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BODIESNOTDISTINCT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBODIESNOTDISTINCT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BODYANDCENTERSAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBODYANDCENTERSAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOGUSENTRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBOGUSENTRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BORESIGHTMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBORESIGHTMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDARYMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBOUNDARYMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDARYTOOBIG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBOUNDARYTOOBIG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDSDISAGREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBOUNDSDISAGREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDSOUTOFORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBOUNDSOUTOFORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUFFEROVERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBUFFEROVERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUFFERSIZESMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBUFFERSIZESMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUFFERTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBUFFERTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBUG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUGWRITEFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceBUGWRITEFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CALLCKBSSFIRST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCALLCKBSSFIRST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CALLEDOUTOFORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCALLEDOUTOFORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CALLZZDSKBSSFIRST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCALLZZDSKBSSFIRST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTFINDGRP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANNOTFINDGRP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTGETPACKET)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANNOTGETPACKET>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTMAKEFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANNOTMAKEFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTPICKFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANNOTPICKFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANTFINDFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANTFINDFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANTGETROTATIONTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANTGETROTATIONTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANTUSEPERIAPEPOCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCANTUSEPERIAPEPOCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CBNOSUCHSTR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCBNOSUCHSTR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CELLARRAYTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCELLARRAYTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CELLTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCELLTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKBOGUSENTRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKBOGUSENTRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKDOESNTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKDOESNTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKNONEXISTREC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKNONEXISTREC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKTOOMANYFILES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKTOOMANYFILES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKUNKNOWNDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKUNKNOWNDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKWRONGDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCKWRONGDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CMDERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCMDERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CMDPARSEERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCMDPARSEERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COARSEGRIDOVERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOARSEGRIDOVERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COLDESCTABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOLDESCTABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COLUMNTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOLUMNTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMMANDTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOMMANDTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMMENTTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOMMENTTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMMFILENOTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOMMFILENOTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMPETINGEPOCHSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOMPETINGEPOCHSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMPETINGFRAMESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOMPETINGFRAMESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COORDSYSNOTREC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOORDSYSNOTREC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COUNTMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOUNTMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COUNTTOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOUNTTOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COVERAGEGAP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCOVERAGEGAP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CROSSANGLEMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceCROSSANGLEMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFBADCRECLEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFBADCRECLEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFBEGGTEND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFBEGGTEND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFCRNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFCRNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFDPWRITEFAIL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFDPWRITEFAIL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFFRNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFFRNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFFTFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFFTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFILLEGWRITE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFILLEGWRITE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFINVALIDACCESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFINVALIDACCESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFINVALIDPARAMS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFINVALIDPARAMS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNEGADDR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNEGADDR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNEWCONFLICT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNEWCONFLICT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOIFNMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOIFNMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNONAMEMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNONAMEMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNORESV)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNORESV>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSEARCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOSEARCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHADDR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOSUCHADDR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOSUCHFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHHANDLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOSUCHHANDLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHUNIT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOSUCHUNIT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOWRITE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFNOWRITE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFOVERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFOVERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFREADFAIL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFREADFAIL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFWRITEFAIL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDAFWRITEFAIL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASFILEREADFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASFILEREADFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASFILEWRITEFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASFILEWRITEFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASFTFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASFTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASINVALIDACCESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASINVALIDACCESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASINVALIDCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASINVALIDCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASINVALIDTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASINVALIDTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHADDRESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASNOSUCHADDRESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASNOSUCHFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHHANDLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASNOSUCHHANDLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHUNIT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASNOSUCHUNIT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOTEMPTY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASNOTEMPTY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASREADFAIL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASREADFAIL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASWRITEFAIL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDASWRITEFAIL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATAITEMLIMITEXCEEDED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDATAITEMLIMITEXCEEDED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATAREADFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDATAREADFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATAWIDTHERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDATAWIDTHERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATEEXPECTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDATEEXPECTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DECODINGERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDECODINGERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGENERATECASE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDEGENERATECASE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGENERATEINTERVAL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDEGENERATEINTERVAL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGENERATESURFACE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDEGENERATESURFACE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGREEOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDEGREEOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEPENDENTVECTORS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDEPENDENTVECTORS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEVICENAMETOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDEVICENAMETOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIFFLINETOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDIFFLINETOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIFFLINETOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDIFFLINETOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIMENSIONTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDIMENSIONTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DISARRAY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDISARRAY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DISORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDISORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIVIDEBYZERO)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDIVIDEBYZERO>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DSKBOGUSENTRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDSKBOGUSENTRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DSKDATANOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDSKDATANOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DSKTOOMANYFILES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDSKTOOMANYFILES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DTOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDTOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DUBIOUSMETHOD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDUBIOUSMETHOD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DUPLICATETIMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceDUPLICATETIMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ECCOUTOFBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceECCOUTOFBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ECCOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceECCOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKCOLATTRTABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKCOLATTRTABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKCOLNUMMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKCOLNUMMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKFILETABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKFILETABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKMISSINGCOLUMN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKMISSINGCOLUMN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKNOSEGMENTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKNOSEGMENTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKSEGTABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKSEGTABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKTABLELISTFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEKTABLELISTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ELEMENTSTOOSHORT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceELEMENTSTOOSHORT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EMPTYINPUTFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEMPTYINPUTFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EMPTYSEGMENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEMPTYSEGMENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ENDOFFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceENDOFFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ENDPOINTSMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceENDPOINTSMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ERROREXIT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceERROREXIT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EVECOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEVECOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EVENHERMITDEGREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEVENHERMITDEGREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EVILBOGUSENTRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEVILBOGUSENTRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EXTERNALOPEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceEXTERNALOPEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FACENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFACENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FAKESCLKEXISTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFAKESCLKEXISTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILARCHMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILARCHMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILARCMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILARCMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEALREADYEXISTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEALREADYEXISTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILECURRENTLYOPEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILECURRENTLYOPEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEDELETEFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEDELETEFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEDOESNOTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEDOESNOTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEEXISTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEEXISTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEISNOTSPK)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEISNOTSPK>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENAMETOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILENAMETOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENOTCONNECTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILENOTCONNECTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENOTOPEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILENOTOPEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENCONFLICT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEOPENCONFLICT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEOPENERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENFAIL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEOPENFAIL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEOPENFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEREADERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEREADERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEREADFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEREADFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILETABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILETABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILETRUNCATED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILETRUNCATED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEWRITEFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFILEWRITEFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FIRSTRECORDMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFIRSTRECORDMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FKDOESNTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFKDOESNTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FMTITEMLIMITEXCEEDED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFMTITEMLIMITEXCEEDED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATDATAMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFORMATDATAMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATDOESNTAPPLY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFORMATDOESNTAPPLY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFORMATERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATNOTAPPLICABLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFORMATNOTAPPLICABLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATSTRINGTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFORMATSTRINGTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FOVTOOWIDE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFOVTOOWIDE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEDATANOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMEDATANOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEDEFERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMEDEFERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEIDNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMEIDNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEINFONOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMEINFONOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMEMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMENAMENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMENAMENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMENOTRECOGNIZED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFRAMENOTRECOGNIZED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FTFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FTPXFERERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceFTPXFERERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(GRIDTOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceGRIDTOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(HANDLENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceHANDLENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(HASHISFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceHASHISFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(HLULOCKFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceHLULOCKFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IDCODENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceIDCODENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IDSTRINGTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceIDSTRINGTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGALCHARACTER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceILLEGALCHARACTER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGALOPTIONNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceILLEGALOPTIONNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGSHIFTDIR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceILLEGSHIFTDIR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGTEMPL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceILLEGTEMPL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IMMUTABLEVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceIMMUTABLEVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IMPROPERFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceIMPROPERFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IMPROPEROPEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceIMPROPEROPEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INACTIVEOBJECT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINACTIVEOBJECT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLEEOL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCOMPATIBLEEOL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLENUMREF)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCOMPATIBLENUMREF>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLESCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCOMPATIBLESCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLEUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCOMPATIBLEUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPLETEFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCOMPLETEFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTCENTERID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCONSISTCENTERID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTENTTIMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCONSISTENTTIMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCONSISTFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTSTARTTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCONSISTSTARTTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTSTOPTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCONSISTSTOPTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCORRECTUSAGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINCORRECTUSAGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDEFINITELOCALSECOND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINDEFINITELOCALSECOND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDEXOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINDEXOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDEXTOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINDEXTOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDICESOUTOFORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINDICESOUTOFORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTDOESNOTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINPUTDOESNOTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTFILENOTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINPUTFILENOTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTOUTOFBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINPUTOUTOFBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTSTOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINPUTSTOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INQUIREERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINQUIREERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INQUIREFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINQUIREFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSIDEBODY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINSIDEBODY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFFICIENTANGLES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINSUFFICIENTANGLES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFFICIENTDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINSUFFICIENTDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFFLEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINSUFFLEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFPTRSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINSUFPTRSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTERVALSTARTNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINTERVALSTARTNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTINDEXTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINTINDEXTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTLENNOTPOS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINTLENNOTPOS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINTOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDACCESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDACCESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDACTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDACTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDADD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDADD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDADDRESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDADDRESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDANGLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDANGLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDARCHTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDARCHTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDARGUMENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDARGUMENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDAXIS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDAXIS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDAXISLENGTH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDAXISLENGTH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCARDINALITY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDCARDINALITY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCASE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDCASE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCOLUMN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDCOLUMN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCONSTSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDCONSTSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDATACOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDATACOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDEGREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDEGREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDESCRTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDESCRTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDIMENSION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDIMENSION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDIRECTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDIRECTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDIVISOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDDIVISOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDELLIPSE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDELLIPSE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDENDPNTSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDENDPNTSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDENDPTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDENDPTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDEPOCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDEPOCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFILETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFILETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFIXREF)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFIXREF>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFOV)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFOV>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFRAMEDEF)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDFRAMEDEF>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDGEOMETRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDGEOMETRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDHANDLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDHANDLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDINDEX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDINDEX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDINTEGER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDINTEGER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLIMBTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDLIMBTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLISTITEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDLISTITEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLOCUS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDLOCUS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLONEXTENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDLONEXTENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDMETADATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDMETADATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDMETHOD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDMETHOD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDMSGTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDMSGTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNODE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNODE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMBEROFINTERVALS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNUMBEROFINTERVALS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMBEROFRECORDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNUMBEROFRECORDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMINT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNUMINT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMINTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNUMINTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMREC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDNUMREC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDOCCTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDOCCTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDOPERATION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDOPERATION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDOPTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDOPTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDPLANE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDPLANE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDRADII)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDRADII>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDRADIUS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDRADIUS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDREFFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDREFFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDREFVAL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDREFVAL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDROLLSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDROLLSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCLKRATE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSCLKRATE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCLKSTRING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSCLKSTRING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCLKTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSCLKTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSEARCHSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSEARCHSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSELECTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSELECTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSHADOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSHADOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSHAPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSHAPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSHAPECOMBO)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSHAPECOMBO>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTARTTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSTARTTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTATE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSTATE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTEPSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSTEPSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSUBLIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSUBLIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSUBTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDSUBTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTABLENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTABLENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTABLESIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTABLESIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTARGET)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTARGET>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTERMTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTERMTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTEXT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTEXT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTIMEFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTIMEFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTIMESTRING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTIMESTRING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTLEORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTLEORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTOLERANCE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTOLERANCE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDVERTEX)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVALIDVERTEX>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVERSTARTSTOPTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceINVERSTARTSTOPTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IRFNOTREC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceIRFNOTREC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ITEMNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceITEMNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ITEMNOTRECOGNIZED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceITEMNOTRECOGNIZED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ITERATIONEXCEEDED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceITERATIONEXCEEDED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERNELNOTLOADED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceKERNELNOTLOADED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERNELPOOLFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceKERNELPOOLFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERNELVARNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceKERNELVARNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERVARSETOVERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceKERVARSETOVERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERVARTOOBIG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceKERVARTOOBIG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KEYWORDNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceKEYWORDNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBCORRUPTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceLBCORRUPTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBLINETOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceLBLINETOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBNOSUCHLINE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceLBNOSUCHLINE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBTOOMANYLINES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceLBTOOMANYLINES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LOWERBOUNDTOOLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceLOWERBOUNDTOOLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LSKDOESNTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceLSKDOESNTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MALFORMEDSEGMENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMALFORMEDSEGMENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MARKERNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMARKERNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MESSAGETOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMESSAGETOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISMATCHFROMTIMETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISMATCHFROMTIMETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISMATCHOUTPUTFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISMATCHOUTPUTFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISMATCHTOTIMETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISMATCHTOTIMETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGARGUMENTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGARGUMENTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCENTER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGCENTER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCOLSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGCOLSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCOORDBOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGCOORDBOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCOORDSYS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGCOORDSYS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATACLASS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGDATACLASS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATAORDERTK)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGDATAORDERTK>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGEOT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGEOT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGEPOCHTOKEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGEPOCHTOKEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGFRAMEVAR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGFRAMEVAR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGGEOCONSTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGGEOCONSTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGHEIGHTREF)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGHEIGHTREF>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGHSCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGHSCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGKPV)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGKPV>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGLEFTCOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGLEFTCOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGLEFTRTFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGLEFTRTFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGNCAPFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGNCAPFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGNCOLS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGNCOLS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGNROWS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGNROWS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGPLATETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGPLATETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGROWMAJFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGROWMAJFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGROWSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGROWSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGSCAPFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGSCAPFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGSURFACE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGSURFACE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTIMEINFO)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGTIMEINFO>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTLEKEYWORD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGTLEKEYWORD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTOPCOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGTOPCOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTOPDOWNFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGTOPDOWNFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGVOXELSCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGVOXELSCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGWRAPFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceMISSINGWRAPFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NAMENOTUNIQUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNAMENOTUNIQUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NAMESNOTRESOLVED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNAMESNOTRESOLVED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NAMETABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNAMETABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NARATESFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNARATESFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NEGATIVETOL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNEGATIVETOL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOACCEPTABLEDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOACCEPTABLEDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOANGULARRATEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOANGULARRATEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOARRAYSTARTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOARRAYSTARTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOATTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOATTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOAVDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOAVDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOBODYID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOBODYID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCANDOSPKSPCKS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCANDOSPKSPCKS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCENTERIDORNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCENTERIDORNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCKSEGMENTTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCKSEGMENTTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCLASS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCLASS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCOMMENTSFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCOMMENTSFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCONVERG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCONVERG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCONVERGENCE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCONVERGENCE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCURRENTARRAY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOCURRENTARRAY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODATAORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNODATAORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODATATYPEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNODATATYPEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODELIMCHARACTER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNODELIMCHARACTER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODETOOFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNODETOOFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODSKSEGMENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNODSKSEGMENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODSKSEGMENTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNODSKSEGMENTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOENVVARIABLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOENVVARIABLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOEULERANGLEUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOEULERANGLEUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFILENAMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFILENAMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFILES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFILES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFILESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFILESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMECONNECT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFRAMECONNECT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMEDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFRAMEDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMEINFO)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFRAMEINFO>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFRAMENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMESKERNELNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFRAMESKERNELNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFREELOGICALUNIT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFREELOGICALUNIT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFREENODES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFREENODES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFROMTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFROMTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFROMTIMESYSTEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOFROMTIMESYSTEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOHEADNODE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOHEADNODE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINPUTDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOINPUTDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINPUTFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOINPUTFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINSTRUMENTID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOINSTRUMENTID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINTERVAL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOINTERVAL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOKERNELLOADED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOKERNELLOADED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLANDINGTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLANDINGTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLEAPSECONDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLEAPSECONDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLINESPERRECCOUNT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLINESPERRECCOUNT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLISTFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLISTFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLOADEDDSKFILES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLOADEDDSKFILES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLOADEDFILES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLOADEDFILES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLSKFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOLSKFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOMOREROOM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOMOREROOM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONCONICMOTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONCONICMOTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONDISTINCTPAIR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONDISTINCTPAIR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONEMPTYENTRY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONEMPTYENTRY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONEMPTYTREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONEMPTYTREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONEXISTELEMENTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONEXISTELEMENTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONINTEGERFIELD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONINTEGERFIELD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONNUMERICSTRING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONNUMERICSTRING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSBUFLENGTH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSBUFLENGTH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVEAXIS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSITIVEAXIS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVEMASS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSITIVEMASS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVERADIUS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSITIVERADIUS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVESCALE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSITIVESCALE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVEVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSITIVEVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSPACKETSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPOSPACKETSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPRINTABLECHARS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPRINTABLECHARS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPRINTINGCHAR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPRINTINGCHAR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPRINTINGCHARS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONPRINTINGCHARS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONUNITNORMAL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNONUNITNORMAL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOBJECTIDORNAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOOBJECTIDORNAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOFFSETANGLEAXES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOOFFSETANGLEAXES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOFFSETANGLEUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOOFFSETANGLEUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOUTPUTFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOOUTPUTFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOUTPUTSPKTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOOUTPUTSPKTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPARTITION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOPARTITION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPICTURE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOPICTURE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPOLYNOMIALDEGREE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOPOLYNOMIALDEGREE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPRECESSIONTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOPRECESSIONTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPRODUCERID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOPRODUCERID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOROTATIONORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOROTATIONORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSCID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSCID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSCLKFILENAMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSCLKFILENAMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSECONDLINE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSECONDLINE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSEGMENTSFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSEGMENTSFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSEPARATION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSEPARATION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSLKFILENAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSLKFILENAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSOLMARKER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSOLMARKER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSPACECRAFTID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSPACECRAFTID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSTARTTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSTARTTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSTOPTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSTOPTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUCHFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSUCHFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUCHHANDLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSUCHHANDLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUCHSYMBOL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSUCHSYMBOL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUNGM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOSUNGM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTABINARYKERNEL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTABINARYKERNEL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTACKFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTACKFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTADAFFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTADAFFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTADASFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTADASFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTADPNUMBER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTADPNUMBER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANDPNUMBER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTANDPNUMBER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANINTEGER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTANINTEGER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANINTEGERNUMBER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTANINTEGERNUMBER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANINTNUMBER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTANINTNUMBER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTAPCKFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTAPCKFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTAROTATION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTAROTATION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTATEXTFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTATEXTFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTATRANSFERFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTATRANSFERFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTCOMPUTABLE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTCOMPUTABLE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTDIMENSIONALLYEQUIV)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTDIMENSIONALLYEQUIV>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTDISJOINT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTDISJOINT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTDISTINCT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTDISTINCT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTENOUGHPEAS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTENOUGHPEAS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTIMETYPEFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTIMETYPEFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTINDEXED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTINDEXED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTINITIALIZED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTINITIALIZED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTINPART)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTINPART>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTLEDATAFOROBJECT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTLEDATAFOROBJECT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTLEGALCB)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTLEGALCB>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTOTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTOTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTOTIMESYSTEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTOTIMESYSTEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTRANSLATION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTRANSLATION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTRECOGNIZED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTRECOGNIZED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTSEMCHECKED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTSEMCHECKED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTSUPPORTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTSUPPORTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTTWOFIELDSCLK)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTTWOFIELDSCLK>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTTWOMODULI)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTTWOMODULI>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTTWOOFFSETS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOTTWOOFFSETS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOUNITSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNOUNITSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMBEREXPECTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNUMBEREXPECTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMCOEFFSNOTPOS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNUMCOEFFSNOTPOS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMERICOVERFLOW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNUMERICOVERFLOW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMPACKETSNOTPOS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNUMPACKETSNOTPOS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMPARTSUNEQUAL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNUMPARTSUNEQUAL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMSTATESNOTPOS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceNUMSTATESNOTPOS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OBJECTLISTFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOBJECTLISTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OBJECTSTOOCLOSE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOBJECTSTOOCLOSE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ORBITDECAY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceORBITDECAY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTOFPLACEDELIMITER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTOFPLACEDELIMITER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTOFROOM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTOFROOM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTFILEEXISTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTPUTFILEEXISTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTISNOTSPK)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTPUTISNOTSPK>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTPUTTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTTOOSHORT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceOUTPUTTOOSHORT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PARSERNOTREADY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePARSERNOTREADY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PARTIALFRAMESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePARTIALFRAMESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PASTENDSTR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePASTENDSTR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PATHMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePATHMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PATHTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePATHTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKDOESNTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePCKDOESNTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePCKFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKFILETABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePCKFILETABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKKRECTOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePCKKRECTOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PLATELISTTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePLATELISTTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTEROUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTEROUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTERSETTOOBIG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTERSETTOOBIG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTERTABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTERTABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTNOTINSEGMENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTNOTINSEGMENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTNOTONSURFACE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTNOTONSURFACE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTOFFSURFACE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTOFFSURFACE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTONZAXIS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTONZAXIS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePOINTTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PTRARRAYTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpicePTRARRAYTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(QPARAMOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceQPARAMOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(QUERYFAILURE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceQUERYFAILURE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(QUERYNOTPARSED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceQUERYNOTPARSED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RADIIOUTOFORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceRADIIOUTOFORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RAYISZEROVECTOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceRAYISZEROVECTOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(READFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceREADFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RECORDNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceRECORDNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RECURSIONTOODEEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceRECURSIONTOODEEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RECURSIVELOADING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceRECURSIVELOADING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REFANGLEMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceREFANGLEMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REFVALNOTINTEGER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceREFVALNOTINTEGER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REFVECTORMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceREFVECTORMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REQUESTOUTOFBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceREQUESTOUTOFBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REQUESTOUTOFORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceREQUESTOUTOFORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RWCONFLICT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceRWCONFLICT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SBINSUFPTRSIZE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSBINSUFPTRSIZE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SBTOOMANYSTRS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSBTOOMANYSTRS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SCLKDOESNTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSCLKDOESNTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SCLKTRUNCATED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSCLKTRUNCATED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGIDTOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSEGIDTOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGMENTNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSEGMENTNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGMENTTABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSEGMENTTABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGTABLETOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSEGTABLETOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGTYPECONFLICT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSEGTYPECONFLICT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SETEXCESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSETEXCESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SETTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSETTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SETUPDOESNOTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSETUPDOESNOTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SHAPEMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSHAPEMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SHAPENOTSUPPORTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSHAPENOTSUPPORTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SIZEMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSIZEMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SIZEOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSIZEOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPACETOONARROW)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPACETOONARROW>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPCRFLNOTCALLED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPCRFLNOTCALLED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPICEISTIRED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPICEISTIRED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKDOESNTEXIST)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKDOESNTEXIST>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKFILETABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKFILETABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKINSUFFDATA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKINSUFFDATA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKINVALIDOPTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKINVALIDOPTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKNOTASUBSET)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKNOTASUBSET>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKRECTOOLARGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKRECTOOLARGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKREFNOTSUPP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKREFNOTSUPP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKSTRUCTUREERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKSTRUCTUREERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKTYPENOTSUPP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKTYPENOTSUPP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKTYPENOTSUPPORTD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPKTYPENOTSUPPORTD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPURIOUSKEYWORD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSPURIOUSKEYWORD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STRINGTOOSHORT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSTRINGTOOSHORT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STRINGTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSTRINGTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STRINGTRUNCATED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSTRINGTRUNCATED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SUBORBITAL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSUBORBITAL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SUBPOINTNOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSUBPOINTNOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SYNTAXERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSYNTAXERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SYSTEMCALLFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceSYSTEMCALLFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TABLENOTLOADED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTABLENOTLOADED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMECONFLICT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTIMECONFLICT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMEOUTOFBOUNDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTIMEOUTOFBOUNDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMESDONTMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTIMESDONTMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMESOUTOFORDER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTIMESOUTOFORDER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMEZONEERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTIMEZONEERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWINPUTLINES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOFEWINPUTLINES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWPACKETS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOFEWPACKETS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWPLATES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOFEWPLATES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWSTATES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOFEWSTATES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWVERTICES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOFEWVERTICES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWWINDOWS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOFEWWINDOWS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYBASEFRAMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYBASEFRAMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYFIELDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYFIELDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYFILESOPEN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYFILESOPEN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYHITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYHITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYITERATIONS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYITERATIONS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYKEYWORDS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYKEYWORDS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPAIRS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYPAIRS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPARTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYPARTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPEAS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYPEAS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPLATES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYPLATES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYSURFACES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYSURFACES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYVERTICES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYVERTICES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYWATCHES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTOOMANYWATCHES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TRANSFERFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTRANSFERFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TRANSFERFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTRANSFERFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TWOSCLKFILENAMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTWOSCLKFILENAMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TYPEMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTYPEMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TYPENOTSUPPORTED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTYPENOTSUPPORTED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TYPESMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceTYPESMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNALLOCATEDNODE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNALLOCATEDNODE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNBALANCEDGROUP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNBALANCEDGROUP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNBALANCEDPAIR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNBALANCEDPAIR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNDEFINEDFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNDEFINEDFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNEQUALTIMESTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNEQUALTIMESTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNINITIALIZED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNINITIALIZED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNINITIALIZEDHASH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNINITIALIZEDHASH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNINITIALIZEDVALUE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNINITIALIZEDVALUE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNITSMISSING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNITSMISSING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNITSNOTREC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNITSNOTREC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNONWNTIMESYSTEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNONWNTIMESYSTEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNBFF)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNBFF>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNCKMETA)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNCKMETA>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNCOMPARE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNCOMPARE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFILARC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNFILARC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFRAMESPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNFRAMESPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFRAMETYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNFRAMETYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNINCLUSION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNINCLUSION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNINDEXTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNINDEXTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNKERNELTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNKERNELTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNKEY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNKEY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNMETAITEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNMETAITEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNMODE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNMODE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNOP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNOP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNPCKTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNPCKTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNREFDIR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNREFDIR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNSPKTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNSPKTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNSYSTEM)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNSYSTEM>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNUNITS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNKNOWNUNITS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNMATCHENDPTS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNMATCHENDPTS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNNATURALACT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNNATURALACT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNNATURALRELATION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNNATURALRELATION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNORDEREDREFS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNORDEREDREFS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNORDEREDTIMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNORDEREDTIMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNPARSEDQUERY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNPARSEDQUERY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNPARSEDTIME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNPARSEDTIME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNAPPFLAG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNAPPFLAG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNDELIMITER)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNDELIMITER>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZABLEFILE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNIZABLEFILE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDACTION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNIZEDACTION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDFORMAT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNIZEDFORMAT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDFRAME)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNIZEDFRAME>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNIZEDTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNPRECTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRECOGNPRECTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRESOLVEDNAMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRESOLVEDNAMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRESOLVEDTIMES)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNRESOLVEDTIMES>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDARCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNSUPPORTEDARCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDBFF)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNSUPPORTEDBFF>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDMETHOD)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNSUPPORTEDMETHOD>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDSPEC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNSUPPORTEDSPEC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNTITLEDHELP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUNTITLEDHELP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UPDATEPENDING)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUPDATEPENDING>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(USAGEERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUSAGEERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UTFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceUTFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VALUEOUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVALUEOUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VALUETABLEFULL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVALUETABLEFULL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VARIABLENOTFOUND)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVARIABLENOTFOUND>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VARNAMETOOLONG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVARNAMETOOLONG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VECTORTOOBIG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVECTORTOOBIG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VERSIONMISMATCH)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVERSIONMISMATCH>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VERTEXNOTINGRID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVERTEXNOTINGRID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VOXELGRIDTOOBIG)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceVOXELGRIDTOOBIG>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WIDTHTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWIDTHTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WINDOWEXCESS)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWINDOWEXCESS>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WINDOWSTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWINDOWSTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WINDOWTOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWINDOWTOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WORKSPACETOOSMALL)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWORKSPACETOOSMALL>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRITEERROR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRITEERROR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRITEFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRITEFAILED>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGARCHITECTURE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRONGARCHITECTURE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGCKTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRONGCKTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGCONIC)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRONGCONIC>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGDATATYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRONGDATATYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGSEGMENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRONGSEGMENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGSPKTYPE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceWRONGSPKTYPE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(YEAROUTOFRANGE)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceYEAROUTOFRANGE>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROBORESIGHT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROBORESIGHT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROBOUNDSEXTENT)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROBOUNDSEXTENT>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROFRAMEID)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROFRAMEID>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROLENGTHCOLUMN)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROLENGTHCOLUMN>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROPOSITION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROPOSITION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROQUATERNION)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROQUATERNION>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROSTEP)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROSTEP>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROVECTOR)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROVECTOR>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROVELOCITY)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZEROVELOCITY>(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZZHOLDDGETFAILED)",{[](const std::string &shortMessage, const std::string &explanation, const std::string &longMessage, const std::string &traceback) {return std::make_unique<tudat::exceptions::SpiceZZHOLDDGETFAILED>(shortMessage,explanation, longMessage,traceback);}}}
    };
    // clang-format on

    auto it = spiceExceptions.find(shortMessage);
    if( it != spiceExceptions.end( ) )
    {
        return it->second(shortMessage, explanation, longMessage, traceback);
    }
    else
    {
        return std::make_unique< tudat::exceptions::SpiceError >( shortMessage, explanation, longMessage, traceback );
    }
}

}  // namespace spice_interface
}  // namespace tudat
