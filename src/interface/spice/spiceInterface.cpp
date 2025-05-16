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
    // error message lengths are defined in cspice/SpiceErr.h
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

    throwSpiceException( shortMessage, explanation, longMessage, traceback );

}

void throwSpiceException( const std::string &shortMessage,
                                   const std::string &explanation,
                                   const std::string &longMessage,
                                   const std::string &traceback )
{

    using namespace tudat::exceptions;
    /* 
    Assuming you have the cspice source code in a directory called cspice,
    SPICE errors can be retrieved using the following commands:
    cd cspice
    grep -rhPo 'sigerr_\(\K"SPICE\([A-Za-z]+\)"' | sort | uniq
    */
   // clang-format off
    std::map< std::string, std::function<void()> > spiceExceptions = {
    {"SPICE(ADDRESSOUTOFBOUNDS)",{[=]() {throw SpiceADDRESSOUTOFBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AGENTLISTOVERFLOW)",{[=]() {throw SpiceAGENTLISTOVERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ALLGONE)",{[=]() {throw SpiceALLGONE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AMBIGTEMPL)",{[=]() {throw SpiceAMBIGTEMPL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ARRAYSIZEMISMATCH)",{[=]() {throw SpiceARRAYSIZEMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ARRAYTOOSMALL)",{[=]() {throw SpiceARRAYTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AVALOUTOFRANGE)",{[=]() {throw SpiceAVALOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(AXISUNDERFLOW)",{[=]() {throw SpiceAXISUNDERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADACTION)",{[=]() {throw SpiceBADACTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADADDRESS)",{[=]() {throw SpiceBADADDRESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGLE)",{[=]() {throw SpiceBADANGLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGLEUNITS)",{[=]() {throw SpiceBADANGLEUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGRATEERROR)",{[=]() {throw SpiceBADANGRATEERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGULARRATE)",{[=]() {throw SpiceBADANGULARRATE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADANGULARRATEFLAG)",{[=]() {throw SpiceBADANGULARRATEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADARCHITECTURE)",{[=]() {throw SpiceBADARCHITECTURE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADARRAYSIZE)",{[=]() {throw SpiceBADARRAYSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADATTIME)",{[=]() {throw SpiceBADATTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADATTRIBUTE)",{[=]() {throw SpiceBADATTRIBUTE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADATTRIBUTES)",{[=]() {throw SpiceBADATTRIBUTES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAUVALUE)",{[=]() {throw SpiceBADAUVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAVFLAG)",{[=]() {throw SpiceBADAVFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAVFRAMEFLAG)",{[=]() {throw SpiceBADAVFRAMEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAXIS)",{[=]() {throw SpiceBADAXIS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAXISLENGTH)",{[=]() {throw SpiceBADAXISLENGTH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADAXISNUMBERS)",{[=]() {throw SpiceBADAXISNUMBERS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBLOCKSIZE)",{[=]() {throw SpiceBADBLOCKSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBODYID)",{[=]() {throw SpiceBADBODYID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBORESIGHTSPEC)",{[=]() {throw SpiceBADBORESIGHTSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADBOUNDARY)",{[=]() {throw SpiceBADBOUNDARY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCATALOGFILE)",{[=]() {throw SpiceBADCATALOGFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCENTERNAME)",{[=]() {throw SpiceBADCENTERNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCHECKFLAG)",{[=]() {throw SpiceBADCHECKFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCKTYPESPEC)",{[=]() {throw SpiceBADCKTYPESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOARSEVOXSCALE)",{[=]() {throw SpiceBADCOARSEVOXSCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOLUMDECL)",{[=]() {throw SpiceBADCOLUMDECL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOLUMNCOUNT)",{[=]() {throw SpiceBADCOLUMNCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOLUMNDECL)",{[=]() {throw SpiceBADCOLUMNDECL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOMMENTAREA)",{[=]() {throw SpiceBADCOMMENTAREA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOMPNUMBER)",{[=]() {throw SpiceBADCOMPNUMBER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOORDBOUNDS)",{[=]() {throw SpiceBADCOORDBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCOORDSYS)",{[=]() {throw SpiceBADCOORDSYS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADCURVETYPE)",{[=]() {throw SpiceBADCURVETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDAFTRANSFERFILE)",{[=]() {throw SpiceBADDAFTRANSFERFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASCOMMENTAREA)",{[=]() {throw SpiceBADDASCOMMENTAREA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASDIRECTORY)",{[=]() {throw SpiceBADDASDIRECTORY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASFILE)",{[=]() {throw SpiceBADDASFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDASTRANSFERFILE)",{[=]() {throw SpiceBADDASTRANSFERFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATALINE)",{[=]() {throw SpiceBADDATALINE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATAORDERTOKEN)",{[=]() {throw SpiceBADDATAORDERTOKEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATATYPE)",{[=]() {throw SpiceBADDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDATATYPEFLAG)",{[=]() {throw SpiceBADDATATYPEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDEFAULTVALUE)",{[=]() {throw SpiceBADDEFAULTVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDESCRTIMES)",{[=]() {throw SpiceBADDESCRTIMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDIMENSION)",{[=]() {throw SpiceBADDIMENSION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDIMENSIONS)",{[=]() {throw SpiceBADDIMENSIONS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDIRECTION)",{[=]() {throw SpiceBADDIRECTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDOUBLEPRECISION)",{[=]() {throw SpiceBADDOUBLEPRECISION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADDOWNSAMPLINGTOL)",{[=]() {throw SpiceBADDOWNSAMPLINGTOL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADECCENTRICITY)",{[=]() {throw SpiceBADECCENTRICITY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADENDPOINTS)",{[=]() {throw SpiceBADENDPOINTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADEULERANGLEUNITS)",{[=]() {throw SpiceBADEULERANGLEUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFILEFORMAT)",{[=]() {throw SpiceBADFILEFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFILENAME)",{[=]() {throw SpiceBADFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFINEVOXELSCALE)",{[=]() {throw SpiceBADFINEVOXELSCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFORMATSPECIFIER)",{[=]() {throw SpiceBADFORMATSPECIFIER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAME)",{[=]() {throw SpiceBADFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAMECLASS)",{[=]() {throw SpiceBADFRAMECLASS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAMECOUNT)",{[=]() {throw SpiceBADFRAMECOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFRAMESPEC)",{[=]() {throw SpiceBADFRAMESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFROMTIME)",{[=]() {throw SpiceBADFROMTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFROMTIMESYSTEM)",{[=]() {throw SpiceBADFROMTIMESYSTEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADFROMTIMETYPE)",{[=]() {throw SpiceBADFROMTIMETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADGEOMETRY)",{[=]() {throw SpiceBADGEOMETRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADGM)",{[=]() {throw SpiceBADGM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADHARDSPACE)",{[=]() {throw SpiceBADHARDSPACE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADHERMITDEGREE)",{[=]() {throw SpiceBADHERMITDEGREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINDEX)",{[=]() {throw SpiceBADINDEX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINITSTATE)",{[=]() {throw SpiceBADINITSTATE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTDATALINE)",{[=]() {throw SpiceBADINPUTDATALINE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTETTIME)",{[=]() {throw SpiceBADINPUTETTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTTYPE)",{[=]() {throw SpiceBADINPUTTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINPUTUTCTIME)",{[=]() {throw SpiceBADINPUTUTCTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINSTRUMENTID)",{[=]() {throw SpiceBADINSTRUMENTID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADINTEGER)",{[=]() {throw SpiceBADINTEGER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADKERNELTYPE)",{[=]() {throw SpiceBADKERNELTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADKERNELVARTYPE)",{[=]() {throw SpiceBADKERNELVARTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLAGRANGEDEGREE)",{[=]() {throw SpiceBADLAGRANGEDEGREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLATITUDEBOUNDS)",{[=]() {throw SpiceBADLATITUDEBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLATITUDERANGE)",{[=]() {throw SpiceBADLATITUDERANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLATUSRECTUM)",{[=]() {throw SpiceBADLATUSRECTUM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLEAPSECONDS)",{[=]() {throw SpiceBADLEAPSECONDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLIMBLOCUSMIX)",{[=]() {throw SpiceBADLIMBLOCUSMIX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLINEPERRECCOUNT)",{[=]() {throw SpiceBADLINEPERRECCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLISTFILENAME)",{[=]() {throw SpiceBADLISTFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADLONGITUDERANGE)",{[=]() {throw SpiceBADLONGITUDERANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMATRIX)",{[=]() {throw SpiceBADMATRIX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMEANMOTION)",{[=]() {throw SpiceBADMEANMOTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMECCENTRICITY)",{[=]() {throw SpiceBADMECCENTRICITY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMETHODSYNTAX)",{[=]() {throw SpiceBADMETHODSYNTAX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMIDNIGHTTYPE)",{[=]() {throw SpiceBADMIDNIGHTTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMSEMIMAJOR)",{[=]() {throw SpiceBADMSEMIMAJOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADMSOPQUATERNION)",{[=]() {throw SpiceBADMSOPQUATERNION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADNOFDIGITS)",{[=]() {throw SpiceBADNOFDIGITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADNOFSTATES)",{[=]() {throw SpiceBADNOFSTATES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADNUMBEROFPOINTS)",{[=]() {throw SpiceBADNUMBEROFPOINTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOBJECTID)",{[=]() {throw SpiceBADOBJECTID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOBJECTNAME)",{[=]() {throw SpiceBADOBJECTNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETANGLES)",{[=]() {throw SpiceBADOFFSETANGLES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETANGUNITS)",{[=]() {throw SpiceBADOFFSETANGUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETAXESFORMAT)",{[=]() {throw SpiceBADOFFSETAXESFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOFFSETAXISXYZ)",{[=]() {throw SpiceBADOFFSETAXISXYZ(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADORBITALPERIOD)",{[=]() {throw SpiceBADORBITALPERIOD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOUTPUTSPKTYPE)",{[=]() {throw SpiceBADOUTPUTSPKTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADOUTPUTTYPE)",{[=]() {throw SpiceBADOUTPUTTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPARTNUMBER)",{[=]() {throw SpiceBADPARTNUMBER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPCKVALUE)",{[=]() {throw SpiceBADPCKVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPECCENTRICITY)",{[=]() {throw SpiceBADPECCENTRICITY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPERIAPSEVALUE)",{[=]() {throw SpiceBADPERIAPSEVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPICTURE)",{[=]() {throw SpiceBADPICTURE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPLATECOUNT)",{[=]() {throw SpiceBADPLATECOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPODLOCATION)",{[=]() {throw SpiceBADPODLOCATION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPRECVALUE)",{[=]() {throw SpiceBADPRECVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADPRIORITYSPEC)",{[=]() {throw SpiceBADPRIORITYSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADQUATSIGN)",{[=]() {throw SpiceBADQUATSIGN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADQUATTHRESHOLD)",{[=]() {throw SpiceBADQUATTHRESHOLD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRADIUS)",{[=]() {throw SpiceBADRADIUS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRADIUSCOUNT)",{[=]() {throw SpiceBADRADIUSCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRATEFRAMEFLAG)",{[=]() {throw SpiceBADRATEFRAMEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRATETHRESHOLD)",{[=]() {throw SpiceBADRATETHRESHOLD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADRECORDCOUNT)",{[=]() {throw SpiceBADRECORDCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADREFVECTORSPEC)",{[=]() {throw SpiceBADREFVECTORSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTATIONAXISXYZ)",{[=]() {throw SpiceBADROTATIONAXISXYZ(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTATIONSORDER)",{[=]() {throw SpiceBADROTATIONSORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTATIONTYPE)",{[=]() {throw SpiceBADROTATIONTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROTAXESFORMAT)",{[=]() {throw SpiceBADROTAXESFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADROWCOUNT)",{[=]() {throw SpiceBADROWCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSCID)",{[=]() {throw SpiceBADSCID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSEMIAXIS)",{[=]() {throw SpiceBADSEMIAXIS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSEMILATUS)",{[=]() {throw SpiceBADSEMILATUS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSHAPE)",{[=]() {throw SpiceBADSHAPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOLDAY)",{[=]() {throw SpiceBADSOLDAY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOLINDEX)",{[=]() {throw SpiceBADSOLINDEX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOLTIME)",{[=]() {throw SpiceBADSOLTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSOURCERADIUS)",{[=]() {throw SpiceBADSOURCERADIUS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSPICEQUATERNION)",{[=]() {throw SpiceBADSPICEQUATERNION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTARINDEX)",{[=]() {throw SpiceBADSTARINDEX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTARTTIME)",{[=]() {throw SpiceBADSTARTTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTDIONAME)",{[=]() {throw SpiceBADSTDIONAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSTOPTIME)",{[=]() {throw SpiceBADSTOPTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSUBSTR)",{[=]() {throw SpiceBADSUBSTR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSUBSTRINGBOUNDS)",{[=]() {throw SpiceBADSUBSTRINGBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADSURFACEMAP)",{[=]() {throw SpiceBADSURFACEMAP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTABLEFLAG)",{[=]() {throw SpiceBADTABLEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTERMLOCUSMIX)",{[=]() {throw SpiceBADTERMLOCUSMIX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEBOUNDS)",{[=]() {throw SpiceBADTIMEBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMECASE)",{[=]() {throw SpiceBADTIMECASE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMECOUNT)",{[=]() {throw SpiceBADTIMECOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEFORMAT)",{[=]() {throw SpiceBADTIMEFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEITEM)",{[=]() {throw SpiceBADTIMEITEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMEOFFSET)",{[=]() {throw SpiceBADTIMEOFFSET(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMESPEC)",{[=]() {throw SpiceBADTIMESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMESTRING)",{[=]() {throw SpiceBADTIMESTRING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMETYPE)",{[=]() {throw SpiceBADTIMETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTIMETYPEFLAG)",{[=]() {throw SpiceBADTIMETYPEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTLE)",{[=]() {throw SpiceBADTLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTLECOVERAGEPAD)",{[=]() {throw SpiceBADTLECOVERAGEPAD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTLEPADS)",{[=]() {throw SpiceBADTLEPADS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTOTIME)",{[=]() {throw SpiceBADTOTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTOTIMESYSTEM)",{[=]() {throw SpiceBADTOTIMESYSTEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTOTIMETYPE)",{[=]() {throw SpiceBADTOTIMETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADTYPESHAPECOMBO)",{[=]() {throw SpiceBADTYPESHAPECOMBO(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARASSIGN)",{[=]() {throw SpiceBADVARASSIGN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARIABLESIZE)",{[=]() {throw SpiceBADVARIABLESIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARIABLETYPE)",{[=]() {throw SpiceBADVARIABLETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVARNAME)",{[=]() {throw SpiceBADVARNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVECTOR)",{[=]() {throw SpiceBADVECTOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVERTEXCOUNT)",{[=]() {throw SpiceBADVERTEXCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADVERTEXINDEX)",{[=]() {throw SpiceBADVERTEXINDEX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BADWINDOWSIZE)",{[=]() {throw SpiceBADWINDOWSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BARRAYTOOSMALL)",{[=]() {throw SpiceBARRAYTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BARYCENTEREPHEM)",{[=]() {throw SpiceBARYCENTEREPHEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BARYCENTERIDCODE)",{[=]() {throw SpiceBARYCENTERIDCODE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BEFOREBEGSTR)",{[=]() {throw SpiceBEFOREBEGSTR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKCOMMANDLINE)",{[=]() {throw SpiceBLANKCOMMANDLINE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKFILENAME)",{[=]() {throw SpiceBLANKFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKFILETYPE)",{[=]() {throw SpiceBLANKFILETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKINPUTFILENAME)",{[=]() {throw SpiceBLANKINPUTFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKINPUTTIME)",{[=]() {throw SpiceBLANKINPUTTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKNAMEASSIGNED)",{[=]() {throw SpiceBLANKNAMEASSIGNED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKOUTPTFILENAME)",{[=]() {throw SpiceBLANKOUTPTFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKSCLKSTRING)",{[=]() {throw SpiceBLANKSCLKSTRING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLANKTIMEFORMAT)",{[=]() {throw SpiceBLANKTIMEFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BLOCKSNOTEVEN)",{[=]() {throw SpiceBLOCKSNOTEVEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BODIESNOTDISTINCT)",{[=]() {throw SpiceBODIESNOTDISTINCT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BODYANDCENTERSAME)",{[=]() {throw SpiceBODYANDCENTERSAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOGUSENTRY)",{[=]() {throw SpiceBOGUSENTRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BORESIGHTMISSING)",{[=]() {throw SpiceBORESIGHTMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDARYMISSING)",{[=]() {throw SpiceBOUNDARYMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDARYTOOBIG)",{[=]() {throw SpiceBOUNDARYTOOBIG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDSDISAGREE)",{[=]() {throw SpiceBOUNDSDISAGREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BOUNDSOUTOFORDER)",{[=]() {throw SpiceBOUNDSOUTOFORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUFFEROVERFLOW)",{[=]() {throw SpiceBUFFEROVERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUFFERSIZESMISMATCH)",{[=]() {throw SpiceBUFFERSIZESMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUFFERTOOSMALL)",{[=]() {throw SpiceBUFFERTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUG)",{[=]() {throw SpiceBUG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(BUGWRITEFAILED)",{[=]() {throw SpiceBUGWRITEFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CALLCKBSSFIRST)",{[=]() {throw SpiceCALLCKBSSFIRST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CALLEDOUTOFORDER)",{[=]() {throw SpiceCALLEDOUTOFORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CALLZZDSKBSSFIRST)",{[=]() {throw SpiceCALLZZDSKBSSFIRST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTFINDGRP)",{[=]() {throw SpiceCANNOTFINDGRP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTGETPACKET)",{[=]() {throw SpiceCANNOTGETPACKET(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTMAKEFILE)",{[=]() {throw SpiceCANNOTMAKEFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANNOTPICKFRAME)",{[=]() {throw SpiceCANNOTPICKFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANTFINDFRAME)",{[=]() {throw SpiceCANTFINDFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANTGETROTATIONTYPE)",{[=]() {throw SpiceCANTGETROTATIONTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CANTUSEPERIAPEPOCH)",{[=]() {throw SpiceCANTUSEPERIAPEPOCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CBNOSUCHSTR)",{[=]() {throw SpiceCBNOSUCHSTR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CELLARRAYTOOSMALL)",{[=]() {throw SpiceCELLARRAYTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CELLTOOSMALL)",{[=]() {throw SpiceCELLTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKBOGUSENTRY)",{[=]() {throw SpiceCKBOGUSENTRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKDOESNTEXIST)",{[=]() {throw SpiceCKDOESNTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKFILE)",{[=]() {throw SpiceCKFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKNONEXISTREC)",{[=]() {throw SpiceCKNONEXISTREC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKTOOMANYFILES)",{[=]() {throw SpiceCKTOOMANYFILES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKUNKNOWNDATATYPE)",{[=]() {throw SpiceCKUNKNOWNDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CKWRONGDATATYPE)",{[=]() {throw SpiceCKWRONGDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CMDERROR)",{[=]() {throw SpiceCMDERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CMDPARSEERROR)",{[=]() {throw SpiceCMDPARSEERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COARSEGRIDOVERFLOW)",{[=]() {throw SpiceCOARSEGRIDOVERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COLDESCTABLEFULL)",{[=]() {throw SpiceCOLDESCTABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COLUMNTOOSMALL)",{[=]() {throw SpiceCOLUMNTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMMANDTOOLONG)",{[=]() {throw SpiceCOMMANDTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMMENTTOOLONG)",{[=]() {throw SpiceCOMMENTTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMMFILENOTEXIST)",{[=]() {throw SpiceCOMMFILENOTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMPETINGEPOCHSPEC)",{[=]() {throw SpiceCOMPETINGEPOCHSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COMPETINGFRAMESPEC)",{[=]() {throw SpiceCOMPETINGFRAMESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COORDSYSNOTREC)",{[=]() {throw SpiceCOORDSYSNOTREC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COUNTMISMATCH)",{[=]() {throw SpiceCOUNTMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COUNTTOOLARGE)",{[=]() {throw SpiceCOUNTTOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(COVERAGEGAP)",{[=]() {throw SpiceCOVERAGEGAP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(CROSSANGLEMISSING)",{[=]() {throw SpiceCROSSANGLEMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFBADCRECLEN)",{[=]() {throw SpiceDAFBADCRECLEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFBEGGTEND)",{[=]() {throw SpiceDAFBEGGTEND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFCRNOTFOUND)",{[=]() {throw SpiceDAFCRNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFDPWRITEFAIL)",{[=]() {throw SpiceDAFDPWRITEFAIL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFFRNOTFOUND)",{[=]() {throw SpiceDAFFRNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFFTFULL)",{[=]() {throw SpiceDAFFTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFILLEGWRITE)",{[=]() {throw SpiceDAFILLEGWRITE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFINVALIDACCESS)",{[=]() {throw SpiceDAFINVALIDACCESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFINVALIDPARAMS)",{[=]() {throw SpiceDAFINVALIDPARAMS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNEGADDR)",{[=]() {throw SpiceDAFNEGADDR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNEWCONFLICT)",{[=]() {throw SpiceDAFNEWCONFLICT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOIFNMATCH)",{[=]() {throw SpiceDAFNOIFNMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNONAMEMATCH)",{[=]() {throw SpiceDAFNONAMEMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNORESV)",{[=]() {throw SpiceDAFNORESV(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSEARCH)",{[=]() {throw SpiceDAFNOSEARCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHADDR)",{[=]() {throw SpiceDAFNOSUCHADDR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHFILE)",{[=]() {throw SpiceDAFNOSUCHFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHHANDLE)",{[=]() {throw SpiceDAFNOSUCHHANDLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOSUCHUNIT)",{[=]() {throw SpiceDAFNOSUCHUNIT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFNOWRITE)",{[=]() {throw SpiceDAFNOWRITE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFOVERFLOW)",{[=]() {throw SpiceDAFOVERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFREADFAIL)",{[=]() {throw SpiceDAFREADFAIL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DAFWRITEFAIL)",{[=]() {throw SpiceDAFWRITEFAIL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASFILEREADFAILED)",{[=]() {throw SpiceDASFILEREADFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASFILEWRITEFAILED)",{[=]() {throw SpiceDASFILEWRITEFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASFTFULL)",{[=]() {throw SpiceDASFTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASINVALIDACCESS)",{[=]() {throw SpiceDASINVALIDACCESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASINVALIDCOUNT)",{[=]() {throw SpiceDASINVALIDCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASINVALIDTYPE)",{[=]() {throw SpiceDASINVALIDTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHADDRESS)",{[=]() {throw SpiceDASNOSUCHADDRESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHFILE)",{[=]() {throw SpiceDASNOSUCHFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHHANDLE)",{[=]() {throw SpiceDASNOSUCHHANDLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOSUCHUNIT)",{[=]() {throw SpiceDASNOSUCHUNIT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASNOTEMPTY)",{[=]() {throw SpiceDASNOTEMPTY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASREADFAIL)",{[=]() {throw SpiceDASREADFAIL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DASWRITEFAIL)",{[=]() {throw SpiceDASWRITEFAIL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATAITEMLIMITEXCEEDED)",{[=]() {throw SpiceDATAITEMLIMITEXCEEDED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATAREADFAILED)",{[=]() {throw SpiceDATAREADFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATAWIDTHERROR)",{[=]() {throw SpiceDATAWIDTHERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DATEEXPECTED)",{[=]() {throw SpiceDATEEXPECTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DECODINGERROR)",{[=]() {throw SpiceDECODINGERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGENERATECASE)",{[=]() {throw SpiceDEGENERATECASE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGENERATEINTERVAL)",{[=]() {throw SpiceDEGENERATEINTERVAL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGENERATESURFACE)",{[=]() {throw SpiceDEGENERATESURFACE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEGREEOUTOFRANGE)",{[=]() {throw SpiceDEGREEOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEPENDENTVECTORS)",{[=]() {throw SpiceDEPENDENTVECTORS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DEVICENAMETOOLONG)",{[=]() {throw SpiceDEVICENAMETOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIFFLINETOOLARGE)",{[=]() {throw SpiceDIFFLINETOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIFFLINETOOSMALL)",{[=]() {throw SpiceDIFFLINETOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIMENSIONTOOSMALL)",{[=]() {throw SpiceDIMENSIONTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DISARRAY)",{[=]() {throw SpiceDISARRAY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DISORDER)",{[=]() {throw SpiceDISORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DIVIDEBYZERO)",{[=]() {throw SpiceDIVIDEBYZERO(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DSKBOGUSENTRY)",{[=]() {throw SpiceDSKBOGUSENTRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DSKDATANOTFOUND)",{[=]() {throw SpiceDSKDATANOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DSKTOOMANYFILES)",{[=]() {throw SpiceDSKTOOMANYFILES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DTOUTOFRANGE)",{[=]() {throw SpiceDTOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DUBIOUSMETHOD)",{[=]() {throw SpiceDUBIOUSMETHOD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(DUPLICATETIMES)",{[=]() {throw SpiceDUPLICATETIMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ECCOUTOFBOUNDS)",{[=]() {throw SpiceECCOUTOFBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ECCOUTOFRANGE)",{[=]() {throw SpiceECCOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKCOLATTRTABLEFULL)",{[=]() {throw SpiceEKCOLATTRTABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKCOLNUMMISMATCH)",{[=]() {throw SpiceEKCOLNUMMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKFILE)",{[=]() {throw SpiceEKFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKFILETABLEFULL)",{[=]() {throw SpiceEKFILETABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKMISSINGCOLUMN)",{[=]() {throw SpiceEKMISSINGCOLUMN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKNOSEGMENTS)",{[=]() {throw SpiceEKNOSEGMENTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKSEGTABLEFULL)",{[=]() {throw SpiceEKSEGTABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EKTABLELISTFULL)",{[=]() {throw SpiceEKTABLELISTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ELEMENTSTOOSHORT)",{[=]() {throw SpiceELEMENTSTOOSHORT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EMPTYINPUTFILE)",{[=]() {throw SpiceEMPTYINPUTFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EMPTYSEGMENT)",{[=]() {throw SpiceEMPTYSEGMENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ENDOFFILE)",{[=]() {throw SpiceENDOFFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ENDPOINTSMATCH)",{[=]() {throw SpiceENDPOINTSMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ERROREXIT)",{[=]() {throw SpiceERROREXIT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EVECOUTOFRANGE)",{[=]() {throw SpiceEVECOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EVENHERMITDEGREE)",{[=]() {throw SpiceEVENHERMITDEGREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EVILBOGUSENTRY)",{[=]() {throw SpiceEVILBOGUSENTRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(EXTERNALOPEN)",{[=]() {throw SpiceEXTERNALOPEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FACENOTFOUND)",{[=]() {throw SpiceFACENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FAKESCLKEXISTS)",{[=]() {throw SpiceFAKESCLKEXISTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILARCHMISMATCH)",{[=]() {throw SpiceFILARCHMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILARCMISMATCH)",{[=]() {throw SpiceFILARCMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEALREADYEXISTS)",{[=]() {throw SpiceFILEALREADYEXISTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILECURRENTLYOPEN)",{[=]() {throw SpiceFILECURRENTLYOPEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEDELETEFAILED)",{[=]() {throw SpiceFILEDELETEFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEDOESNOTEXIST)",{[=]() {throw SpiceFILEDOESNOTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEEXISTS)",{[=]() {throw SpiceFILEEXISTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEISNOTSPK)",{[=]() {throw SpiceFILEISNOTSPK(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENAMETOOLONG)",{[=]() {throw SpiceFILENAMETOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENOTCONNECTED)",{[=]() {throw SpiceFILENOTCONNECTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENOTFOUND)",{[=]() {throw SpiceFILENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILENOTOPEN)",{[=]() {throw SpiceFILENOTOPEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENCONFLICT)",{[=]() {throw SpiceFILEOPENCONFLICT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENERROR)",{[=]() {throw SpiceFILEOPENERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENFAIL)",{[=]() {throw SpiceFILEOPENFAIL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEOPENFAILED)",{[=]() {throw SpiceFILEOPENFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEREADERROR)",{[=]() {throw SpiceFILEREADERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEREADFAILED)",{[=]() {throw SpiceFILEREADFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILETABLEFULL)",{[=]() {throw SpiceFILETABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILETRUNCATED)",{[=]() {throw SpiceFILETRUNCATED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FILEWRITEFAILED)",{[=]() {throw SpiceFILEWRITEFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FIRSTRECORDMISMATCH)",{[=]() {throw SpiceFIRSTRECORDMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FKDOESNTEXIST)",{[=]() {throw SpiceFKDOESNTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FMTITEMLIMITEXCEEDED)",{[=]() {throw SpiceFMTITEMLIMITEXCEEDED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATDATAMISMATCH)",{[=]() {throw SpiceFORMATDATAMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATDOESNTAPPLY)",{[=]() {throw SpiceFORMATDOESNTAPPLY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATERROR)",{[=]() {throw SpiceFORMATERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATNOTAPPLICABLE)",{[=]() {throw SpiceFORMATNOTAPPLICABLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FORMATSTRINGTOOLONG)",{[=]() {throw SpiceFORMATSTRINGTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FOVTOOWIDE)",{[=]() {throw SpiceFOVTOOWIDE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEDATANOTFOUND)",{[=]() {throw SpiceFRAMEDATANOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEDEFERROR)",{[=]() {throw SpiceFRAMEDEFERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEIDNOTFOUND)",{[=]() {throw SpiceFRAMEIDNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEINFONOTFOUND)",{[=]() {throw SpiceFRAMEINFONOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMEMISSING)",{[=]() {throw SpiceFRAMEMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMENAMENOTFOUND)",{[=]() {throw SpiceFRAMENAMENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMENOTFOUND)",{[=]() {throw SpiceFRAMENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FRAMENOTRECOGNIZED)",{[=]() {throw SpiceFRAMENOTRECOGNIZED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FTFULL)",{[=]() {throw SpiceFTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(FTPXFERERROR)",{[=]() {throw SpiceFTPXFERERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(GRIDTOOLARGE)",{[=]() {throw SpiceGRIDTOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(HANDLENOTFOUND)",{[=]() {throw SpiceHANDLENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(HASHISFULL)",{[=]() {throw SpiceHASHISFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(HLULOCKFAILED)",{[=]() {throw SpiceHLULOCKFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IDCODENOTFOUND)",{[=]() {throw SpiceIDCODENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IDSTRINGTOOLONG)",{[=]() {throw SpiceIDSTRINGTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGALCHARACTER)",{[=]() {throw SpiceILLEGALCHARACTER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGALOPTIONNAME)",{[=]() {throw SpiceILLEGALOPTIONNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGSHIFTDIR)",{[=]() {throw SpiceILLEGSHIFTDIR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ILLEGTEMPL)",{[=]() {throw SpiceILLEGTEMPL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IMMUTABLEVALUE)",{[=]() {throw SpiceIMMUTABLEVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IMPROPERFILE)",{[=]() {throw SpiceIMPROPERFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IMPROPEROPEN)",{[=]() {throw SpiceIMPROPEROPEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INACTIVEOBJECT)",{[=]() {throw SpiceINACTIVEOBJECT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLEEOL)",{[=]() {throw SpiceINCOMPATIBLEEOL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLENUMREF)",{[=]() {throw SpiceINCOMPATIBLENUMREF(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLESCALE)",{[=]() {throw SpiceINCOMPATIBLESCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPATIBLEUNITS)",{[=]() {throw SpiceINCOMPATIBLEUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCOMPLETEFRAME)",{[=]() {throw SpiceINCOMPLETEFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTCENTERID)",{[=]() {throw SpiceINCONSISTCENTERID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTENTTIMES)",{[=]() {throw SpiceINCONSISTENTTIMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTFRAME)",{[=]() {throw SpiceINCONSISTFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTSTARTTIME)",{[=]() {throw SpiceINCONSISTSTARTTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCONSISTSTOPTIME)",{[=]() {throw SpiceINCONSISTSTOPTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INCORRECTUSAGE)",{[=]() {throw SpiceINCORRECTUSAGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDEFINITELOCALSECOND)",{[=]() {throw SpiceINDEFINITELOCALSECOND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDEXOUTOFRANGE)",{[=]() {throw SpiceINDEXOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDEXTOOLARGE)",{[=]() {throw SpiceINDEXTOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INDICESOUTOFORDER)",{[=]() {throw SpiceINDICESOUTOFORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTDOESNOTEXIST)",{[=]() {throw SpiceINPUTDOESNOTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTFILENOTEXIST)",{[=]() {throw SpiceINPUTFILENOTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTOUTOFBOUNDS)",{[=]() {throw SpiceINPUTOUTOFBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INPUTSTOOLARGE)",{[=]() {throw SpiceINPUTSTOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INQUIREERROR)",{[=]() {throw SpiceINQUIREERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INQUIREFAILED)",{[=]() {throw SpiceINQUIREFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSIDEBODY)",{[=]() {throw SpiceINSIDEBODY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFFICIENTANGLES)",{[=]() {throw SpiceINSUFFICIENTANGLES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFFICIENTDATA)",{[=]() {throw SpiceINSUFFICIENTDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFFLEN)",{[=]() {throw SpiceINSUFFLEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INSUFPTRSIZE)",{[=]() {throw SpiceINSUFPTRSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTERVALSTARTNOTFOUND)",{[=]() {throw SpiceINTERVALSTARTNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTINDEXTOOSMALL)",{[=]() {throw SpiceINTINDEXTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTLENNOTPOS)",{[=]() {throw SpiceINTLENNOTPOS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INTOUTOFRANGE)",{[=]() {throw SpiceINTOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDACCESS)",{[=]() {throw SpiceINVALIDACCESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDACTION)",{[=]() {throw SpiceINVALIDACTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDADD)",{[=]() {throw SpiceINVALIDADD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDADDRESS)",{[=]() {throw SpiceINVALIDADDRESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDANGLE)",{[=]() {throw SpiceINVALIDANGLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDARCHTYPE)",{[=]() {throw SpiceINVALIDARCHTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDARGUMENT)",{[=]() {throw SpiceINVALIDARGUMENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDAXIS)",{[=]() {throw SpiceINVALIDAXIS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDAXISLENGTH)",{[=]() {throw SpiceINVALIDAXISLENGTH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDBOUNDS)",{[=]() {throw SpiceINVALIDBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCARDINALITY)",{[=]() {throw SpiceINVALIDCARDINALITY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCASE)",{[=]() {throw SpiceINVALIDCASE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCOLUMN)",{[=]() {throw SpiceINVALIDCOLUMN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCONSTSTEP)",{[=]() {throw SpiceINVALIDCONSTSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDCOUNT)",{[=]() {throw SpiceINVALIDCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDATA)",{[=]() {throw SpiceINVALIDDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDATACOUNT)",{[=]() {throw SpiceINVALIDDATACOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDATATYPE)",{[=]() {throw SpiceINVALIDDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDEGREE)",{[=]() {throw SpiceINVALIDDEGREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDESCRTIME)",{[=]() {throw SpiceINVALIDDESCRTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDIMENSION)",{[=]() {throw SpiceINVALIDDIMENSION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDIRECTION)",{[=]() {throw SpiceINVALIDDIRECTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDDIVISOR)",{[=]() {throw SpiceINVALIDDIVISOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDELLIPSE)",{[=]() {throw SpiceINVALIDELLIPSE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDENDPNTSPEC)",{[=]() {throw SpiceINVALIDENDPNTSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDENDPTS)",{[=]() {throw SpiceINVALIDENDPTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDEPOCH)",{[=]() {throw SpiceINVALIDEPOCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFILETYPE)",{[=]() {throw SpiceINVALIDFILETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFIXREF)",{[=]() {throw SpiceINVALIDFIXREF(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFLAG)",{[=]() {throw SpiceINVALIDFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFORMAT)",{[=]() {throw SpiceINVALIDFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFOV)",{[=]() {throw SpiceINVALIDFOV(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFRAME)",{[=]() {throw SpiceINVALIDFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDFRAMEDEF)",{[=]() {throw SpiceINVALIDFRAMEDEF(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDGEOMETRY)",{[=]() {throw SpiceINVALIDGEOMETRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDHANDLE)",{[=]() {throw SpiceINVALIDHANDLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDINDEX)",{[=]() {throw SpiceINVALIDINDEX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDINTEGER)",{[=]() {throw SpiceINVALIDINTEGER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLIMBTYPE)",{[=]() {throw SpiceINVALIDLIMBTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLISTITEM)",{[=]() {throw SpiceINVALIDLISTITEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLOCUS)",{[=]() {throw SpiceINVALIDLOCUS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDLONEXTENT)",{[=]() {throw SpiceINVALIDLONEXTENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDMETADATA)",{[=]() {throw SpiceINVALIDMETADATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDMETHOD)",{[=]() {throw SpiceINVALIDMETHOD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDMSGTYPE)",{[=]() {throw SpiceINVALIDMSGTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNAME)",{[=]() {throw SpiceINVALIDNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNODE)",{[=]() {throw SpiceINVALIDNODE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMBEROFINTERVALS)",{[=]() {throw SpiceINVALIDNUMBEROFINTERVALS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMBEROFRECORDS)",{[=]() {throw SpiceINVALIDNUMBEROFRECORDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMINT)",{[=]() {throw SpiceINVALIDNUMINT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMINTS)",{[=]() {throw SpiceINVALIDNUMINTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDNUMREC)",{[=]() {throw SpiceINVALIDNUMREC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDOCCTYPE)",{[=]() {throw SpiceINVALIDOCCTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDOPERATION)",{[=]() {throw SpiceINVALIDOPERATION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDOPTION)",{[=]() {throw SpiceINVALIDOPTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDPLANE)",{[=]() {throw SpiceINVALIDPLANE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDRADII)",{[=]() {throw SpiceINVALIDRADII(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDRADIUS)",{[=]() {throw SpiceINVALIDRADIUS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDREFFRAME)",{[=]() {throw SpiceINVALIDREFFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDREFVAL)",{[=]() {throw SpiceINVALIDREFVAL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDROLLSTEP)",{[=]() {throw SpiceINVALIDROLLSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCALE)",{[=]() {throw SpiceINVALIDSCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCLKRATE)",{[=]() {throw SpiceINVALIDSCLKRATE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCLKSTRING)",{[=]() {throw SpiceINVALIDSCLKSTRING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSCLKTIME)",{[=]() {throw SpiceINVALIDSCLKTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSEARCHSTEP)",{[=]() {throw SpiceINVALIDSEARCHSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSELECTION)",{[=]() {throw SpiceINVALIDSELECTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSHADOW)",{[=]() {throw SpiceINVALIDSHADOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSHAPE)",{[=]() {throw SpiceINVALIDSHAPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSHAPECOMBO)",{[=]() {throw SpiceINVALIDSHAPECOMBO(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSIZE)",{[=]() {throw SpiceINVALIDSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTARTTIME)",{[=]() {throw SpiceINVALIDSTARTTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTATE)",{[=]() {throw SpiceINVALIDSTATE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTEP)",{[=]() {throw SpiceINVALIDSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSTEPSIZE)",{[=]() {throw SpiceINVALIDSTEPSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSUBLIST)",{[=]() {throw SpiceINVALIDSUBLIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDSUBTYPE)",{[=]() {throw SpiceINVALIDSUBTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTABLENAME)",{[=]() {throw SpiceINVALIDTABLENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTABLESIZE)",{[=]() {throw SpiceINVALIDTABLESIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTARGET)",{[=]() {throw SpiceINVALIDTARGET(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTERMTYPE)",{[=]() {throw SpiceINVALIDTERMTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTEXT)",{[=]() {throw SpiceINVALIDTEXT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTIMEFORMAT)",{[=]() {throw SpiceINVALIDTIMEFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTIMESTRING)",{[=]() {throw SpiceINVALIDTIMESTRING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTLEORDER)",{[=]() {throw SpiceINVALIDTLEORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTOLERANCE)",{[=]() {throw SpiceINVALIDTOLERANCE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDTYPE)",{[=]() {throw SpiceINVALIDTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDVALUE)",{[=]() {throw SpiceINVALIDVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVALIDVERTEX)",{[=]() {throw SpiceINVALIDVERTEX(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(INVERSTARTSTOPTIME)",{[=]() {throw SpiceINVERSTARTSTOPTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(IRFNOTREC)",{[=]() {throw SpiceIRFNOTREC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ITEMNOTFOUND)",{[=]() {throw SpiceITEMNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ITEMNOTRECOGNIZED)",{[=]() {throw SpiceITEMNOTRECOGNIZED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ITERATIONEXCEEDED)",{[=]() {throw SpiceITERATIONEXCEEDED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERNELNOTLOADED)",{[=]() {throw SpiceKERNELNOTLOADED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERNELPOOLFULL)",{[=]() {throw SpiceKERNELPOOLFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERNELVARNOTFOUND)",{[=]() {throw SpiceKERNELVARNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERVARSETOVERFLOW)",{[=]() {throw SpiceKERVARSETOVERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KERVARTOOBIG)",{[=]() {throw SpiceKERVARTOOBIG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(KEYWORDNOTFOUND)",{[=]() {throw SpiceKEYWORDNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBCORRUPTED)",{[=]() {throw SpiceLBCORRUPTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBLINETOOLONG)",{[=]() {throw SpiceLBLINETOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBNOSUCHLINE)",{[=]() {throw SpiceLBNOSUCHLINE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LBTOOMANYLINES)",{[=]() {throw SpiceLBTOOMANYLINES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LOWERBOUNDTOOLOW)",{[=]() {throw SpiceLOWERBOUNDTOOLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(LSKDOESNTEXIST)",{[=]() {throw SpiceLSKDOESNTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MALFORMEDSEGMENT)",{[=]() {throw SpiceMALFORMEDSEGMENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MARKERNOTFOUND)",{[=]() {throw SpiceMARKERNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MESSAGETOOLONG)",{[=]() {throw SpiceMESSAGETOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISMATCHFROMTIMETYPE)",{[=]() {throw SpiceMISMATCHFROMTIMETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISMATCHOUTPUTFORMAT)",{[=]() {throw SpiceMISMATCHOUTPUTFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISMATCHTOTIMETYPE)",{[=]() {throw SpiceMISMATCHTOTIMETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGARGUMENTS)",{[=]() {throw SpiceMISSINGARGUMENTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCENTER)",{[=]() {throw SpiceMISSINGCENTER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCOLSTEP)",{[=]() {throw SpiceMISSINGCOLSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCOORDBOUND)",{[=]() {throw SpiceMISSINGCOORDBOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGCOORDSYS)",{[=]() {throw SpiceMISSINGCOORDSYS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATA)",{[=]() {throw SpiceMISSINGDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATACLASS)",{[=]() {throw SpiceMISSINGDATACLASS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATAORDERTK)",{[=]() {throw SpiceMISSINGDATAORDERTK(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGDATATYPE)",{[=]() {throw SpiceMISSINGDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGEOT)",{[=]() {throw SpiceMISSINGEOT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGEPOCHTOKEN)",{[=]() {throw SpiceMISSINGEPOCHTOKEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGFRAME)",{[=]() {throw SpiceMISSINGFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGFRAMEVAR)",{[=]() {throw SpiceMISSINGFRAMEVAR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGGEOCONSTS)",{[=]() {throw SpiceMISSINGGEOCONSTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGHEIGHTREF)",{[=]() {throw SpiceMISSINGHEIGHTREF(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGHSCALE)",{[=]() {throw SpiceMISSINGHSCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGKPV)",{[=]() {throw SpiceMISSINGKPV(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGLEFTCOR)",{[=]() {throw SpiceMISSINGLEFTCOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGLEFTRTFLAG)",{[=]() {throw SpiceMISSINGLEFTRTFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGNCAPFLAG)",{[=]() {throw SpiceMISSINGNCAPFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGNCOLS)",{[=]() {throw SpiceMISSINGNCOLS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGNROWS)",{[=]() {throw SpiceMISSINGNROWS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGPLATETYPE)",{[=]() {throw SpiceMISSINGPLATETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGROWMAJFLAG)",{[=]() {throw SpiceMISSINGROWMAJFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGROWSTEP)",{[=]() {throw SpiceMISSINGROWSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGSCAPFLAG)",{[=]() {throw SpiceMISSINGSCAPFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGSURFACE)",{[=]() {throw SpiceMISSINGSURFACE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTIMEINFO)",{[=]() {throw SpiceMISSINGTIMEINFO(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTLEKEYWORD)",{[=]() {throw SpiceMISSINGTLEKEYWORD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTOPCOR)",{[=]() {throw SpiceMISSINGTOPCOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGTOPDOWNFLAG)",{[=]() {throw SpiceMISSINGTOPDOWNFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGVALUE)",{[=]() {throw SpiceMISSINGVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGVOXELSCALE)",{[=]() {throw SpiceMISSINGVOXELSCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(MISSINGWRAPFLAG)",{[=]() {throw SpiceMISSINGWRAPFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NAMENOTUNIQUE)",{[=]() {throw SpiceNAMENOTUNIQUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NAMESNOTRESOLVED)",{[=]() {throw SpiceNAMESNOTRESOLVED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NAMETABLEFULL)",{[=]() {throw SpiceNAMETABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NARATESFLAG)",{[=]() {throw SpiceNARATESFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NEGATIVETOL)",{[=]() {throw SpiceNEGATIVETOL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOACCEPTABLEDATA)",{[=]() {throw SpiceNOACCEPTABLEDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOANGULARRATEFLAG)",{[=]() {throw SpiceNOANGULARRATEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOARRAYSTARTED)",{[=]() {throw SpiceNOARRAYSTARTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOATTIME)",{[=]() {throw SpiceNOATTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOAVDATA)",{[=]() {throw SpiceNOAVDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOBODYID)",{[=]() {throw SpiceNOBODYID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCANDOSPKSPCKS)",{[=]() {throw SpiceNOCANDOSPKSPCKS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCENTERIDORNAME)",{[=]() {throw SpiceNOCENTERIDORNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCKSEGMENTTYPE)",{[=]() {throw SpiceNOCKSEGMENTTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCLASS)",{[=]() {throw SpiceNOCLASS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCOMMENTSFILE)",{[=]() {throw SpiceNOCOMMENTSFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCONVERG)",{[=]() {throw SpiceNOCONVERG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCONVERGENCE)",{[=]() {throw SpiceNOCONVERGENCE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOCURRENTARRAY)",{[=]() {throw SpiceNOCURRENTARRAY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODATAORDER)",{[=]() {throw SpiceNODATAORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODATATYPEFLAG)",{[=]() {throw SpiceNODATATYPEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODELIMCHARACTER)",{[=]() {throw SpiceNODELIMCHARACTER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODETOOFULL)",{[=]() {throw SpiceNODETOOFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODSKSEGMENT)",{[=]() {throw SpiceNODSKSEGMENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NODSKSEGMENTS)",{[=]() {throw SpiceNODSKSEGMENTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOENVVARIABLE)",{[=]() {throw SpiceNOENVVARIABLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOEULERANGLEUNITS)",{[=]() {throw SpiceNOEULERANGLEUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFILENAMES)",{[=]() {throw SpiceNOFILENAMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFILES)",{[=]() {throw SpiceNOFILES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFILESPEC)",{[=]() {throw SpiceNOFILESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAME)",{[=]() {throw SpiceNOFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMECONNECT)",{[=]() {throw SpiceNOFRAMECONNECT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMEDATA)",{[=]() {throw SpiceNOFRAMEDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMEINFO)",{[=]() {throw SpiceNOFRAMEINFO(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMENAME)",{[=]() {throw SpiceNOFRAMENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFRAMESKERNELNAME)",{[=]() {throw SpiceNOFRAMESKERNELNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFREELOGICALUNIT)",{[=]() {throw SpiceNOFREELOGICALUNIT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFREENODES)",{[=]() {throw SpiceNOFREENODES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFROMTIME)",{[=]() {throw SpiceNOFROMTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOFROMTIMESYSTEM)",{[=]() {throw SpiceNOFROMTIMESYSTEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOHEADNODE)",{[=]() {throw SpiceNOHEADNODE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINPUTDATATYPE)",{[=]() {throw SpiceNOINPUTDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINPUTFILENAME)",{[=]() {throw SpiceNOINPUTFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINSTRUMENTID)",{[=]() {throw SpiceNOINSTRUMENTID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOINTERVAL)",{[=]() {throw SpiceNOINTERVAL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOKERNELLOADED)",{[=]() {throw SpiceNOKERNELLOADED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLANDINGTIME)",{[=]() {throw SpiceNOLANDINGTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLEAPSECONDS)",{[=]() {throw SpiceNOLEAPSECONDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLINESPERRECCOUNT)",{[=]() {throw SpiceNOLINESPERRECCOUNT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLISTFILENAME)",{[=]() {throw SpiceNOLISTFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLOADEDDSKFILES)",{[=]() {throw SpiceNOLOADEDDSKFILES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLOADEDFILES)",{[=]() {throw SpiceNOLOADEDFILES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOLSKFILENAME)",{[=]() {throw SpiceNOLSKFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOMOREROOM)",{[=]() {throw SpiceNOMOREROOM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONCONICMOTION)",{[=]() {throw SpiceNONCONICMOTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONDISTINCTPAIR)",{[=]() {throw SpiceNONDISTINCTPAIR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONEMPTYENTRY)",{[=]() {throw SpiceNONEMPTYENTRY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONEMPTYTREE)",{[=]() {throw SpiceNONEMPTYTREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONEXISTELEMENTS)",{[=]() {throw SpiceNONEXISTELEMENTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONINTEGERFIELD)",{[=]() {throw SpiceNONINTEGERFIELD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONNUMERICSTRING)",{[=]() {throw SpiceNONNUMERICSTRING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSBUFLENGTH)",{[=]() {throw SpiceNONPOSBUFLENGTH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVEAXIS)",{[=]() {throw SpiceNONPOSITIVEAXIS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVEMASS)",{[=]() {throw SpiceNONPOSITIVEMASS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVERADIUS)",{[=]() {throw SpiceNONPOSITIVERADIUS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVESCALE)",{[=]() {throw SpiceNONPOSITIVESCALE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSITIVEVALUE)",{[=]() {throw SpiceNONPOSITIVEVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPOSPACKETSIZE)",{[=]() {throw SpiceNONPOSPACKETSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPRINTABLECHARS)",{[=]() {throw SpiceNONPRINTABLECHARS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPRINTINGCHAR)",{[=]() {throw SpiceNONPRINTINGCHAR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONPRINTINGCHARS)",{[=]() {throw SpiceNONPRINTINGCHARS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NONUNITNORMAL)",{[=]() {throw SpiceNONUNITNORMAL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOBJECTIDORNAME)",{[=]() {throw SpiceNOOBJECTIDORNAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOFFSETANGLEAXES)",{[=]() {throw SpiceNOOFFSETANGLEAXES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOFFSETANGLEUNITS)",{[=]() {throw SpiceNOOFFSETANGLEUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOUTPUTFILENAME)",{[=]() {throw SpiceNOOUTPUTFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOOUTPUTSPKTYPE)",{[=]() {throw SpiceNOOUTPUTSPKTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPARTITION)",{[=]() {throw SpiceNOPARTITION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPICTURE)",{[=]() {throw SpiceNOPICTURE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPOLYNOMIALDEGREE)",{[=]() {throw SpiceNOPOLYNOMIALDEGREE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPRECESSIONTYPE)",{[=]() {throw SpiceNOPRECESSIONTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOPRODUCERID)",{[=]() {throw SpiceNOPRODUCERID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOROTATIONORDER)",{[=]() {throw SpiceNOROTATIONORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSCID)",{[=]() {throw SpiceNOSCID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSCLKFILENAMES)",{[=]() {throw SpiceNOSCLKFILENAMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSECONDLINE)",{[=]() {throw SpiceNOSECONDLINE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSEGMENTSFOUND)",{[=]() {throw SpiceNOSEGMENTSFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSEPARATION)",{[=]() {throw SpiceNOSEPARATION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSLKFILENAME)",{[=]() {throw SpiceNOSLKFILENAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSOLMARKER)",{[=]() {throw SpiceNOSOLMARKER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSPACECRAFTID)",{[=]() {throw SpiceNOSPACECRAFTID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSTARTTIME)",{[=]() {throw SpiceNOSTARTTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSTOPTIME)",{[=]() {throw SpiceNOSTOPTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUCHFILE)",{[=]() {throw SpiceNOSUCHFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUCHHANDLE)",{[=]() {throw SpiceNOSUCHHANDLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUCHSYMBOL)",{[=]() {throw SpiceNOSUCHSYMBOL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOSUNGM)",{[=]() {throw SpiceNOSUNGM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTABINARYKERNEL)",{[=]() {throw SpiceNOTABINARYKERNEL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTACKFILE)",{[=]() {throw SpiceNOTACKFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTADAFFILE)",{[=]() {throw SpiceNOTADAFFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTADASFILE)",{[=]() {throw SpiceNOTADASFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTADPNUMBER)",{[=]() {throw SpiceNOTADPNUMBER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANDPNUMBER)",{[=]() {throw SpiceNOTANDPNUMBER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANINTEGER)",{[=]() {throw SpiceNOTANINTEGER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANINTEGERNUMBER)",{[=]() {throw SpiceNOTANINTEGERNUMBER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTANINTNUMBER)",{[=]() {throw SpiceNOTANINTNUMBER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTAPCKFILE)",{[=]() {throw SpiceNOTAPCKFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTAROTATION)",{[=]() {throw SpiceNOTAROTATION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTATEXTFILE)",{[=]() {throw SpiceNOTATEXTFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTATRANSFERFILE)",{[=]() {throw SpiceNOTATRANSFERFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTCOMPUTABLE)",{[=]() {throw SpiceNOTCOMPUTABLE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTDIMENSIONALLYEQUIV)",{[=]() {throw SpiceNOTDIMENSIONALLYEQUIV(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTDISJOINT)",{[=]() {throw SpiceNOTDISJOINT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTDISTINCT)",{[=]() {throw SpiceNOTDISTINCT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTENOUGHPEAS)",{[=]() {throw SpiceNOTENOUGHPEAS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTIMETYPEFLAG)",{[=]() {throw SpiceNOTIMETYPEFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTINDEXED)",{[=]() {throw SpiceNOTINDEXED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTINITIALIZED)",{[=]() {throw SpiceNOTINITIALIZED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTINPART)",{[=]() {throw SpiceNOTINPART(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTLEDATAFOROBJECT)",{[=]() {throw SpiceNOTLEDATAFOROBJECT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTLEGALCB)",{[=]() {throw SpiceNOTLEGALCB(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTOTIME)",{[=]() {throw SpiceNOTOTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTOTIMESYSTEM)",{[=]() {throw SpiceNOTOTIMESYSTEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTRANSLATION)",{[=]() {throw SpiceNOTRANSLATION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTRECOGNIZED)",{[=]() {throw SpiceNOTRECOGNIZED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTSEMCHECKED)",{[=]() {throw SpiceNOTSEMCHECKED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTSUPPORTED)",{[=]() {throw SpiceNOTSUPPORTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTTWOFIELDSCLK)",{[=]() {throw SpiceNOTTWOFIELDSCLK(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTTWOMODULI)",{[=]() {throw SpiceNOTTWOMODULI(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOTTWOOFFSETS)",{[=]() {throw SpiceNOTTWOOFFSETS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NOUNITSPEC)",{[=]() {throw SpiceNOUNITSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMBEREXPECTED)",{[=]() {throw SpiceNUMBEREXPECTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMCOEFFSNOTPOS)",{[=]() {throw SpiceNUMCOEFFSNOTPOS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMERICOVERFLOW)",{[=]() {throw SpiceNUMERICOVERFLOW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMPACKETSNOTPOS)",{[=]() {throw SpiceNUMPACKETSNOTPOS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMPARTSUNEQUAL)",{[=]() {throw SpiceNUMPARTSUNEQUAL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(NUMSTATESNOTPOS)",{[=]() {throw SpiceNUMSTATESNOTPOS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OBJECTLISTFULL)",{[=]() {throw SpiceOBJECTLISTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OBJECTSTOOCLOSE)",{[=]() {throw SpiceOBJECTSTOOCLOSE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ORBITDECAY)",{[=]() {throw SpiceORBITDECAY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTOFPLACEDELIMITER)",{[=]() {throw SpiceOUTOFPLACEDELIMITER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTOFRANGE)",{[=]() {throw SpiceOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTOFROOM)",{[=]() {throw SpiceOUTOFROOM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTFILEEXISTS)",{[=]() {throw SpiceOUTPUTFILEEXISTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTISNOTSPK)",{[=]() {throw SpiceOUTPUTISNOTSPK(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTTOOLONG)",{[=]() {throw SpiceOUTPUTTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(OUTPUTTOOSHORT)",{[=]() {throw SpiceOUTPUTTOOSHORT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PARSERNOTREADY)",{[=]() {throw SpicePARSERNOTREADY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PARTIALFRAMESPEC)",{[=]() {throw SpicePARTIALFRAMESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PASTENDSTR)",{[=]() {throw SpicePASTENDSTR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PATHMISMATCH)",{[=]() {throw SpicePATHMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PATHTOOLONG)",{[=]() {throw SpicePATHTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKDOESNTEXIST)",{[=]() {throw SpicePCKDOESNTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKFILE)",{[=]() {throw SpicePCKFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKFILETABLEFULL)",{[=]() {throw SpicePCKFILETABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PCKKRECTOOLARGE)",{[=]() {throw SpicePCKKRECTOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PLATELISTTOOSMALL)",{[=]() {throw SpicePLATELISTTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTEROUTOFRANGE)",{[=]() {throw SpicePOINTEROUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTERSETTOOBIG)",{[=]() {throw SpicePOINTERSETTOOBIG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTERTABLEFULL)",{[=]() {throw SpicePOINTERTABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTNOTFOUND)",{[=]() {throw SpicePOINTNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTNOTINSEGMENT)",{[=]() {throw SpicePOINTNOTINSEGMENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTNOTONSURFACE)",{[=]() {throw SpicePOINTNOTONSURFACE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTOFFSURFACE)",{[=]() {throw SpicePOINTOFFSURFACE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTONZAXIS)",{[=]() {throw SpicePOINTONZAXIS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(POINTTOOSMALL)",{[=]() {throw SpicePOINTTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(PTRARRAYTOOSMALL)",{[=]() {throw SpicePTRARRAYTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(QPARAMOUTOFRANGE)",{[=]() {throw SpiceQPARAMOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(QUERYFAILURE)",{[=]() {throw SpiceQUERYFAILURE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(QUERYNOTPARSED)",{[=]() {throw SpiceQUERYNOTPARSED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RADIIOUTOFORDER)",{[=]() {throw SpiceRADIIOUTOFORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RAYISZEROVECTOR)",{[=]() {throw SpiceRAYISZEROVECTOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(READFAILED)",{[=]() {throw SpiceREADFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RECORDNOTFOUND)",{[=]() {throw SpiceRECORDNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RECURSIONTOODEEP)",{[=]() {throw SpiceRECURSIONTOODEEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RECURSIVELOADING)",{[=]() {throw SpiceRECURSIVELOADING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REFANGLEMISSING)",{[=]() {throw SpiceREFANGLEMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REFVALNOTINTEGER)",{[=]() {throw SpiceREFVALNOTINTEGER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REFVECTORMISSING)",{[=]() {throw SpiceREFVECTORMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REQUESTOUTOFBOUNDS)",{[=]() {throw SpiceREQUESTOUTOFBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(REQUESTOUTOFORDER)",{[=]() {throw SpiceREQUESTOUTOFORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(RWCONFLICT)",{[=]() {throw SpiceRWCONFLICT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SBINSUFPTRSIZE)",{[=]() {throw SpiceSBINSUFPTRSIZE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SBTOOMANYSTRS)",{[=]() {throw SpiceSBTOOMANYSTRS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SCLKDOESNTEXIST)",{[=]() {throw SpiceSCLKDOESNTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SCLKTRUNCATED)",{[=]() {throw SpiceSCLKTRUNCATED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGIDTOOLONG)",{[=]() {throw SpiceSEGIDTOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGMENTNOTFOUND)",{[=]() {throw SpiceSEGMENTNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGMENTTABLEFULL)",{[=]() {throw SpiceSEGMENTTABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGTABLETOOSMALL)",{[=]() {throw SpiceSEGTABLETOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SEGTYPECONFLICT)",{[=]() {throw SpiceSEGTYPECONFLICT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SETEXCESS)",{[=]() {throw SpiceSETEXCESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SETTOOSMALL)",{[=]() {throw SpiceSETTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SETUPDOESNOTEXIST)",{[=]() {throw SpiceSETUPDOESNOTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SHAPEMISSING)",{[=]() {throw SpiceSHAPEMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SHAPENOTSUPPORTED)",{[=]() {throw SpiceSHAPENOTSUPPORTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SIZEMISMATCH)",{[=]() {throw SpiceSIZEMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SIZEOUTOFRANGE)",{[=]() {throw SpiceSIZEOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPACETOONARROW)",{[=]() {throw SpiceSPACETOONARROW(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPCRFLNOTCALLED)",{[=]() {throw SpiceSPCRFLNOTCALLED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPICEISTIRED)",{[=]() {throw SpiceSPICEISTIRED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKDOESNTEXIST)",{[=]() {throw SpiceSPKDOESNTEXIST(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKFILE)",{[=]() {throw SpiceSPKFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKFILETABLEFULL)",{[=]() {throw SpiceSPKFILETABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKINSUFFDATA)",{[=]() {throw SpiceSPKINSUFFDATA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKINVALIDOPTION)",{[=]() {throw SpiceSPKINVALIDOPTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKNOTASUBSET)",{[=]() {throw SpiceSPKNOTASUBSET(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKRECTOOLARGE)",{[=]() {throw SpiceSPKRECTOOLARGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKREFNOTSUPP)",{[=]() {throw SpiceSPKREFNOTSUPP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKSTRUCTUREERROR)",{[=]() {throw SpiceSPKSTRUCTUREERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKTYPENOTSUPP)",{[=]() {throw SpiceSPKTYPENOTSUPP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPKTYPENOTSUPPORTD)",{[=]() {throw SpiceSPKTYPENOTSUPPORTD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SPURIOUSKEYWORD)",{[=]() {throw SpiceSPURIOUSKEYWORD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STFULL)",{[=]() {throw SpiceSTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STRINGTOOSHORT)",{[=]() {throw SpiceSTRINGTOOSHORT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STRINGTOOSMALL)",{[=]() {throw SpiceSTRINGTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(STRINGTRUNCATED)",{[=]() {throw SpiceSTRINGTRUNCATED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SUBORBITAL)",{[=]() {throw SpiceSUBORBITAL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SUBPOINTNOTFOUND)",{[=]() {throw SpiceSUBPOINTNOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SYNTAXERROR)",{[=]() {throw SpiceSYNTAXERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(SYSTEMCALLFAILED)",{[=]() {throw SpiceSYSTEMCALLFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TABLENOTLOADED)",{[=]() {throw SpiceTABLENOTLOADED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMECONFLICT)",{[=]() {throw SpiceTIMECONFLICT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMEOUTOFBOUNDS)",{[=]() {throw SpiceTIMEOUTOFBOUNDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMESDONTMATCH)",{[=]() {throw SpiceTIMESDONTMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMESOUTOFORDER)",{[=]() {throw SpiceTIMESOUTOFORDER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TIMEZONEERROR)",{[=]() {throw SpiceTIMEZONEERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWINPUTLINES)",{[=]() {throw SpiceTOOFEWINPUTLINES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWPACKETS)",{[=]() {throw SpiceTOOFEWPACKETS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWPLATES)",{[=]() {throw SpiceTOOFEWPLATES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWSTATES)",{[=]() {throw SpiceTOOFEWSTATES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWVERTICES)",{[=]() {throw SpiceTOOFEWVERTICES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOFEWWINDOWS)",{[=]() {throw SpiceTOOFEWWINDOWS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYBASEFRAMES)",{[=]() {throw SpiceTOOMANYBASEFRAMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYFIELDS)",{[=]() {throw SpiceTOOMANYFIELDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYFILESOPEN)",{[=]() {throw SpiceTOOMANYFILESOPEN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYHITS)",{[=]() {throw SpiceTOOMANYHITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYITERATIONS)",{[=]() {throw SpiceTOOMANYITERATIONS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYKEYWORDS)",{[=]() {throw SpiceTOOMANYKEYWORDS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPAIRS)",{[=]() {throw SpiceTOOMANYPAIRS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPARTS)",{[=]() {throw SpiceTOOMANYPARTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPEAS)",{[=]() {throw SpiceTOOMANYPEAS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYPLATES)",{[=]() {throw SpiceTOOMANYPLATES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYSURFACES)",{[=]() {throw SpiceTOOMANYSURFACES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYVERTICES)",{[=]() {throw SpiceTOOMANYVERTICES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TOOMANYWATCHES)",{[=]() {throw SpiceTOOMANYWATCHES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TRANSFERFILE)",{[=]() {throw SpiceTRANSFERFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TRANSFERFORMAT)",{[=]() {throw SpiceTRANSFERFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TWOSCLKFILENAMES)",{[=]() {throw SpiceTWOSCLKFILENAMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TYPEMISMATCH)",{[=]() {throw SpiceTYPEMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TYPENOTSUPPORTED)",{[=]() {throw SpiceTYPENOTSUPPORTED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(TYPESMISMATCH)",{[=]() {throw SpiceTYPESMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNALLOCATEDNODE)",{[=]() {throw SpiceUNALLOCATEDNODE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNBALANCEDGROUP)",{[=]() {throw SpiceUNBALANCEDGROUP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNBALANCEDPAIR)",{[=]() {throw SpiceUNBALANCEDPAIR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNDEFINEDFRAME)",{[=]() {throw SpiceUNDEFINEDFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNEQUALTIMESTEP)",{[=]() {throw SpiceUNEQUALTIMESTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNINITIALIZED)",{[=]() {throw SpiceUNINITIALIZED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNINITIALIZEDHASH)",{[=]() {throw SpiceUNINITIALIZEDHASH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNINITIALIZEDVALUE)",{[=]() {throw SpiceUNINITIALIZEDVALUE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNITSMISSING)",{[=]() {throw SpiceUNITSMISSING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNITSNOTREC)",{[=]() {throw SpiceUNITSNOTREC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNONWNTIMESYSTEM)",{[=]() {throw SpiceUNKNONWNTIMESYSTEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNBFF)",{[=]() {throw SpiceUNKNOWNBFF(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNCKMETA)",{[=]() {throw SpiceUNKNOWNCKMETA(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNCOMPARE)",{[=]() {throw SpiceUNKNOWNCOMPARE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNDATATYPE)",{[=]() {throw SpiceUNKNOWNDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFILARC)",{[=]() {throw SpiceUNKNOWNFILARC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFRAME)",{[=]() {throw SpiceUNKNOWNFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFRAMESPEC)",{[=]() {throw SpiceUNKNOWNFRAMESPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNFRAMETYPE)",{[=]() {throw SpiceUNKNOWNFRAMETYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNID)",{[=]() {throw SpiceUNKNOWNID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNINCLUSION)",{[=]() {throw SpiceUNKNOWNINCLUSION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNINDEXTYPE)",{[=]() {throw SpiceUNKNOWNINDEXTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNKERNELTYPE)",{[=]() {throw SpiceUNKNOWNKERNELTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNKEY)",{[=]() {throw SpiceUNKNOWNKEY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNMETAITEM)",{[=]() {throw SpiceUNKNOWNMETAITEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNMODE)",{[=]() {throw SpiceUNKNOWNMODE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNOP)",{[=]() {throw SpiceUNKNOWNOP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNPCKTYPE)",{[=]() {throw SpiceUNKNOWNPCKTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNREFDIR)",{[=]() {throw SpiceUNKNOWNREFDIR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNSPKTYPE)",{[=]() {throw SpiceUNKNOWNSPKTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNSYSTEM)",{[=]() {throw SpiceUNKNOWNSYSTEM(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNTYPE)",{[=]() {throw SpiceUNKNOWNTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNKNOWNUNITS)",{[=]() {throw SpiceUNKNOWNUNITS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNMATCHENDPTS)",{[=]() {throw SpiceUNMATCHENDPTS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNNATURALACT)",{[=]() {throw SpiceUNNATURALACT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNNATURALRELATION)",{[=]() {throw SpiceUNNATURALRELATION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNORDEREDREFS)",{[=]() {throw SpiceUNORDEREDREFS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNORDEREDTIMES)",{[=]() {throw SpiceUNORDEREDTIMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNPARSEDQUERY)",{[=]() {throw SpiceUNPARSEDQUERY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNPARSEDTIME)",{[=]() {throw SpiceUNPARSEDTIME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNAPPFLAG)",{[=]() {throw SpiceUNRECOGNAPPFLAG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNDATATYPE)",{[=]() {throw SpiceUNRECOGNDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNDELIMITER)",{[=]() {throw SpiceUNRECOGNDELIMITER(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZABLEFILE)",{[=]() {throw SpiceUNRECOGNIZABLEFILE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDACTION)",{[=]() {throw SpiceUNRECOGNIZEDACTION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDFORMAT)",{[=]() {throw SpiceUNRECOGNIZEDFORMAT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDFRAME)",{[=]() {throw SpiceUNRECOGNIZEDFRAME(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNIZEDTYPE)",{[=]() {throw SpiceUNRECOGNIZEDTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRECOGNPRECTYPE)",{[=]() {throw SpiceUNRECOGNPRECTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRESOLVEDNAMES)",{[=]() {throw SpiceUNRESOLVEDNAMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNRESOLVEDTIMES)",{[=]() {throw SpiceUNRESOLVEDTIMES(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDARCH)",{[=]() {throw SpiceUNSUPPORTEDARCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDBFF)",{[=]() {throw SpiceUNSUPPORTEDBFF(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDMETHOD)",{[=]() {throw SpiceUNSUPPORTEDMETHOD(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNSUPPORTEDSPEC)",{[=]() {throw SpiceUNSUPPORTEDSPEC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UNTITLEDHELP)",{[=]() {throw SpiceUNTITLEDHELP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UPDATEPENDING)",{[=]() {throw SpiceUPDATEPENDING(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(USAGEERROR)",{[=]() {throw SpiceUSAGEERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(UTFULL)",{[=]() {throw SpiceUTFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VALUEOUTOFRANGE)",{[=]() {throw SpiceVALUEOUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VALUETABLEFULL)",{[=]() {throw SpiceVALUETABLEFULL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VARIABLENOTFOUND)",{[=]() {throw SpiceVARIABLENOTFOUND(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VARNAMETOOLONG)",{[=]() {throw SpiceVARNAMETOOLONG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VECTORTOOBIG)",{[=]() {throw SpiceVECTORTOOBIG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VERSIONMISMATCH)",{[=]() {throw SpiceVERSIONMISMATCH(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VERTEXNOTINGRID)",{[=]() {throw SpiceVERTEXNOTINGRID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(VOXELGRIDTOOBIG)",{[=]() {throw SpiceVOXELGRIDTOOBIG(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WIDTHTOOSMALL)",{[=]() {throw SpiceWIDTHTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WINDOWEXCESS)",{[=]() {throw SpiceWINDOWEXCESS(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WINDOWSTOOSMALL)",{[=]() {throw SpiceWINDOWSTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WINDOWTOOSMALL)",{[=]() {throw SpiceWINDOWTOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WORKSPACETOOSMALL)",{[=]() {throw SpiceWORKSPACETOOSMALL(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRITEERROR)",{[=]() {throw SpiceWRITEERROR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRITEFAILED)",{[=]() {throw SpiceWRITEFAILED(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGARCHITECTURE)",{[=]() {throw SpiceWRONGARCHITECTURE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGCKTYPE)",{[=]() {throw SpiceWRONGCKTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGCONIC)",{[=]() {throw SpiceWRONGCONIC(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGDATATYPE)",{[=]() {throw SpiceWRONGDATATYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGSEGMENT)",{[=]() {throw SpiceWRONGSEGMENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(WRONGSPKTYPE)",{[=]() {throw SpiceWRONGSPKTYPE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(YEAROUTOFRANGE)",{[=]() {throw SpiceYEAROUTOFRANGE(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROBORESIGHT)",{[=]() {throw SpiceZEROBORESIGHT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROBOUNDSEXTENT)",{[=]() {throw SpiceZEROBOUNDSEXTENT(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROFRAMEID)",{[=]() {throw SpiceZEROFRAMEID(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROLENGTHCOLUMN)",{[=]() {throw SpiceZEROLENGTHCOLUMN(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROPOSITION)",{[=]() {throw SpiceZEROPOSITION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROQUATERNION)",{[=]() {throw SpiceZEROQUATERNION(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROSTEP)",{[=]() {throw SpiceZEROSTEP(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROVECTOR)",{[=]() {throw SpiceZEROVECTOR(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZEROVELOCITY)",{[=]() {throw SpiceZEROVELOCITY(shortMessage, explanation, longMessage, traceback);}}},
    {"SPICE(ZZHOLDDGETFAILED)",{[=]() {throw SpiceZZHOLDDGETFAILED(shortMessage,explanation, longMessage,traceback);}}}
    };
    // clang-format on

    auto it = spiceExceptions.find( shortMessage );
    if( it != spiceExceptions.end( ) )
    {
        it->second( );
    }
    else
    {
        throw SpiceError( shortMessage, explanation, longMessage, traceback );
    }
}

}  // namespace spice_interface
}  // namespace tudat
