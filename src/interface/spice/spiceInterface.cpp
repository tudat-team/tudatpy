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
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/paths.hpp"

#include <math.h>

namespace tudat {
namespace spice_interface {
using Eigen::Vector6d;

std::string getCorrectedTargetBodyName(
        const std::string &targetBodyName )
{
    std::string correctedTargetBodyName;
    if( targetBodyName == "Mercury" ||
            targetBodyName == "Venus" ||
            targetBodyName == "MERCURY" ||
            targetBodyName == "VENUS" ||
            targetBodyName == "mercury" ||
            targetBodyName == "venus" )
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
double convertJulianDateToEphemerisTime(const double julianDate) {
    return (julianDate - j2000_c()) * spd_c();
}

//! Convert ephemeris time (equivalent to TDB) to a Julian date.
double convertEphemerisTimeToJulianDate(const double ephemerisTime) {
    return j2000_c() + (ephemerisTime) / spd_c();
}

//! Converts a date string to ephemeris time.
double convertDateStringToEphemerisTime(const std::string &dateString) {
    double ephemerisTime = 0.0;
    str2et_c(dateString.c_str(), &ephemerisTime);
    return ephemerisTime;
}

//! Get Cartesian state of a body, as observed from another body.
Vector6d getBodyCartesianStateAtEpoch(
        const std::string &targetBodyName, const std::string &observerBodyName,
        const std::string &referenceFrameName, const std::string &aberrationCorrections,
        const double ephemerisTime) {


    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving Cartesian state from Spice, input time is " + std::to_string(ephemerisTime) );
    }
    // Declare variables for cartesian state and light-time to be determined by Spice.
    double stateAtEpoch[6];
    double lightTime;

    // Call Spice function to calculate state and light-time.
    spkezr_c(getCorrectedTargetBodyName( targetBodyName ).c_str(), ephemerisTime, referenceFrameName.c_str(),
             aberrationCorrections.c_str(),
             getCorrectedTargetBodyName( observerBodyName ).c_str(), stateAtEpoch,
             &lightTime);

    // Put result in Eigen Vector.
    Vector6d cartesianStateVector;
    for (unsigned int i = 0; i < 6; i++) {
        cartesianStateVector(i) = stateAtEpoch[i];
    }

    // Convert from km(/s) to m(/s).
    return unit_conversions::convertKilometersToMeters<Vector6d>(
                cartesianStateVector);
}

//! Get Cartesian position of a body, as observed from another body.
Eigen::Vector3d getBodyCartesianPositionAtEpoch(const std::string &targetBodyName,
                                                const std::string &observerBodyName,
                                                const std::string &referenceFrameName,
                                                const std::string &aberrationCorrections,
                                                const double ephemerisTime) {

    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving Cartesian position from Spice, input time is " + std::to_string(ephemerisTime) );
    }
    // Declare variables for cartesian position and light-time to be determined by Spice.
    double positionAtEpoch[3];
    double lightTime;

    // Call Spice function to calculate position and light-time.
    spkpos_c(getCorrectedTargetBodyName( targetBodyName ).c_str(), ephemerisTime, referenceFrameName.c_str(),
             aberrationCorrections.c_str(),
             getCorrectedTargetBodyName( observerBodyName ).c_str(), positionAtEpoch,
             &lightTime);

    // Put result in Eigen Vector.
    Eigen::Vector3d cartesianPositionVector;
    for (unsigned int i = 0; i < 3; i++) {
        cartesianPositionVector(i) = positionAtEpoch[i];
    }

    // Convert from km to m.
    return unit_conversions::convertKilometersToMeters<Eigen::Vector3d>(
                cartesianPositionVector);
}

//! Get Cartesian state of a satellite from its two-line element set at a specified epoch.
Vector6d getCartesianStateFromTleAtEpoch(double epoch, std::shared_ptr<ephemerides::Tle> tle) {
    if( !( epoch == epoch ))
    {
        throw std::invalid_argument( "Error when retrieving TLE from Spice, input time is " + std::to_string(epoch) );
    }

    // Physical constants used by CSpice's implementation of SGP4.
    double physicalConstants[8] = {1.082616E-3, -2.53881E-6, -1.65597E-6, 7.43669161e-2, 120.0, 78.0, 6378.135, 1.0};

    // Declare variable that will hold the state as returned by Spice.
    double stateAtEpoch[6];

    // TODO: convert elements to units required by CSpice (?)
    double elements[10];
    elements[0] = 0.0;// This element is mandatory as input to ev2lin_ but not used internally (used to be accessed in SGP).
    elements[1] = 0.0;// Idem dito.
    elements[2] = tle->getBStar();
    elements[3] = tle->getInclination();
    elements[4] = tle->getRightAscension();
    elements[5] = tle->getEccentricity();
    elements[6] = tle->getArgOfPerigee();
    elements[7] = tle->getMeanAnomaly();
    elements[8] = tle->getMeanMotion();
    elements[9] = tle->getEpoch();// TLE ephemeris epoch in seconds since J2000

    // Call Spice function. Return value is always 0, so no need to save it.
    ev2lin_(&epoch, physicalConstants, elements, stateAtEpoch);

    // Put result in Eigen Vector.
    Vector6d cartesianStateVector;
    for (unsigned int i = 0; i < 6; i++) {
        cartesianStateVector(i) = stateAtEpoch[i];
    }

    // Convert from km to m.
    return unit_conversions::convertKilometersToMeters<Vector6d>(cartesianStateVector);
}

//! Compute quaternion of rotation between two frames.
Eigen::Quaterniond computeRotationQuaternionBetweenFrames(const std::string &originalFrame,
                                                          const std::string &newFrame,
                                                          const double ephemerisTime) {
    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving rotation quaternion from Spice, input time is " + std::to_string(ephemerisTime) );
    }

    // Declare rotation matrix.
    double rotationArray[3][3];

    // Calculate rotation matrix.
    pxform_c(originalFrame.c_str(), newFrame.c_str(), ephemerisTime, rotationArray);

    // Put rotation matrix in Eigen Matrix3d.
    Eigen::Matrix3d rotationMatrix;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            rotationMatrix(i, j) = rotationArray[i][j];
        }
    }

    // Convert matrix3d to Quaternion.
    return Eigen::Quaterniond(rotationMatrix);
}

Eigen::Matrix3d computeRotationMatrixBetweenFrames(const std::string &originalFrame,
                                                   const std::string &newFrame,
                                                   const double ephemerisTime)
{
    return Eigen::Matrix3d( computeRotationQuaternionBetweenFrames(
                                originalFrame, newFrame, ephemerisTime ) );
}

//! Compute rotation matrix for state vector between two frames.
Eigen::Matrix6d computeStateRotationMatrixBetweenFrames(const std::string &originalFrame,
                                                          const std::string &newFrame,
                                                          const double ephemerisTime) {
    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving state rotation matrix from Spice, input time is " + std::to_string(ephemerisTime) );
    }

    double stateTransition[6][6];

    // Calculate state transition matrix.
    sxform_c(originalFrame.c_str(), newFrame.c_str(), ephemerisTime, stateTransition);

    // Put rotation matrix in Eigen Matrix6d
    Eigen::Matrix6d stateTransitionMatrix = Eigen::Matrix6d::Zero();
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            stateTransitionMatrix(i, j) = stateTransition[i][j];
        }
    }

    return stateTransitionMatrix;
}

//! Computes time derivative of rotation matrix between two frames.
Eigen::Matrix3d computeRotationMatrixDerivativeBetweenFrames(const std::string &originalFrame,
                                                             const std::string &newFrame,
                                                             const double ephemerisTime) {

    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving rotation matrix derivative from Spice, input time is " + std::to_string(ephemerisTime) );
    }

    double stateTransition[6][6];

    // Calculate state transition matrix.
    sxform_c(originalFrame.c_str(), newFrame.c_str(), ephemerisTime, stateTransition);

    // Put rotation matrix derivative in Eigen Matrix3d
    Eigen::Matrix3d matrixDerivative = Eigen::Matrix3d::Zero();
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            matrixDerivative(i, j) = stateTransition[i + 3][j];
        }
    }

    return matrixDerivative;
}

//! Computes the angular velocity of one frame w.r.t. to another frame.
Eigen::Vector3d getAngularVelocityVectorOfFrameInOriginalFrame(const std::string &originalFrame,
                                                               const std::string &newFrame,
                                                               const double ephemerisTime) {

    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving angular velocity from Spice, input time is " + std::to_string(ephemerisTime) );
    }

    double stateTransition[6][6];

    // Calculate state transition matrix.
    sxform_c(originalFrame.c_str(), newFrame.c_str(), ephemerisTime, stateTransition);

    double rotation[3][3];
    double angularVelocity[3];

    // Calculate angular velocity vector.
    xf2rav_c(stateTransition, rotation, angularVelocity);

    return (Eigen::Vector3d() << angularVelocity[0], angularVelocity[1], angularVelocity[2]).finished();
}

std::pair<Eigen::Quaterniond, Eigen::Matrix3d> computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames(
        const std::string &originalFrame, const std::string &newFrame, const double ephemerisTime) {
    double stateTransition[6][6];

    if( !( ephemerisTime == ephemerisTime )  )
    {
        throw std::invalid_argument( "Error when retrieving rotational state from Spice, input time is " + std::to_string(ephemerisTime) );
    }

    sxform_c(originalFrame.c_str(), newFrame.c_str(), ephemerisTime, stateTransition);

    Eigen::Matrix3d matrixDerivative;
    Eigen::Matrix3d rotationMatrix;

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            rotationMatrix(i, j) = stateTransition[i][j];
            matrixDerivative(i, j) = stateTransition[i + 3][j];
        }
    }
    return std::make_pair(Eigen::Quaterniond(rotationMatrix), matrixDerivative);
}

//! Get property of a body from Spice.
std::vector<double> getBodyProperties(const std::string &body, const std::string &property,
                                      const int maximumNumberOfValues) {
    // Delcare variable in which raw result is to be put by Spice function.
    double* propertyArray = new double[maximumNumberOfValues];

    // Call Spice function to retrieve property.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c(body.c_str(), property.c_str(), maximumNumberOfValues, &numberOfReturnedParameters,
             propertyArray);

    // Put result in STL vector.
    std::vector<double> bodyProperties;
    bodyProperties.resize(numberOfReturnedParameters);
    for (int i = 0; i < numberOfReturnedParameters; i++) {
        bodyProperties.at(i) = propertyArray[i];
    }
    delete [] propertyArray;
    return bodyProperties;
}

//! Get gravitational parameter of a body.
double getBodyGravitationalParameter(const std::string &body) {
    // Delcare variable in which raw result is to be put by Spice function.
    double gravitationalParameter[1];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c(body.c_str(), "GM", 1, &numberOfReturnedParameters, gravitationalParameter);

    // Convert from km^3/s^2 to m^3/s^2
    return unit_conversions::convertKilometersToMeters<double>(
                unit_conversions::convertKilometersToMeters<double>(
                    unit_conversions::convertKilometersToMeters<double>(
                        gravitationalParameter[0])));
}

//! Get the (arithmetic) mean of the three principal axes of the tri-axial ellipsoid shape.
double getAverageRadius(const std::string &body) {
    // Delcare variable in which raw result is to be put by Spice function.
    double radii[3];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c(body.c_str(), "RADII", 3, &numberOfReturnedParameters, radii);

    // Compute average and convert from km to m.
    return unit_conversions::convertKilometersToMeters<double>(
                radii[0] + radii[1] + radii[2])
            / 3.0;
}

//! Get the (arithmetic) mean of the two equatorial axes of the tri-axial ellipsoid shape.
double getAverageEquatorialRadius( const std::string& body )
{
    // Declare variable in which raw result is to be put by Spice function.
    double radii[3];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str(), "RADII", 3, &numberOfReturnedParameters, radii );

    // Compute average and convert from km to m.
    return unit_conversions::convertKilometersToMeters< double >(
            radii[0] + radii[1] ) / 2.0;
}

//! Get the polar radius of the tri-axial ellipsoid shape.
double getPolarRadius( const std::string& body )
{
    // Declare variable in which raw result is to be put by Spice function.
    double radii[3];

    // Call Spice function to retrieve gravitational parameter.
    SpiceInt numberOfReturnedParameters;
    bodvrd_c( body.c_str(), "RADII", 3, &numberOfReturnedParameters, radii );

    // Compute average and convert from km to m.
    return unit_conversions::convertKilometersToMeters< double >(radii[2] );
}

//! Convert a body name to its NAIF identification number.
int convertBodyNameToNaifId(const std::string &bodyName) {
    // Convert body name to NAIF ID number.
    SpiceInt bodyNaifId;
    SpiceBoolean isIdFound;
    bods2c_c(bodyName.c_str(), &bodyNaifId, &isIdFound);

    // Convert SpiceInt (typedef for long) to int and return.
    return static_cast<int>(bodyNaifId);
}

//! Convert a NAIF identification number to its body name.
std::string convertNaifIdToBodyName( int bodyNaifId )
{
    // Maximum SPICE name length is 32. Therefore, a name length of 33 is used (+1 for null terminator)
    SpiceChar bodyName[33];

    bodc2s_c( bodyNaifId, 33, bodyName );

    // Convert SpiceChar to std::string
    return static_cast< std::string >( bodyName );
}

//! Check if a certain property of a body is in the kernel pool.
bool checkBodyPropertyInKernelPool(const std::string &bodyName, const std::string &bodyProperty) {
    // Convert body name to NAIF ID.
    const int naifId = convertBodyNameToNaifId(bodyName);

    // Determine if property is in pool.
    SpiceBoolean isPropertyInPool = bodfnd_c(naifId, bodyProperty.c_str());
    return static_cast<bool>(isPropertyInPool);
}

//! Load a Spice kernel.
void loadSpiceKernelInTudat(const std::string &fileName) {
    furnsh_c(fileName.c_str());
}

//! Get the amount of loaded Spice kernels.
int getTotalCountOfKernelsLoaded() {
    SpiceInt count;
    ktotal_c("ALL", &count);
    return count;
}

//! Clear all Spice kernels.
void clearSpiceKernels() { kclear_c(); }

//! Get all standard Spice kernels used in tudat.
std::vector<std::string> getStandardSpiceKernels(const std::vector<std::string> alternativeEphemerisKernels) {
    std::vector<std::string> standardSpiceKernels;

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

void loadStandardSpiceKernels(const std::vector<std::string> alternativeEphemerisKernels) {

    std::string kernelPath = paths::getSpiceKernelPath();
    loadSpiceKernelInTudat(kernelPath + "/pck00010.tpc");
//    loadSpiceKernelInTudat(kernelPath + "/gm_de431.tpc");
    loadSpiceKernelInTudat(kernelPath + "/inpop19a_TDB_m100_p100_spice.tpc");
    loadSpiceKernelInTudat(kernelPath + "/NOE-4-2020.tpc");
    loadSpiceKernelInTudat(kernelPath + "/NOE-5-2021.tpc");
    loadSpiceKernelInTudat(kernelPath + "/NOE-6-2018-MAIN-v2.tpc");

    if (alternativeEphemerisKernels.size() == 0)
    {

        loadSpiceKernelInTudat(kernelPath + "/codes_300ast_20100725.bsp");
        loadSpiceKernelInTudat(kernelPath + "/codes_300ast_20100725.tf");
        loadSpiceKernelInTudat(kernelPath + "/inpop19a_TDB_m100_p100_spice.bsp");
        loadSpiceKernelInTudat(kernelPath + "/NOE-4-2020.bsp");
        loadSpiceKernelInTudat(kernelPath + "/NOE-5-2021.bsp");
        loadSpiceKernelInTudat(kernelPath + "/NOE-6-2018-MAIN-v2.bsp");
        loadSpiceKernelInTudat(kernelPath + "/juice_mat_crema_4_0_20220601_20330626_v01.bsp");

    }
    else
    {
        for (unsigned int i = 0; i < alternativeEphemerisKernels.size(); i++)
        {
            loadSpiceKernelInTudat(alternativeEphemerisKernels.at(i));
        }
    }
    loadSpiceKernelInTudat(kernelPath + "/naif0012.tls");
}

}// namespace spice_interface
}// namespace tudat
