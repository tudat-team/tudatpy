/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEROTATIONMODEL_H
#define TUDAT_CREATEROTATIONMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <memory>

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/interface/sofa/earthOrientation.h"

namespace tudat
{

namespace simulation_setup
{

//List of rotation models available in simulations
/*
 *  List of rotation models available in simulations. Rotation models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
//! @get_docstring(RotationModelType.__docstring__)
enum RotationModelType
{
    simple_rotation_model,
    spice_rotation_model,
    gcrs_to_itrs_rotation_model,
    synchronous_rotation_model,
    planetary_rotation_model,
    tabulated_rotation_model,
    aerodynamic_angle_based_rotation_model,
    pitch_trim_rotation_model,
    body_fixed_direction_based_rotation_model,
    orbital_state_based_rotation_model
};

//Class for providing settings for rotation model.
/*
 *  Class for providing settings for automatic rotation model creation. This class is a
 *  functional (base) class for settings of rotation models that require no information in
 *  addition to their type. Rotation model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */

//! @get_docstring(RotationModelSettings.__docstring__)
class RotationModelSettings
{
public:

    //Constructor, sets type of rotation model.
    /*
     *  Constructor, sets type of rotation model and base and target frame identifiers.
     *  Settings for rotation models requiring additional information should be defined in a
     *  derived class.
     *  \param rotationType Type of rotation model that is to be created.
     *  \param originalFrame Base frame of rotation model.
     *  \param targetFrame Target frame of rotation model.
     */

    RotationModelSettings( const RotationModelType rotationType,
                           const std::string& originalFrame,
                           const std::string& targetFrame ):
        rotationType_( rotationType ), originalFrame_( originalFrame ),
        targetFrame_( targetFrame ){ }

    //Destructor.
    virtual ~RotationModelSettings( ){ }

    //Function to return the type of rotation model that is to be created.
    /*
     *  Function to return the type of rotation model that is to be created.
     *  \return Type of rotation model that is to be created.
     */

    RotationModelType getRotationType( ){ return rotationType_; }

    //Function to return the base frame of rotation model.
    /*
     *  Function to return the base frame of rotation model.
     *  \return Base frame of rotation model.
     */

    std::string getOriginalFrame( ){ return originalFrame_; }

    //Function to return the target frame of rotation model.
    /*
     *  Function to return the target frame of rotation model.
     *  \return Target frame of rotation model.
     */

    std::string getTargetFrame( ){ return targetFrame_; }


    //Function to reset the orientation of the base frame.
    /*
     * Function to reset the orientation of the base frame.
     * \param originalFrame New base frame orientation
     */

    void resetOriginalFrame( const std::string& originalFrame )
    {
        originalFrame_ = originalFrame;
    }

protected:

    //Type of rotation model that is to be created.
    RotationModelType rotationType_;

    //Target frame of rotation model.
    std::string originalFrame_;

    //Base frame of rotation model.
    std::string targetFrame_;

};

//RotationModelSettings derived class for defining settings of a simple rotational ephemeris.
class SimpleRotationModelSettings: public RotationModelSettings
{
public:
    //Constructor,
    /*
     *  Constructor, sets simple rotational ephemeris properties.
     *  \param originalFrame Base frame of rotation model.
     *  \param targetFrame Target frame of rotation model.
     *  \param initialOrientation Rotation from base to target frame at initialTime.
     *  \param initialTime Time at which initialOrientation represents the instantaneous rotation.
     *  \param rotationRate Rotation rate of body about its local z-axis.
     */
    SimpleRotationModelSettings( const std::string& originalFrame,
                                 const std::string& targetFrame,
                                 const Eigen::Quaterniond& initialOrientation,
                                 const double initialTime,
                                 const double rotationRate ):
        RotationModelSettings( simple_rotation_model, originalFrame, targetFrame ),
        initialOrientation_( initialOrientation ),
        initialTime_( initialTime ), rotationRate_( rotationRate ){ }

    //Function to return rotation from base to target frame at initialTime.
    /*
     *  Function to return rotation from base to target frame at initialTime.
     *  \return Rotation from base to target frame at initialTime.
     */
    Eigen::Quaterniond getInitialOrientation( ){ return initialOrientation_; }

    //Function to return time at which initialOrientation represents the instantaneous rotation.
    /*
     *  Function to return time at which initialOrientation represents the instantaneous rotation.
     *  \return Time at which initialOrientation represents the instantaneous rotation.
     */
    double getInitialTime( ){ return initialTime_; }

    //Function to return rotation rate of body about its local z-axis.
    /*
     *  Function to return rotation rate of body about its local z-axis.
     *  \return Rotation rate of body about its local z-axis.
     */
    double getRotationRate( ){ return rotationRate_; }

private:

    // Rotation from base to target frame at initialTime.
    Eigen::Quaterniond initialOrientation_;

    //Time at which initialOrientation represents the instantaneous rotation.
    double initialTime_;

    //Rotation rate of body about its local z-axis.
    double rotationRate_;
};

//#ifdef TUDAT_BUILD_WITH_SOFA_INTERFACE

//Struct that holds settings for EOP short-period variation
struct EopCorrectionSettings
{
    //Constructor
    /*
     *  Constructor
     *  \param conversionFactor Conversion factor to be used for amplitudes to multiply input values, typically for unit
     *  conversion purposes.
     *  \param minimumAmplitude Minimum amplitude that is read from files and considered in calculations.
     *  \param amplitudesFiles List of files with amplitudes for corrections
     *  \param argumentMultipliersFile Fundamental argument multiplier for corrections
     */
    EopCorrectionSettings(
            const double conversionFactor,
            const double minimumAmplitude,
            const std::vector< std::string >& amplitudesFiles,
            const std::vector< std::string >& argumentMultipliersFile ):
        conversionFactor_( conversionFactor ), minimumAmplitude_( minimumAmplitude ),
        amplitudesFiles_( amplitudesFiles ), argumentMultipliersFile_( argumentMultipliersFile ){ }

    //Conversion factor to be used for amplitudes to multiply input values
    double conversionFactor_;

    //Minimum amplitude that is read from files and considered in calculations.
    double minimumAmplitude_;

    //List of files with amplitudes for corrections
    std::vector< std::string > amplitudesFiles_;

    //Fundamental argument multiplier for corrections
    std::vector< std::string > argumentMultipliersFile_;
};

//Settings for creating a GCRS<->ITRS rotation model
class GcrsToItrsRotationModelSettings: public RotationModelSettings
{
public:

    //Constructor
    /*
     * \param baseFrameName Name of base frame (typically GCRS, which is default)
     * \param timeScale Time scale in which input to the rotation model class is provided, default TDB
     * \param nutationTheory IAU precession-nutation theory that is to be used.
     * \param eopFile Name of EOP file that is to be used
     * \param eopFileFormat Identifier for file format that is provided
     * \param ut1CorrectionSettings Settings for short-period UT1-UTC variations
     * \param polarMotionCorrectionSettings Settings for short-period polar motion variations
     */
    GcrsToItrsRotationModelSettings(
            const basic_astrodynamics::IAUConventions nutationTheory = basic_astrodynamics::iau_2006,
            const std::string baseFrameName = "GCRS",
            const std::string& eopFile = paths::getEarthOrientationDataFilesPath( ) + "/eopc04_08_IAU2000.62-now.txt",
            const basic_astrodynamics::TimeScales inputTimeScale = basic_astrodynamics::tdb_scale,
            const std::shared_ptr< EopCorrectionSettings > ut1CorrectionSettings =
            std::make_shared< EopCorrectionSettings >(
                1.0E-6, 0.0, std::vector< std::string >{
                    paths::getEarthOrientationDataFilesPath( ) + "/utcLibrationAmplitudes.txt",
                    paths::getEarthOrientationDataFilesPath( ) + "/utcOceanTidesAmplitudes.txt" },
                std::vector< std::string >{
                    paths::getEarthOrientationDataFilesPath( ) +
                    "/utcLibrationFundamentalArgumentMultipliers.txt",
                    paths::getEarthOrientationDataFilesPath( ) +
                    "/utcOceanTidesFundamentalArgumentMultipliers.txt" } ),
            const std::shared_ptr< EopCorrectionSettings > polarMotionCorrectionSettings =
            std::make_shared< EopCorrectionSettings >(
                unit_conversions::convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0, std::vector< std::string >{
                    paths::getEarthOrientationDataFilesPath( ) +
                    "/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
                    paths::getEarthOrientationDataFilesPath( ) +
                    "/polarMotionOceanTidesAmplitudes.txt", },
                std::vector< std::string >{
                    paths::getEarthOrientationDataFilesPath( ) +
                    "/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
                    paths::getEarthOrientationDataFilesPath( ) +
                    "/polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } ) ):
        RotationModelSettings( gcrs_to_itrs_rotation_model, baseFrameName, "ITRS" ),
        inputTimeScale_( inputTimeScale ), nutationTheory_( nutationTheory ), eopFile_( eopFile ),
        eopFileFormat_( "C04" ), ut1CorrectionSettings_( ut1CorrectionSettings ),
        polarMotionCorrectionSettings_( polarMotionCorrectionSettings ){ }

    //Destructor
    ~GcrsToItrsRotationModelSettings( ){ }

    //Function to retrieve the time scale in which input to the rotation model class is provided
    /*
     * Function to retrieve the time scale in which input to the rotation model class is provided
     * \return Time scale in which input to the rotation model class is provided
     */
    basic_astrodynamics::TimeScales getInputTimeScale( )
    {
        return inputTimeScale_;
    }
    //Function to retrieve the IAU precession-nutation theory that is to be used
    /*
     * Function to retrieve the IAU precession-nutation theory that is to be used
     * \return IAU precession-nutation theory that is to be used
     */
    basic_astrodynamics::IAUConventions getNutationTheory( )
    {
        return nutationTheory_;
    }

    //Function to retrieve the name of EOP file that is to be used
    /*
     * Function to retrieve the name of EOP file that is to be used
     * \return Name of EOP file that is to be used
     */
    std::string getEopFile( )
    {
        return eopFile_;
    }

    //Function to retrieve the identifier for file format that is provided
    /*
     * Function to retrieve the identifier for file format that is provided
     * \return Identifier for file format that is provided
     */
    std::string getEopFileFormat( )
    {
        return eopFileFormat_;
    }

    //Function to retrieve the settings for short-period UT1-UTC variations
    /*
     * Function to retrieve the settings for short-period UT1-UTC variations
     * \return Settings for short-period UT1-UTC variations
     */
    std::shared_ptr< EopCorrectionSettings > getUt1CorrectionSettings( )
    {
        return ut1CorrectionSettings_;
    }

    //Function to retrieve the settings for short-period polar motion variations
    /*
     * Function to retrieve the Settings for short-period polar motion variations
     * \return settings for short-period polar motion variations
     */
    std::shared_ptr< EopCorrectionSettings > getPolarMotionCorrectionSettings( )
    {
        return polarMotionCorrectionSettings_;
    }

private:

    //Time scale in which input to the rotation model class is provided
    basic_astrodynamics::TimeScales inputTimeScale_;

    //IAU precession-nutation theory that is to be used
    basic_astrodynamics::IAUConventions nutationTheory_;

    //Name of EOP file that is to be used
    std::string eopFile_;

    //Identifier for file format that is provided
    std::string eopFileFormat_;

    //Settings for short-period UT1-UTC variations
    std::shared_ptr< EopCorrectionSettings > ut1CorrectionSettings_;

    //Settings for short-period polar motion variations
    std::shared_ptr< EopCorrectionSettings > polarMotionCorrectionSettings_;

};
//#endif


//RotationModelSettings derived class for defining settings of a synchronous rotational ephemeris (body-fixed x-axis always
//pointing to central body; z-axis along r x v (with r and v the position and velocity w.r.t. central body)
class SynchronousRotationModelSettings: public RotationModelSettings
{
public:

    //Constructor
    /*
     * Constructor
     * \param centralBodyName Name of central body to which this body is locked.
     * \param baseFrameOrientation Base frame of rotation model.
     * \param targetFrameOrientation Target frame of rotation model.
     */
    SynchronousRotationModelSettings(
            const std::string& centralBodyName,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation ):
        RotationModelSettings( synchronous_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        centralBodyName_( centralBodyName ){ }

    //Function to retrieve name of central body to which this body is locked.
    /*
     * Function to retrieve name of central body to which this body is locked.
     * \return  Name of central body to which this body is locked.
     */
    std::string getCentralBodyName( )
    {
        return centralBodyName_;
    }

private:

    // Name of central body to which this body is locked.
    std::string centralBodyName_;
};

class TabulatedRotationSettings: public RotationModelSettings
{
public:

    TabulatedRotationSettings(
            const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ):
        RotationModelSettings( tabulated_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        rotationalStateHistory_( rotationalStateHistory ), interpolatorSettings_( interpolatorSettings ){ }

    std::map< double, Eigen::Vector7d > getBodyStateHistory( )
    { return rotationalStateHistory_; }

    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( )
    { return interpolatorSettings_; }

private:

    std::map< double, Eigen::Vector7d > rotationalStateHistory_;

    const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;
};


class AerodynamicAngleRotationSettings: public RotationModelSettings
{
public:

    AerodynamicAngleRotationSettings(
            const std::string& centralBody,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction = nullptr  ):
        RotationModelSettings( aerodynamic_angle_based_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        centralBody_( centralBody ), aerodynamicAngleFunction_( aerodynamicAngleFunction ){ }

    std::string centralBody_;

    std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction_;
};

class PitchTrimRotationSettings: public RotationModelSettings
{
public:

    PitchTrimRotationSettings(
            const std::string& centralBody,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction = nullptr  ):
        RotationModelSettings( pitch_trim_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        centralBody_( centralBody ), sideslipAndBankAngleFunction_( sideslipAndBankAngleFunction ){ }

    std::string centralBody_;

    std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction_;
};

enum InertialDirectionTypes
{
    custom_inertial_direction,
    state_based_inertial_direction,
    bilinear_tangent_inertial_direction
};


class InertialDirectionSettings
{
public:

    InertialDirectionSettings(
            const InertialDirectionTypes inertialDirectionType,
            std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame =
            std::make_pair( ephemerides::inertial_satellite_based_frame, "" ) ):
    inertialDirectionType_( inertialDirectionType ),
    directionFrame_( directionFrame )
    {
        if( inertialDirectionType_ == state_based_inertial_direction &&
                directionFrame.first != ephemerides::inertial_satellite_based_frame )
        {
            throw std::runtime_error( "Error when making inertial direction settings for direction-based rotation model, state-based direction can only be provided in inertial coordinates" );
        }
    }

    virtual ~InertialDirectionSettings( ){ }

    const InertialDirectionTypes inertialDirectionType_;

    std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame_;

};

class CustomInertialDirectionSettings: public InertialDirectionSettings
{
public:

    CustomInertialDirectionSettings(
            const std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction,
            std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame =
            std::make_pair( ephemerides::inertial_satellite_based_frame, "" ) ):
        InertialDirectionSettings( custom_inertial_direction, directionFrame ),
    inertialBodyAxisDirectionFunction_( inertialBodyAxisDirectionFunction ){ }

    ~CustomInertialDirectionSettings( ){ }

    std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction_;
};


class StateBasedInertialDirectionSettings: public InertialDirectionSettings
{
public:

    StateBasedInertialDirectionSettings(
            const std::string& centralBody,
            const bool isColinearWithVelocity,
            const bool directionIsOppositeToVector,
            std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame =
            std::make_pair( ephemerides::inertial_satellite_based_frame, "" )  ):
        InertialDirectionSettings( state_based_inertial_direction, directionFrame ),
        centralBody_( centralBody ),
        isColinearWithVelocity_( isColinearWithVelocity ),
        directionIsOppositeToVector_( directionIsOppositeToVector ){ }

    ~StateBasedInertialDirectionSettings( ){ }

    std::string centralBody_;

    bool isColinearWithVelocity_;

    bool directionIsOppositeToVector_;
};


class BilinearTangentInertialDirectionSettings: public InertialDirectionSettings
{
public:

    BilinearTangentInertialDirectionSettings(
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcDurations,
            const std::vector< Eigen::Vector3d >& constantDirectionTerms,
            const std::vector< Eigen::Vector3d >& linearDirectionTerms,
            std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame =
            std::make_pair( ephemerides::inertial_satellite_based_frame, "" ) ):
        InertialDirectionSettings( bilinear_tangent_inertial_direction, directionFrame ),
    arcStartTimes_( arcStartTimes ),
    arcDurations_( arcDurations ),
    constantDirectionTerms_( constantDirectionTerms ),
    linearDirectionTerms_( linearDirectionTerms ){ }

    ~BilinearTangentInertialDirectionSettings( ){ }

    std::vector< double > arcStartTimes_;
    std::vector< double > arcDurations_;
    std::vector< Eigen::Vector3d > constantDirectionTerms_;
    std::vector< Eigen::Vector3d > linearDirectionTerms_;

};


class BodyFixedDirectionBasedRotationSettings: public RotationModelSettings
{
public:

    BodyFixedDirectionBasedRotationSettings(
            const std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< double( const double ) > freeRotationAngleFunction = nullptr,
            const std::pair< ephemerides::SatelliteBasedFrames, std::string > directionFrame =
            std::make_pair( ephemerides::inertial_satellite_based_frame, "" ) ):
        RotationModelSettings( body_fixed_direction_based_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        inertialDirectionSettings_( std::make_shared< CustomInertialDirectionSettings >(
                                        inertialBodyAxisDirectionFunction, directionFrame ) ),
        freeRotationAngleFunction_( freeRotationAngleFunction ){ }

    BodyFixedDirectionBasedRotationSettings(
            const std::shared_ptr< InertialDirectionSettings > inertialDirectionSettings,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< double( const double ) > freeRotationAngleFunction = nullptr ):
        RotationModelSettings( body_fixed_direction_based_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        inertialDirectionSettings_( inertialDirectionSettings ),
        freeRotationAngleFunction_( freeRotationAngleFunction ){ }

    std::shared_ptr< InertialDirectionSettings > inertialDirectionSettings_;

    std::function< double( const double ) > freeRotationAngleFunction_;
};

class OrbitalStateBasedRotationSettings: public RotationModelSettings
{
public:

    OrbitalStateBasedRotationSettings(
            const std::string& centralBody,
            const bool isColinearWithVelocity,
            const bool directionIsOppositeToVector,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< double( const double ) > freeRotationAngleFunction = nullptr ):
        RotationModelSettings( orbital_state_based_rotation_model, baseFrameOrientation, targetFrameOrientation ),
        centralBody_( centralBody ),
        isColinearWithVelocity_( isColinearWithVelocity ),
        directionIsOppositeToVector_( directionIsOppositeToVector ),
        freeRotationAngleFunction_( freeRotationAngleFunction ){ }

    std::string centralBody_;

    bool isColinearWithVelocity_;

    bool directionIsOppositeToVector_;

    std::function< double( const double ) > freeRotationAngleFunction_;
};

std::function< Eigen::Matrix3d( const double ) > getRotationFunctionFromSatelliteBasedFrame(
        const ephemerides::SatelliteBasedFrames frameId,
        const SystemOfBodies& bodies,
        const std::string& body,
        const std::string& centralBody = "" );

//Function to retrieve a state from one of two functions
/*
 *  Function to retrieve a state from one of two functions, typically from an Ephemeris or a Body object.
 *  \param currentTime Time at which state function is to be evaluated
 *  \param useFirstFunction Boolean defining whether stateFunction1 or stateFunction2 is used
 *  \param stateFunction1 First function returning state as function of time
 *  \param stateFunction2 Second function returning state as function of time
 *  \return Selected function, evaluated at given time
 */
Eigen::Vector6d getStateFromSelectedStateFunction(
        const double currentTime,
        const bool useFirstFunction,
        const std::function< Eigen::Vector6d( const double ) > stateFunction1,
        const std::function< Eigen::Vector6d( const double ) > stateFunction2 );

//Function to create a state function for a body, valid both during propagation, and outside propagation
/*
 * Function to create a state function for a body, valid both during propagation, and outside propagation
 * \param bodies List of body objects
 * \param orbitingBody Body for which state function is to be created
 * \param centralBody Central body w.r.t. which state function is to be created
 * \return Required state function
 */
std::function< Eigen::Vector6d( const double, bool ) > createRelativeStateFunction(
        const SystemOfBodies& bodies,
        const std::string orbitingBody,
        const std::string centralBody );

class PlanetaryRotationModelSettings: public RotationModelSettings
{
public:
    PlanetaryRotationModelSettings( const double angleN,
                                    const double angleJ,
                                    const double anglePsiAtEpoch,
                                    const double anglePsiRateAtEpoch,
                                    const double angleIAtEpoch,
                                    const double angleIRateAtEpoch,
                                    const double anglePhiAtEpoch,
                                    const double anglePhiRateAtEpoch,
                                    const double coreFactor,
                                    const double freeCoreNutationRate,
                                    const std::string originalFrame,
                                    const std::string targetFrame,
                                    const std::string centralBody,
                                    const std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
            ( std::map< double, std::pair< double, double > >( ) ),
                                    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
            ( std::vector< std::map< double, std::pair< double, double > > > ( ) ),
                                    std::vector< std::function< double( const double ) > > timeDependentPhaseCorrectionFunctions =
            ( std::vector< std::function< double( const double ) > >( ) ),
                                    const std::map< double, std::pair< double, double > > rotationRateCorrections =
            ( std::map< double, std::pair< double, double > >( ) ),
                                    const std::map< double, std::pair< double, double > > xPolarMotionCoefficients =
            ( std::map< double, std::pair< double, double > >( ) ),
                                    const std::map< double, std::pair< double, double > > yPolarMotionCoefficients =
            ( std::map< double, std::pair< double, double > >( ) ) ):

        RotationModelSettings( planetary_rotation_model, originalFrame, targetFrame ),
        angleN_( angleN ), angleJ_( angleJ ), anglePsiAtEpoch_( anglePsiAtEpoch ), anglePsiRateAtEpoch_( anglePsiRateAtEpoch ),
        angleIAtEpoch_( angleIAtEpoch ), angleIRateAtEpoch_( angleIRateAtEpoch ), anglePhiAtEpoch_( anglePhiAtEpoch ),
        anglePhiRateAtEpoch_( anglePhiRateAtEpoch ), coreFactor_( coreFactor ), freeCoreNutationRate_( freeCoreNutationRate ),
        meanMotionDirectNutationCorrections_( meanMotionDirectNutationCorrections ),
        meanMotionTimeDependentPhaseNutationCorrections_( meanMotionTimeDependentPhaseNutationCorrections ),
        timeDependentPhaseCorrectionFunctions_( timeDependentPhaseCorrectionFunctions ),
        rotationRateCorrections_( rotationRateCorrections ),
        xPolarMotionCoefficients_( xPolarMotionCoefficients ), yPolarMotionCoefficients_( yPolarMotionCoefficients ),
        centralBody_( centralBody ){ }

    void updateAnglesAtEpoch( Eigen::Vector3d anglesAtEpoch )
    {
        anglePsiAtEpoch_ = anglesAtEpoch.x( );
        angleIAtEpoch_ = anglesAtEpoch.y( );
        anglePhiAtEpoch_ = anglesAtEpoch.z( );
    }

    double getAngleN( )
    {
        return angleN_;
    }

    double getAngleJ( )
    {
        return angleJ_;
    }

    double getAnglePsiAtEpoch( )
    {
        return anglePsiAtEpoch_;
    }
    double getAnglePsiRateAtEpoch( )
    {
        return anglePsiRateAtEpoch_;
    }

    double getAngleIAtEpoch( )
    {
        return angleIAtEpoch_;
    }

    double getAngleIRateAtEpoch( )
    {
        return angleIRateAtEpoch_;
    }

    double getAnglePhiAtEpoch( )
    {
        return anglePhiAtEpoch_;
    }

    double getAnglePhiRateAtEpoch( )
    {
        return anglePhiRateAtEpoch_;
    }

    double getCoreFactor()
    {
        return coreFactor_;
    }

    double getFreeCoreNutationRate()
    {
        return freeCoreNutationRate_;
    }

    std::map< double, std::pair< double, double > > getMeanMotionDirectNutationCorrections( )
    {
        return meanMotionDirectNutationCorrections_;
    }

    std::vector< std::map< double, std::pair< double, double > > > getMeanMotionTimeDependentPhaseNutationCorrections( )
    {
        return meanMotionTimeDependentPhaseNutationCorrections_;
    }

    std::vector< std::function< double( const double ) > > getTimeDependentPhaseCorrectionFunctions( )
    {
        return timeDependentPhaseCorrectionFunctions_;
    }

    std::map< double, std::pair< double, double > > getRotationRateCorrections( )
    {
        return rotationRateCorrections_;
    }

    std::map< double, std::pair< double, double > > getxPolarMotionCoefficients( )
    {
        return xPolarMotionCoefficients_;
    }

    std::map< double, std::pair< double, double > > getyPolarMotionCoefficients( )
    {
        return yPolarMotionCoefficients_;
    }


    std::string getCentralBody( )
    {
        return centralBody_;
    }

    void setPeriodTermsToZero( )
    {
        meanMotionDirectNutationCorrections_.clear( );
        meanMotionTimeDependentPhaseNutationCorrections_.clear( );
        timeDependentPhaseCorrectionFunctions_.clear( );
        rotationRateCorrections_.clear( );
        xPolarMotionCoefficients_.clear( );
        yPolarMotionCoefficients_.clear( );
    }

    void setPeriodTerms( std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections,
                         std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrection,
                         std::map< double, std::pair< double, double > > rotationRateCorrections,
                         std::map< double, std::pair< double, double > > xPolarMotionCoefficients,
                         std::map< double, std::pair< double, double > > yPolarMotionCoefficients)
    {
        meanMotionDirectNutationCorrections_ = meanMotionDirectNutationCorrections;
        meanMotionTimeDependentPhaseNutationCorrections_ = meanMotionTimeDependentPhaseNutationCorrection;
        rotationRateCorrections_ = rotationRateCorrections;
        xPolarMotionCoefficients_ = xPolarMotionCoefficients;
        yPolarMotionCoefficients_ = yPolarMotionCoefficients;
    }

    void setCoreFactorAndFreeCoreNutation ( double coreFactor, double freeCoreNutationRate )
    {
        coreFactor_ = coreFactor;
        freeCoreNutationRate_ = freeCoreNutationRate;
    }

private:
    double angleN_;
    double angleJ_;
    double anglePsiAtEpoch_;
    double anglePsiRateAtEpoch_;
    double angleIAtEpoch_;
    double angleIRateAtEpoch_;
    double anglePhiAtEpoch_;
    double anglePhiRateAtEpoch_;
    double coreFactor_;
    double freeCoreNutationRate_;
    std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections_;
    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections_;
    std::vector< std::function< double( const double ) > > timeDependentPhaseCorrectionFunctions_;
    std::map< double, std::pair< double, double > > rotationRateCorrections_;
    std::map< double, std::pair< double, double > > xPolarMotionCoefficients_;
    std::map< double, std::pair< double, double > > yPolarMotionCoefficients_;

    std::string centralBody_;
};


//Function to create a state function for a body, valid both during propagation, and outside propagation
/*
 * Function to create a state function for a body, valid both during propagation, and outside propagation
 * \param bodies List of body objects
 * \param orbitingBody Body for which state function is to be created
 * \param centralBody Central body w.r.t. which state function is to be created
 * \return Required state function
 */
std::function< Eigen::Vector6d( const double, bool ) > createRelativeStateFunction(
        const SystemOfBodies& bodies,
        const std::string orbitingBody,
        const std::string centralBody );

//Function to create a rotation model.
/*
 *  Function to create a rotation model based on model-specific settings for the rotation.
 *  \param rotationModelSettings Settings for the rotation model that is to be created, defined
 *  a pointer to an object of class (derived from) RotationSettings.
 *  \param body Name of the body for which the rotation model is to be created.
 * \param bodies List of body objects
 *  \return Rotation model created according to settings in rotationModelSettings.
 */


//! Function to set the angle of attack to trimmed conditions.
/*!
 * Function to set the angle of attack to trimmed conditions. Using this function requires the aerodynamic coefficient
 * interface to be dependent on the angle of attack.
 * \param bodyWithFlightConditions Body for which trimmed conditions are to be imposed.
 */
std::shared_ptr< aerodynamics::TrimOrientationCalculator > setTrimmedConditions(
        const std::shared_ptr< Body > bodyWithFlightConditions );

std::shared_ptr< ephemerides::InertialBodyFixedDirectionCalculator > createInertialDirectionCalculator(
        const std::shared_ptr< InertialDirectionSettings > directionSettings,
        const std::string& body,
        const SystemOfBodies& bodies = SystemOfBodies( ) );

std::shared_ptr< ephemerides::DirectionBasedRotationalEphemeris > createStateDirectionBasedRotationModel(
        const std::string& body,
        const std::string& centralBody,
        const SystemOfBodies& bodies,
        const Eigen::Vector3d& associatedBodyFixedDirection,
        const std::string& originalFrame,
        const std::string& targetFrame,
        const bool isColinearWithVelocity,
        const bool directionIsOppositeToVector,
        const std::function< double( const double ) > freeRotationAngleFunction = nullptr );

std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > createAerodynamicAngleBasedRotationModel(
        const std::string& body,
        const std::string& centralBody,
        const SystemOfBodies& bodies ,
        const std::string& originalFrame,
        const std::string& targetFrame );

std::shared_ptr< ephemerides::RotationalEphemeris > createTrimmedAerodynamicAngleBasedRotationModel(
        const std::string& body,
        const std::string& centralBody,
        const SystemOfBodies& bodies ,
        const std::string& originalFrame,
        const std::string& targetFrame,
        const std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction = nullptr );

std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body,
        const SystemOfBodies& bodies = SystemOfBodies( ) );

//! @get_docstring(simpleRotationModelSettings)
inline std::shared_ptr< RotationModelSettings > simpleRotationModelSettings(
        const std::string& originalFrame,
        const std::string& targetFrame,
        const Eigen::Quaterniond& initialOrientation,
        const double initialTime,
        const double rotationRate )
{
    return std::make_shared< SimpleRotationModelSettings >(
                originalFrame, targetFrame, initialOrientation, initialTime, rotationRate );
}

//! @get_docstring(simpleRotationModelSettings, 1)
inline std::shared_ptr< RotationModelSettings > simpleRotationModelSettings(
        const std::string& originalFrame,
        const std::string& targetFrame,
        const Eigen::Matrix3d& initialOrientation,
        const double initialTime,
        const double rotationRate )
{
    return std::make_shared< SimpleRotationModelSettings >(
                originalFrame, targetFrame, Eigen::Quaterniond( initialOrientation ), initialTime, rotationRate );
}

//! @get_docstring(simpleRotationModelFromSpiceSettings)
inline std::shared_ptr< RotationModelSettings > simpleRotationModelFromSpiceSettings(
        const std::string& originalFrame,
        const std::string& targetFrame,
        const std::string& targetFrameSpice,
        const double initialTime )
{
    return std::make_shared< SimpleRotationModelSettings >(
                originalFrame, targetFrame, spice_interface::computeRotationQuaternionBetweenFrames(
                    originalFrame, targetFrameSpice, initialTime ), initialTime,
                spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame(
                    originalFrame, targetFrameSpice, initialTime ).norm( ) );
}

//! @get_docstring(constantRotationModelSettings)
inline std::shared_ptr< RotationModelSettings > constantRotationModelSettings(
        const std::string& originalFrame,
        const std::string& targetFrame,
        const Eigen::Quaterniond& initialOrientation )
{
    return std::make_shared< SimpleRotationModelSettings >( originalFrame, targetFrame, initialOrientation,
                                                            0.0, 0.0 );
}

//! @get_docstring(constantRotationModelSettings, 1)
inline std::shared_ptr< RotationModelSettings > constantRotationModelSettings(
        const std::string& originalFrame,
        const std::string& targetFrame,
        const Eigen::Matrix3d& initialOrientation )
{
    return std::make_shared< SimpleRotationModelSettings >(
                originalFrame, targetFrame, Eigen::Quaterniond( initialOrientation ), 0.0, 0.0 );
}

//! @get_docstring(spiceRotationModelSettings)
inline std::shared_ptr< RotationModelSettings > spiceRotationModelSettings(
        const std::string& originalFrame,
        const std::string& targetFrame
        )
{
    return std::make_shared< RotationModelSettings >(
                spice_rotation_model, originalFrame, targetFrame );
}

//! @get_docstring(gcrsToItrsRotationModelSettings)
inline std::shared_ptr< RotationModelSettings > gcrsToItrsRotationModelSettings(
        const basic_astrodynamics::IAUConventions nutationTheory = basic_astrodynamics::iau_2006,
        const std::string baseFrameName = "GCRS" )
{
    return std::make_shared< GcrsToItrsRotationModelSettings >(
                nutationTheory, baseFrameName
                );
}

//! @get_docstring(synchronousRotationModelSettings)
inline std::shared_ptr< RotationModelSettings > synchronousRotationModelSettings(
        const std::string& centralBodyName,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation )
{
    return std::make_shared< SynchronousRotationModelSettings >(
                centralBodyName, baseFrameOrientation, targetFrameOrientation );
}

inline std::shared_ptr< RotationModelSettings > tabulatedRotationSettings(
        const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{
    return std::make_shared< TabulatedRotationSettings >(
                rotationalStateHistory, baseFrameOrientation, targetFrameOrientation, interpolatorSettings );
}

inline std::shared_ptr< RotationModelSettings > aerodynamicAngleRotationSettings(
        const std::string& centralBody,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction = nullptr )
{
    return std::make_shared< AerodynamicAngleRotationSettings >(
                centralBody, baseFrameOrientation, targetFrameOrientation, aerodynamicAngleFunction );
}

inline std::shared_ptr< RotationModelSettings > pitchTrimRotationSettings(
        const std::string& centralBody,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const std::function< Eigen::Vector2d( const double ) > aerodynamicAngleFunction = nullptr )
{
    return std::make_shared< PitchTrimRotationSettings >(
                centralBody, baseFrameOrientation, targetFrameOrientation, aerodynamicAngleFunction );
}


inline std::shared_ptr< RotationModelSettings > bodyFixedDirectionBasedRotationSettings(
        const std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const std::function< double( const double ) > freeRotationAngleFunction = nullptr )
{
    return std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                inertialBodyAxisDirectionFunction, baseFrameOrientation, targetFrameOrientation, freeRotationAngleFunction );
}

inline std::shared_ptr< RotationModelSettings > orbitalStateBasedRotationSettings(
        const std::string& centralBody,
        const bool isColinearWithVelocity,
        const bool directionIsOppositeToVector,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const std::function< double( const double ) > freeRotationAngleFunction = nullptr )
{
    return std::make_shared< OrbitalStateBasedRotationSettings >(
                centralBody, isColinearWithVelocity, directionIsOppositeToVector, baseFrameOrientation, targetFrameOrientation, freeRotationAngleFunction );
}


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEROTATIONMODEL_H
