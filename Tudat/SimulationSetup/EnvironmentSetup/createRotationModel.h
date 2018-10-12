/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"

namespace tudat
{

namespace simulation_setup
{

//! List of rotation models available in simulations
/*!
 *  List of rotation models available in simulations. Rotation models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum RotationModelType
{
    simple_rotation_model,
    spice_rotation_model,
    gcrs_to_itrs_rotation_model
};

//! Class for providing settings for rotation model.
/*!
 *  Class for providing settings for automatic rotation model creation. This class is a
 *  functional (base) class for settings of rotation models that require no information in
 *  addition to their type. Rotation model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
class RotationModelSettings
{
public:

    //! Constructor, sets type of rotation model.
    /*!
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

    //! Destructor.
    virtual ~RotationModelSettings( ){ }

    //! Function to return the type of rotation model that is to be created.
    /*!
     *  Function to return the type of rotation model that is to be created.
     *  \return Type of rotation model that is to be created.
     */
    RotationModelType getRotationType( ){ return rotationType_; }

    //! Function to return the base frame of rotation model.
    /*!
     *  Function to return the base frame of rotation model.
     *  \return Base frame of rotation model.
     */
    std::string getOriginalFrame( ){ return originalFrame_; }

    //! Function to return the target frame of rotation model.
    /*!
     *  Function to return the target frame of rotation model.
     *  \return Target frame of rotation model.
     */
    std::string getTargetFrame( ){ return targetFrame_; }


    //! Function to rese the orientation of the base frame.
    /*!
     * Function to reset the orientation of the base frame.
     * \param originalFrame New base frame orientation
     */
    void resetOriginalFrame( const std::string& originalFrame )
    {
        originalFrame_ = originalFrame;
    }

protected:

    //! Type of rotation model that is to be created.
    RotationModelType rotationType_;

    //! Target frame of rotation model.
    std::string originalFrame_;

    //! Base frame of rotation model.
    std::string targetFrame_;

};

//! RotationModelSettings derived class for defining settings of a simple rotational ephemeris.
class SimpleRotationModelSettings: public RotationModelSettings
{
public:
    //! Constructor,
    /*!
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

    //! Function to return rotation from base to target frame at initialTime.
    /*!
     *  Function to return rotation from base to target frame at initialTime.
     *  \return Rotation from base to target frame at initialTime.
     */
    Eigen::Quaterniond getInitialOrientation( ){ return initialOrientation_; }

    //! Function to return time at which initialOrientation represents the instantaneous rotation.
    /*!
     *  Function to return time at which initialOrientation represents the instantaneous rotation.
     *  \return Time at which initialOrientation represents the instantaneous rotation.
     */
    double getInitialTime( ){ return initialTime_; }

    //! Function to return rotation rate of body about its local z-axis.
    /*!
     *  Function to return rotation rate of body about its local z-axis.
     *  \return Rotation rate of body about its local z-axis.
     */
    double getRotationRate( ){ return rotationRate_; }

private:

    //!  Rotation from base to target frame at initialTime.
    Eigen::Quaterniond initialOrientation_;

    //! Time at which initialOrientation represents the instantaneous rotation.
    double initialTime_;

    //! Rotation rate of body about its local z-axis.
    double rotationRate_;
};

#if USE_SOFA

//! Struct that holds settings for EOP short-period variation
struct EopCorrectionSettings
{
    //! Constructor
    /*!
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

    //! Conversion factor to be used for amplitudes to multiply input values
    double conversionFactor_;

    //! Minimum amplitude that is read from files and considered in calculations.
    double minimumAmplitude_;

    //! List of files with amplitudes for corrections
    std::vector< std::string > amplitudesFiles_;

    //! Fundamental argument multiplier for corrections
    std::vector< std::string > argumentMultipliersFile_;
};

//! Settings for creating a GCRS<->ITRS rotation model
class GcrsToItrsRotationModelSettings: public RotationModelSettings
{
public:

    //! Constructor
    /*!
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
            const std::string& eopFile = input_output::getEarthOrientationDataFilesPath( ) + "eopc04_08_IAU2000.62-now.txt",
            const basic_astrodynamics::TimeScales inputTimeScale = basic_astrodynamics::tdb_scale,
            const std::shared_ptr< EopCorrectionSettings > ut1CorrectionSettings =
            std::make_shared< EopCorrectionSettings >(
                1.0E-6, 0.0, std::vector< std::string >{
                    input_output::getEarthOrientationDataFilesPath( ) + "utcLibrationAmplitudes.txt",
                    input_output::getEarthOrientationDataFilesPath( ) + "utcOceanTidesAmplitudes.txt" },
                std::vector< std::string >{
                    input_output::getEarthOrientationDataFilesPath( ) +
                    "utcLibrationFundamentalArgumentMultipliers.txt",
                    input_output::getEarthOrientationDataFilesPath( ) +
                    "utcOceanTidesFundamentalArgumentMultipliers.txt" } ),
            const std::shared_ptr< EopCorrectionSettings > polarMotionCorrectionSettings =
            std::make_shared< EopCorrectionSettings >(
                unit_conversions::convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0, std::vector< std::string >{
                    input_output::getEarthOrientationDataFilesPath( ) +
                    "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
                    input_output::getEarthOrientationDataFilesPath( ) +
                    "polarMotionOceanTidesAmplitudes.txt", },
                std::vector< std::string >{
                    input_output::getEarthOrientationDataFilesPath( ) +
                    "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
                    input_output::getEarthOrientationDataFilesPath( ) +
                    "polarMotionOceanTidesFundamentalArgumentMultipliers.txt" } ) ):
        RotationModelSettings( gcrs_to_itrs_rotation_model, baseFrameName, "ITRS" ),
        inputTimeScale_( inputTimeScale ), nutationTheory_( nutationTheory ), eopFile_( eopFile ),
        eopFileFormat_( "C04" ), ut1CorrectionSettings_( ut1CorrectionSettings ),
        polarMotionCorrectionSettings_( polarMotionCorrectionSettings ){ }

    //! Destructor
    ~GcrsToItrsRotationModelSettings( ){ }

    //! Function to retrieve the time scale in which input to the rotation model class is provided
    /*!
     * Function to retrieve the time scale in which input to the rotation model class is provided
     * \return Time scale in which input to the rotation model class is provided
     */
    basic_astrodynamics::TimeScales getInputTimeScale( )
    {
        return inputTimeScale_;
    }
    //! Function to retrieve the IAU precession-nutation theory that is to be used
    /*!
     * Function to retrieve the IAU precession-nutation theory that is to be used
     * \return IAU precession-nutation theory that is to be used
     */
    basic_astrodynamics::IAUConventions getNutationTheory( )
    {
        return nutationTheory_;
    }

    //! Function to retrieve the name of EOP file that is to be used
    /*!
     * Function to retrieve the name of EOP file that is to be used
     * \return Name of EOP file that is to be used
     */
    std::string getEopFile( )
    {
        return eopFile_;
    }

    //! Function to retrieve the identifier for file format that is provided
    /*!
     * Function to retrieve the identifier for file format that is provided
     * \return Identifier for file format that is provided
     */
    std::string getEopFileFormat( )
    {
        return eopFileFormat_;
    }

    //! Function to retrieve the settings for short-period UT1-UTC variations
    /*!
     * Function to retrieve the settings for short-period UT1-UTC variations
     * \return Settings for short-period UT1-UTC variations
     */
    std::shared_ptr< EopCorrectionSettings > getUt1CorrectionSettings( )
    {
        return ut1CorrectionSettings_;
    }

    //! Function to retrieve the settings for short-period polar motion variations
    /*!
     * Function to retrieve the Settings for short-period polar motion variations
     * \return settings for short-period polar motion variations
     */
    std::shared_ptr< EopCorrectionSettings > getPolarMotionCorrectionSettings( )
    {
        return polarMotionCorrectionSettings_;
    }

private:

    //! Time scale in which input to the rotation model class is provided
    basic_astrodynamics::TimeScales inputTimeScale_;

    //! IAU precession-nutation theory that is to be used
    basic_astrodynamics::IAUConventions nutationTheory_;

    //! Name of EOP file that is to be used
    std::string eopFile_;

    //! Identifier for file format that is provided
    std::string eopFileFormat_;

    //! Settings for short-period UT1-UTC variations
    std::shared_ptr< EopCorrectionSettings > ut1CorrectionSettings_;

    //! Settings for short-period polar motion variations
    std::shared_ptr< EopCorrectionSettings > polarMotionCorrectionSettings_;

};
#endif

//! Function to create a rotation model.
/*!
 *  Function to create a rotation model based on model-specific settings for the rotation.
 *  \param rotationModelSettings Settings for the rotation model that is to be created, defined
 *  a pointer to an object of class (derived from) RotationSettings.
 *  \param body Name of the body for which the rotation model is to be created.
 *  \return Rotation model created according to settings in rotationModelSettings.
 */
std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEROTATIONMODEL_H
