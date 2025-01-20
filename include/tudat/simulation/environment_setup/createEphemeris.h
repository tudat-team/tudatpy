/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEEPHEMERIS_H
#define TUDAT_CREATEEPHEMERIS_H

#include <string>
#include <map>

#include <memory>

#include "tudat/io/matrixTextFileReader.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositionsCircularCoplanar.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/spice/spiceEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

// List of ephemeris models available in simulations
/*
 *  List of ephemeris models available in simulations. Ephemeris models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
//! @get_docstring(EphemerisType.__docstring__)
enum EphemerisType
{
    approximate_planet_positions,
    direct_spice_ephemeris,
    tabulated_ephemeris,
    auto_generated_tabulated_ephemeris,
    interpolated_spice,
    constant_ephemeris,
    kepler_ephemeris,
    custom_ephemeris,
    direct_tle_ephemeris,
    interpolated_tle_ephemeris,
    scaled_ephemeris
};

// Class for providing settings for ephemeris model.
/*
 *  Class for providing settings for automatic ephemeris model creation. This class is a
 *  functional (base) class for settings of ephemeris models that require no information in
 *  addition to their type (and frame origin and orientation). Ephemeris model classes defining
 *  requiring additional information must be created using an object derived from this class.
 */

//! @get_docstring(EphemerisSettings.__docstring__)
class EphemerisSettings
{
public:

    // Constructor, sets type of ephemeris model.
    /*
     *  Constructor, sets type of ephemeris model and frame origin and orientation.
     *  Settings for ephemeris models requiring additional information should be defined in a
     *  derived class.
     *  \param ephemerisType Type of ephemeris model that is to be created.
     *  \param frameOrigin Origin of frame in which ephemeris data is defined
     *         (optional "SSB" by default).
     *  \param frameOrientation Orientation of frame in which ephemeris data is defined
     *         (optional "ECLIPJ2000" by default).
     */
    EphemerisSettings( const EphemerisType ephemerisType,
                       const std::string& frameOrigin = "SSB",
                       const  std::string& frameOrientation = "ECLIPJ2000" ):
        ephemerisType_( ephemerisType ),
        frameOrigin_( frameOrigin ),
        frameOrientation_( frameOrientation ),
        makeMultiArcEphemeris_( false ){ }

    // Destructor
    virtual ~EphemerisSettings( ){ }

    // Function to return type of ephemeris that is to be created.
    /*
     *  Function to return type of ephemeris that is to be created.
     *  \return Type of ephemeris that is to be created.
     */
    EphemerisType getEphemerisType( ){ return ephemerisType_; }

    // Function to return the origin of frame in which ephemeris data is defined.
    /*
     *  Function to return the origin of frame in which ephemeris data is defined.
     *  \return Origin of frame in which ephemeris data is defined.
     */
    std::string getFrameOrigin( ){ return frameOrigin_; }

    // Function to return the orientation of frame in which ephemeris data is defined.
    /*
     *  Function to return the orientation of frame in which ephemeris data is defined.
     *  \return Orientation of frame in which ephemeris data is defined.
     */
    std::string getFrameOrientation( ){ return frameOrientation_;}

    // Function to retrieve boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
    /*
     * Function to retrieve boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
     * \return Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
     */
    bool getMakeMultiArcEphemeris( )
    {
        return makeMultiArcEphemeris_;
    }

    // Function to reset the origin of the frame.
    /*
     * Function to reset the origin of the frame.
     * \param frameOrigin New frame origin
     */
    void resetFrameOrigin( const std::string& frameOrigin ){ frameOrigin_ = frameOrigin; }

    // Function to rese the orientation of the frame.
    /*
     * Function to reset the orientation of the frame.
     * \param frameOrientation New frame orientation
     */
    void resetFrameOrientation( const std::string& frameOrientation ){ frameOrientation_ = frameOrientation; }

    // Function to reset boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
    /*
     * Function to reset boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
     * \param makeMultiArcEphemeris New boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
     */
    void resetMakeMultiArcEphemeris( const bool makeMultiArcEphemeris )
    {
        makeMultiArcEphemeris_ = makeMultiArcEphemeris;
    }

protected:

    // Type of ephemeris model that is to be created.
    EphemerisType ephemerisType_;

    // Origin of frame in which ephemeris data is defined.
    std::string frameOrigin_;

    // Orientation of frame in which ephemeris data is defined.
    std::string frameOrientation_;

    // Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris
    /*
     *  Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris. If true, the createEphemeris
     *  function creates a multi-arc ephemeris with a single arc spanning all time, created according to the contents of the
     *  EphemerisSettings object.
     */
    bool makeMultiArcEphemeris_;
};


    //! @get_docstring(ScaledEphemerisSettings.__docstring__)
class ScaledEphemerisSettings: public EphemerisSettings
{
public:

    ScaledEphemerisSettings(
            const std::shared_ptr< EphemerisSettings > baseSettings,
            const double scaling,
            const bool isScalingAbsolute ):
        EphemerisSettings( scaled_ephemeris, baseSettings->getFrameOrigin( ), baseSettings->getFrameOrientation( ) ),
        baseSettings_( baseSettings ), scaling_( [=]( const double ){ return Eigen::Vector6d::Constant( scaling ); } ), isScalingAbsolute_( isScalingAbsolute ){ }


    ScaledEphemerisSettings(
            const std::shared_ptr< EphemerisSettings > baseSettings,
            const Eigen::Vector6d scaling,
            const bool isScalingAbsolute ):
        EphemerisSettings( scaled_ephemeris, baseSettings->getFrameOrigin( ), baseSettings->getFrameOrientation( ) ),
        baseSettings_( baseSettings ), scaling_( [=]( const double ){ return scaling; } ), isScalingAbsolute_( isScalingAbsolute ){ }

    ScaledEphemerisSettings(
            const std::shared_ptr< EphemerisSettings > baseSettings,
            const std::function< Eigen::Vector6d( const double ) > scaling,
            const bool isScalingAbsolute ):
        EphemerisSettings( scaled_ephemeris, baseSettings->getFrameOrigin( ), baseSettings->getFrameOrientation( ) ),
    baseSettings_( baseSettings ), scaling_( scaling ), isScalingAbsolute_( isScalingAbsolute ){ }

    std::shared_ptr< EphemerisSettings > getBaseSettings( )
    {
        return baseSettings_;
    }

    std::function< Eigen::Vector6d( const double ) > getScaling( )
    {
        return scaling_;
    }

    bool getIsScalingAbsolute( )
    {
        return isScalingAbsolute_;
    }

protected:

    std::shared_ptr< EphemerisSettings > baseSettings_;

    std::function< Eigen::Vector6d( const double ) > scaling_;

    bool isScalingAbsolute_;
};

// EphemerisSettings derived class for defining settings of an ephemeris linked directly to Spice.
//! @get_docstring(DirectSpiceEphemerisSettings.__docstring__)
class DirectSpiceEphemerisSettings: public EphemerisSettings
{
public:

    // Constructor.
    /* Constructor, sets the properties from which the Spice ephemeris is to be retrieved.
     * \param frameOrigin Name of body relative to which the ephemeris is to be calculated
     *        (optional "SSB" by default).
     * \param correctForStellarAberration Boolean whether to correct for stellar Aberration in
     *          retrieved values of (observed state) (optional "ECLIPJ2000" by defalut).
     * \param correctForLightTimeAberration Boolean whether to correct for light time in
     *          retrieved values of (observed state) (optional false by default).
     * \param convergeLighTimeAberration Boolean whether to use single iteration or max. 3
     *          iterations for calculating light time (optional false by default).
     * \param frameOrientation Orientatioan of the reference frame in which the epehemeris is to be
     *          calculated (optional false by default).
     * \param ephemerisType Type of ephemeris that is to be created, always set to
     *          direct_spice_ephemeris when using this class directly. The derived class
     *          InterpolatedSpiceEphemerisSettings sets this parameter to a different value. Not
     *          to be changed by used (optional direct_spice_ephemeris by default).
     */
    DirectSpiceEphemerisSettings( const std::string frameOrigin = "SSB",
                                  const std::string frameOrientation = "ECLIPJ2000",
                                  const bool correctForStellarAberration = false,
                                  const bool correctForLightTimeAberration = false,
                                  const bool convergeLighTimeAberration = false,
                                  const EphemerisType ephemerisType = direct_spice_ephemeris ):
        EphemerisSettings( ephemerisType, frameOrigin, frameOrientation ),
        correctForStellarAberration_( correctForStellarAberration ),
        correctForLightTimeAberration_( correctForLightTimeAberration ),
        convergeLighTimeAberration_( convergeLighTimeAberration ),
        bodyNameOverride_( "" ){ }

    DirectSpiceEphemerisSettings( const std::string frameOrigin,
                                  const std::string frameOrientation,
                                  const std::string bodyNameOverride,
                                  const EphemerisType ephemerisType = direct_spice_ephemeris ):
        EphemerisSettings( ephemerisType, frameOrigin, frameOrientation ),
        correctForStellarAberration_( false ),
        correctForLightTimeAberration_( false ),
        convergeLighTimeAberration_( false ),
        bodyNameOverride_( bodyNameOverride ){ }


    // Destructor
    virtual ~DirectSpiceEphemerisSettings( ){ }

    // Returns whether to correct for stellar aberration in retrieved values of (observed state).
    /*
     *  Returns whether to correct for stellar aberration in retrieved values of (observed state).
     *  \return Boolean defining whether to correct for stellar aberration in retrieved
     *  values of (observed state).
     */
    bool getCorrectForStellarAberration( ){ return correctForStellarAberration_; }

    // Returns whether to correct for light time in retrieved values of (observed state).
    /*
     *  Returns whether to correct for light time in retrieved values of (observed state).
     *  \return Boolean defining whether to correct for light time in retrieved values of
     *  (observed state).
     */
    bool getCorrectForLightTimeAberration( ){ return correctForLightTimeAberration_; }

    // Returns whether to use single iteration or max. 3 iterations for calculating light time.
    /*
     *  Returns whether to use single iteration or max. 3 iterations for calculating light time.
     *  \return Boolean defining whether to use single iteration or max. 3 iterations for
     *  calculating light time.
     */
    bool getConvergeLighTimeAberration( ){ return convergeLighTimeAberration_; }

    std::string getBodyNameOverride( ){ return bodyNameOverride_; }

protected:

    // Boolean whether to correct for stellar aberration in retrieved values of (observed state).
    bool correctForStellarAberration_;

    // Boolean whether to correct for light time in retrieved values of (observed state).
    bool correctForLightTimeAberration_;

    // Boolean whether to use single iteration or max. 3 iterations for calculating light time.
    bool convergeLighTimeAberration_;

    std::string bodyNameOverride_;
};

// EphemerisSettings derived class for defining settings of a ephemeris interpolated from Spice
// data.
/*
 *  EphemerisSettings derived class for defining settings of a ephemeris interpolated from Spice
 *  data. Retrieving a state from Spice can be very computationally intensive. Using the
 *  settings of this class, ephemeris data for a body is pre-computed using a limited number
 *  of calls to Spice, which is then used to create an interpolator (6th order Lagrange). For
 *  many numerical integration scenarios, this approach may be faster than using
 *  DirectSpiceEphemerisSettings, with negligible influence on accuracy.
 */

//! @get_docstring(InterpolatedSpiceEphemerisSettings.__docstring__)
class InterpolatedSpiceEphemerisSettings: public DirectSpiceEphemerisSettings
{
public:

    // Constructor.
    /* Constructor, sets the properties from which the tabulated spice data is to be created
     *  from which an ephemeris is to be created.
     * \param initialTime Initial time from which interpolated data from Spice should be created.
     * \param finalTime Final time from which interpolated data from Spice should be created.
     * \param timeStep Time step with which interpolated data from Spice should be created.
     * \param frameOrigin Name of body relative to which the ephemeris is to be calculated
     *        (optional "SSB" by default).
     * \param frameOrientation Orientatioan of the reference frame in which the epehemeris is to be
     *          calculated (optional "ECLIPJ2000" by default).
     * \param interpolatorSettings Settings to be used for the state interpolation.
     */
    InterpolatedSpiceEphemerisSettings( double initialTime,
                                        double finalTime,
                                        double timeStep,
                                        std::string frameOrigin = "SSB",
                                        std::string frameOrientation = "ECLIPJ2000",
                                        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ),
                                        const std::string bodyNameOverride = "" ):
        DirectSpiceEphemerisSettings( frameOrigin, frameOrientation, bodyNameOverride,
                                      interpolated_spice ),
        initialTime_( initialTime ), finalTime_( finalTime ), timeStep_( timeStep ),
        interpolatorSettings_( interpolatorSettings ){ }

    // Function to return initial time from which interpolated data from Spice should be created.
    /*
     *  Function to return initial time from which interpolated data from Spice should be created.
     *  \return Initial time from which interpolated data from Spice should be created.
     */
    double getInitialTime( ){ return initialTime_; }

    // Function to return final time from which interpolated data from Spice should be created.
    /*
     *  Function to return final time from which interpolated data from Spice should be created.
     *  \return Final time from which interpolated data from Spice should be created.
     */
    double getFinalTime( ){ return finalTime_; }

    // Function to return time step with which interpolated data from Spice should be created.
    /*
     *  Function to return time step with which interpolated data from Spice should be created.
     *  \return Time step with which interpolated data from Spice should be created.
     */
    double getTimeStep( ){ return timeStep_; }

    // Function to return settings to be used for the state interpolation.
    /*
     *  Function to return settings to be used for the state interpolation.
     *  \return Settings to be used for the state interpolation.
     */

    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( )
    {
        return interpolatorSettings_;
    }


private:

    // Initial time from which interpolated data from Spice should be created.
    double initialTime_;

    // Final time until which interpolated data from Spice should be created.
    double finalTime_;

    // Time step with which interpolated data from Spice should be created.
    double timeStep_;

    // Settings to be used for the state interpolation.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;
};

// EphemerisSettings derived class for defining settings of an approximate ephemeris for major
// planets.
/*
 *  EphemerisSettings derived class for defining settings of an approximate ephemeris for major
 *  planets, as inplemented in ApproximateJplEphemeris class and derived class,
 *  described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf.
 */
//! @get_docstring(ApproximateJplEphemerisSettings.__docstring__)
class ApproximateJplEphemerisSettings: public EphemerisSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyIdentifier Parameter identifying for which body an ephemeris is to be created.
     *  \param useCircularCoplanarApproximation Boolean defining whether a circular, coplanar
     *  orbit of the body is to be assumed, or whether a non-zero inclination and long-period
     *  changes in the orbit are to be included.
     */
    ApproximateJplEphemerisSettings(
            const std::string bodyName,
            const bool useCircularCoplanarApproximation = false ):
        EphemerisSettings( approximate_planet_positions ),
        bodyName_( bodyName ),
        useCircularCoplanarApproximation_( useCircularCoplanarApproximation ){ }

    // Function to return parameter identifying for which body an ephemeris is to be created.
    /*
     *  Function to return parameter identifying for which body an ephemeris is to be created.
     *  \return Parameter identifying for which body an ephemeris is to be created.
     */
    std::string getBodyName( )
    {
        return bodyName_;
    }

    // Function to return whether a circular, coplanar orbit of the body is to be assumed.
    /*
     *  Function to return whether a circular, coplanar orbit of the body is to be assumed.
     *  \return Boolean defining whether a circular, coplanar orbit of the body is to be assumed.
     */
    bool getUseCircularCoplanarApproximation( )
    {
        return useCircularCoplanarApproximation_;
    }

private:

    // Parameter identifying for which body an ephemeris is to be created.
    std::string bodyName_;

    //  Boolean defining whether a circular, coplanar orbit of the body is to be assumed,
    /*
     *  Boolean defining whether a circular, coplanar orbit of the body is to be assumed
     *  (creating an ApproximateJplCircularCoplanarEphemeris object), or whether
     *  a non-zero inclination and long-period changes in the orbit are to be included
     *  (creating an ApproximateJplEphemeris object).
     */
     bool useCircularCoplanarApproximation_;
};


// EphemerisSettings derived class for defining settings of an ephemeris producing a constant
// (time-independent) state

//! @get_docstring(ConstantEphemerisSettings.__docstring__)
class ConstantEphemerisSettings: public EphemerisSettings
{
public:

    // Constructor of settings for an ephemeris producing a constant (time-independent) state.
    /*
     * Constructor of settings for an ephemeris producing a constant (time-independent) state.
     * \param constantState Constant state that will be provided as output of the ephemeris at all times.
     * \param frameOrigin Origin of frame in which ephemeris data is defined.
     * \param frameOrientation Orientation of frame in which ephemeris data is defined.
     */
    ConstantEphemerisSettings( const Eigen::Vector6d& constantState,
                               const std::string& frameOrigin = "SSB",
                               const std::string& frameOrientation = "ECLIPJ2000" ):
        EphemerisSettings( constant_ephemeris,
                           frameOrigin,
                           frameOrientation ), constantState_( constantState ){ }


    Eigen::Vector6d getConstantState( ){ return constantState_; }

private:

    // Constant state that will be provided as output of the ephemeris at all times.
    Eigen::Vector6d constantState_;
};

// EphemerisSettings derived class for defining settings of an ephemeris producing a custom
// state (e.g. arbitrary state as a function of time)
class CustomEphemerisSettings: public EphemerisSettings
{
public:

    // Constructor of settings for an ephemeris producing a constant (time-independent) state.
    /*
     * Constructor of settings for an ephemeris producing a constant (time-independent) state.
     * \param customStateFunction Function returning the state as a function of time
     * \param frameOrigin Origin of frame in which ephemeris data is defined.
     * \param frameOrientation Orientation of frame in which ephemeris data is defined.
     */
    CustomEphemerisSettings( const std::function< Eigen::Vector6d( const double ) > customStateFunction,
                               const std::string& frameOrigin = "SSB",
                               const std::string& frameOrientation = "ECLIPJ2000" ):
        EphemerisSettings( custom_ephemeris,
                           frameOrigin,
                           frameOrientation ), customStateFunction_( customStateFunction ){ }

    // Function to return the function returning the state as a function of time
    /*
     *  Function to return the function returning the state as a function of time
     *  \return  Function returning the state as a function of time
     */
    std::function< Eigen::Vector6d( const double ) > getCustomStateFunction( )
    {
        return customStateFunction_;
    }

private:

    // Function returning the state as a function of time
    std::function< Eigen::Vector6d( const double ) > customStateFunction_;
};

// EphemerisSettings derived class for defining settings of an ephemeris representing an ideal
// Kepler orbit.

//! @get_docstring(KeplerEphemerisSettings.__docstring__)
class KeplerEphemerisSettings: public EphemerisSettings
{
public:
    // Constructor
    /*
    *  Constructor
    *  \param initialStateInKeplerianElements Kepler elements at time epochOfInitialState.
    *  \param epochOfInitialState Time at which initialStateInKeplerianElements represents
    *  the Keplerian state.
    *  \param centralBodyGravitationalParameter Gravitational parameter of the central body
    *  that is used in the computations.
    *  \param referenceFrameOrigin Origin of reference frame (string identifier).
    *  \param referenceFrameOrientation Orientation of reference frame (string identifier)
    *  \param rootFinderAbsoluteTolerance Convergence tolerance for root finder used to
    *  convert mean to eccentric anomaly on each call to getCartesianState.
    *  \param rootFinderMaximumNumberOfIterations Maximum iteration for root finder used to
    *  convert mean to eccentric anomaly on each call to getCartesianState.
    */
    KeplerEphemerisSettings( const Eigen::Vector6d& initialStateInKeplerianElements,
                             const double epochOfInitialState,
                             const double centralBodyGravitationalParameter,
                             const std::string& referenceFrameOrigin = "SSB",
                             const std::string& referenceFrameOrientation = "ECLIPJ2000",
                             const double rootFinderAbsoluteTolerance =
            200.0 * std::numeric_limits< double >::epsilon( ),
                             const double rootFinderMaximumNumberOfIterations = 1000.0 ):
        EphemerisSettings( kepler_ephemeris, referenceFrameOrigin, referenceFrameOrientation ),
        initialStateInKeplerianElements_( initialStateInKeplerianElements ),
        epochOfInitialState_( epochOfInitialState ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        rootFinderAbsoluteTolerance_( rootFinderAbsoluteTolerance ),
        rootFinderMaximumNumberOfIterations_( rootFinderMaximumNumberOfIterations ){ }


    // Function to return the kepler elements at time epochOfInitialState.
    /*
     *  Function to return the kepler elements at time epochOfInitialState.
     *  \return Kepler elements at time epochOfInitialState.
     */
    Eigen::Vector6d getInitialStateInKeplerianElements( )
    {
        return initialStateInKeplerianElements_;
    }

    // Function to return the initial epoch from which propagation of Kepler orbit is performed.
    /*
     *  Function to return the initial epoch from which propagation of Kepler orbit is performed.
     *  \return  Initial epoch from which propagation of Kepler orbit is performed.
     */
    double getEpochOfInitialState( )
    {
        return epochOfInitialState_;
    }

    // Function to return the gravitational parameter of central body about which the Kepler orbit
    // is defined.
    /*
     *  Function to return the gravitational parameter of central body about which the Kepler orbit
     *  is defined.
     *  \return Gravitational parameter of central body about which the Kepler orbit is defined.
     */
    double getCentralBodyGravitationalParameter( )
    {
        return centralBodyGravitationalParameter_;
    }

    // Function to return the convergence tolerance for root finder used to convert mean to
    // eccentric anomaly
    /*
     *  Function to return the convergence tolerance for root finder used to convert mean to
     *  eccentric anomaly
     *  \return Convergence tolerance for root finder used to convert mean to eccentric anomaly
     */
    double getRootFinderAbsoluteTolerance( )
    {
        return rootFinderAbsoluteTolerance_;
    }

    // Function to return the maximum iteration for root finder used to convert mean to eccentric
    // anomaly
    /*
     *  Function to return the maximum iteration for root finder used to convert mean to eccentric
     *  anomaly
     *  \return Maximum iteration for root finder used to convert mean to eccentric anomaly
     */
    double getRootFinderMaximumNumberOfIterations( )
    {
        return rootFinderMaximumNumberOfIterations_;
    }

private:

    // Kepler elements at time epochOfInitialState.
    Eigen::Vector6d initialStateInKeplerianElements_;

    // Initial epoch from which propagation of Kepler orbit is performed.
    double epochOfInitialState_;

    // Gravitational parameter of central body about which the Kepler orbit is defined.
    double centralBodyGravitationalParameter_;

    // Convergence tolerance for root finder used to convert mean to eccentric anomaly.
    double rootFinderAbsoluteTolerance_;

    // Maximum iteration for root finder used to convert mean to eccentric anomaly
    double rootFinderMaximumNumberOfIterations_;
};

// EphemerisSettings derived class for defining settings of an ephemeris created from tabulated
// data.
/*
 *  EphemerisSettings derived class for defining settings of an ephemeris created from tabulated
 *  data. Currently the use of an 6th order Lagrange interpolator is hardcoded, which is created
 *  from the data that is provided. Note that at the edges of the interpolation interval, a
 *  Cubic spline interpolator is used to suppres the influence of Runge's phenomenon.
 */

//! @get_docstring(TabulatedEphemerisSettings.__docstring__)
class TabulatedEphemerisSettings: public EphemerisSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyStateHistory Data map (time as key, Cartesian state as values) defining data
     *  from which an interpolated ephemeris is to be created.
     * \param frameOrigin Name of body relative to which the ephemeris is to be calculated
     *        (optional "SSB" by default).
     * \param frameOrientation Orientatioan of the reference frame in which the epehemeris is to be
     *          calculated (optional, "ECLIPJ2000" by default).
     */
    TabulatedEphemerisSettings(
            const std::map< double, Eigen::Vector6d >& bodyStateHistory,
            std::string frameOrigin = "SSB",
            std::string frameOrientation = "ECLIPJ2000" ):
        EphemerisSettings( tabulated_ephemeris, frameOrigin, frameOrientation ),
        bodyStateHistory_( bodyStateHistory ){ }

    // Function returning data map defining discrete data from which an ephemeris is to be created.
    /*
     *  Function returning data map defining discrete data from which an ephemeris is to be created.
     *  \return Data map defining discrete data from which an ephemeris is to be created.
     */
    std::map< double, Eigen::Vector6d > getBodyStateHistory( )
    { return bodyStateHistory_; }

private:

    // Data map defining discrete data from which an ephemeris is to be created.
    /*
     *  Data map (time as key, Cartesian state as values) defining data from which an interpolated
     *  ephemeris is to be created.
     */
    std::map< double, Eigen::Vector6d > bodyStateHistory_;
};

class AutoGeneratedTabulatedEphemerisSettings: public EphemerisSettings
{
public:

    AutoGeneratedTabulatedEphemerisSettings(
            const std::shared_ptr< EphemerisSettings > ephemerisSettings,
            const double startTime,
            const double endTime,
            const double timeStep,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ):
        EphemerisSettings( auto_generated_tabulated_ephemeris, ephemerisSettings->getFrameOrigin( ), ephemerisSettings->getFrameOrientation( ) ),
        ephemerisSettings_( ephemerisSettings ),
        startTime_( startTime ), endTime_( endTime ), timeStep_( timeStep ), interpolatorSettings_( interpolatorSettings ){ }

    std::shared_ptr< EphemerisSettings > getEphemerisSettings( )
    {
        return ephemerisSettings_;
    }

    double getStartTime( )
    {
        return startTime_;
    }

    double getEndTime( )
    {
        return endTime_;
    }

    double getTimeStep( )
    {
        return timeStep_;
    }

    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( )
    {
        return interpolatorSettings_;
    }

private:

    std::shared_ptr< EphemerisSettings > ephemerisSettings_;

    double startTime_;

    double endTime_;

    double timeStep_;

    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;

};


class DirectTleEphemerisSettings: public EphemerisSettings
{
public:

	DirectTleEphemerisSettings( std::shared_ptr< ephemerides::Tle > tle, const std::string frameOrigin = "Earth",
			const std::string frameOrientation = "J2000" ):
		EphemerisSettings( direct_tle_ephemeris, frameOrigin, frameOrientation ){ }

	const std::shared_ptr<ephemerides::Tle> getTle( ) const
	{
		return tle_;
	}

private:

	std::shared_ptr< ephemerides::Tle > tle_;

};

class InterpolatedTleEphemerisSettings : public EphemerisSettings
{
public:

	InterpolatedTleEphemerisSettings( const double initialTime, const double finalTime,
			const double timeStep, std::shared_ptr< ephemerides::Tle > tle,
			std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
					std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ),
			const std::string& frameOrigin = "Earth", const std::string& frameOrientation = "J2000") :
		EphemerisSettings( interpolated_tle_ephemeris, frameOrigin, frameOrientation ),
		initialTime_( initialTime ), finalTime_( finalTime ), timeStep_( timeStep ),
        interpolatorSettings_( interpolatorSettings ),
		tle_( tle ){ }

	double getInitialTime( ) const
	{
		return initialTime_;
	}

	double getFinalTime( ) const
	{
		return finalTime_;
	}

	double getTimeStep( ) const
	{
		return timeStep_;
	}

	const std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( ) const
	{
		return interpolatorSettings_;
	}

	const std::shared_ptr<ephemerides::Tle> getTle( ) const
	{
		return tle_;
	}

private:

	// Initial time from which interpolated data from TLE should be created.
	double initialTime_;

	// Final time until which interpolated data from TLE should be created.
	double finalTime_;

	// Time step with which interpolated data from TLE should be created.
	double timeStep_;

	// Settings to be used for the state interpolation.
	std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;

	std::shared_ptr< ephemerides::Tle > tle_;
};


// Function to create a tabulated ephemeris using data from Spice.
/*
 *  Function to create a tabulated ephemeris using data from Spice.
 *  Retrieving a state from Spice can be very computationally intensive. Using this function,
 *  ephemeris data for a body is pre-computed using a limited number
 *  of calls to Spice, which is then used to create an interpolator (6th order Lagrange). For
 *  many numerical integration scenarios, this approach may be faster than using
 *  DirectSpiceEphemerisSettings, with negligible influence on accuracy.
 * \param body Name of body for which ephemeris data is to be retrieved.
 * \param initialTime Initial time from which interpolated data from Spice should be created.
 * \param endTime Final time from which interpolated data from Spice should be created.
 * \param timeStep Time step with which interpolated data from Spice should be created.
 * \param observerName Name of body relative to which the ephemeris is to be calculated.
 * \param referenceFrameName Orientatioan of the reference frame in which the epehemeris is to be
 *          calculated.
 * \return Tabulated ephemeris using data from Spice.
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< ephemerides::Ephemeris > createTabulatedEphemerisFromSpice(
        const std::string& body,
        const TimeType initialTime,
        const TimeType endTime,
        const TimeType timeStep,
        const std::string& observerName,
        const std::string& referenceFrameName,
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{
    using namespace interpolators;

    std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > timeHistoryOfState;

    // Calculate state from spice at given time intervals and store in timeHistoryOfState.
    TimeType currentTime = initialTime;
    while( currentTime < endTime )
    {
        timeHistoryOfState[ currentTime ] = spice_interface::getBodyCartesianStateAtEpoch(
                    body, observerName, referenceFrameName, "none", static_cast< double >( currentTime ) ).
                template cast< StateScalarType >( );
        currentTime += timeStep;
    }

    // Create interpolator.
    std::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > > interpolator =
            interpolators::createOneDimensionalInterpolator(
                timeHistoryOfState, interpolatorSettings );

    // Create ephemeris and return.
    return std::make_shared< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                interpolator, observerName, referenceFrameName );
}

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< ephemerides::Ephemeris > createTabulatedEphemerisFromTLE(
		const std::string& body,
		const TimeType initialTime,
		const TimeType endTime,
		const TimeType timeStep,
		const std::string& observerName,
		const std::string& referenceFrameName,
		std::shared_ptr< ephemerides::Tle > tle,
		std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
		std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{
	using namespace interpolators;

	std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > timeHistoryOfState;

	// Calculate state from spice at given time intervals and store in timeHistoryOfState.
	TimeType currentTime = initialTime;
	while( currentTime < endTime )
	{
        timeHistoryOfState[ currentTime ] = spice_interface::getCartesianStateFromTleAtEpoch( static_cast< double >( currentTime ), tle );
        currentTime += timeStep;
	}

	// Create interpolator.
	std::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > > interpolator =
			interpolators::createOneDimensionalInterpolator(
					timeHistoryOfState, interpolatorSettings );

	// Create ephemeris and return.
	return std::make_shared< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
			interpolator, observerName, referenceFrameName );
}

//! @get_docstring(keplerEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > keplerEphemerisSettings(
        const Eigen::Vector6d& initialStateInKeplerianElements,
        const double epochOfInitialState,
        const double centralBodyGravitationalParameter,
        const std::string& referenceFrameOrigin = "SSB",
        const std::string& referenceFrameOrientation = "ECLIPJ2000",
        const double rootFinderAbsoluteTolerance =
        200.0 * std::numeric_limits< double >::epsilon( ),
        const double rootFinderMaximumNumberOfIterations = 1000.0 )
{
    return std::make_shared< KeplerEphemerisSettings >(
                initialStateInKeplerianElements, epochOfInitialState, centralBodyGravitationalParameter,
                referenceFrameOrigin, referenceFrameOrientation, rootFinderAbsoluteTolerance,
                rootFinderMaximumNumberOfIterations );
}

//! @get_docstring(keplerEphemerisFromSpiceSettings)
inline std::shared_ptr< EphemerisSettings > keplerEphemerisFromSpiceSettings(
        const std::string body,
        const double epochOfInitialState,
        const double centralBodyGravitationalParameter,
        const std::string& referenceFrameOrigin = "SSB",
        const std::string& referenceFrameOrientation = "ECLIPJ2000",
        const double rootFinderAbsoluteTolerance =
        200.0 * std::numeric_limits< double >::epsilon( ),
        const double rootFinderMaximumNumberOfIterations = 1000.0 )
{
    Eigen::Vector6d initialCartesianState = spice_interface::getBodyCartesianStateAtEpoch(
                body, referenceFrameOrigin, referenceFrameOrientation, "None", epochOfInitialState );
    Eigen::Vector6d initialKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialCartesianState, centralBodyGravitationalParameter );
    return std::make_shared< KeplerEphemerisSettings >(
                initialKeplerianState, epochOfInitialState, centralBodyGravitationalParameter,
                referenceFrameOrigin, referenceFrameOrientation, rootFinderAbsoluteTolerance,
                rootFinderMaximumNumberOfIterations );
}


//! @get_docstring(approximateJplEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > approximateJplEphemerisSettings(
        const std::string bodyName )
{
    return std::make_shared< ApproximateJplEphemerisSettings >(
            bodyName, false );
}

//inline std::shared_ptr< EphemerisSettings > approximatePlanetPositionsSettings( )
//{
//    return std::make_shared< ApproximateJplEphemerisSettings >(
//                ephemerides::undefined_body, false );
//}

//! @get_docstring(directSpiceEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > directSpiceEphemerisSettings(
		const std::string frameOrigin = "SSB",
		const std::string frameOrientation = "ECLIPJ2000",
		const bool correctForStellarAberration = false,
		const bool correctForLightTimeAberration = false,
        const bool convergeLightTimeAberration = false )
{
	return std::make_shared< DirectSpiceEphemerisSettings >(
			frameOrigin, frameOrientation, correctForStellarAberration,
            correctForLightTimeAberration, convergeLightTimeAberration );
}

//! @get_docstring(directSpiceEphemerisSettings,1)
inline std::shared_ptr< EphemerisSettings > directSpiceEphemerisSettings(
        const std::string frameOrigin = "SSB",
        const std::string frameOrientation = "ECLIPJ2000",
        const std::string bodyNameOverride = "" )
{
    return std::make_shared< DirectSpiceEphemerisSettings >(
            frameOrigin, frameOrientation, bodyNameOverride);
}

//! @get_docstring(interpolatedSpiceEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > interpolatedSpiceEphemerisSettings(
		double initialTime,
		double finalTime,
		double timeStep,
		std::string frameOrigin = "SSB",
		std::string frameOrientation = "ECLIPJ2000",
		std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ),
        const std::string bodyNameOverride = "" )
{
	return std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialTime, finalTime, timeStep, frameOrigin, frameOrientation, interpolatorSettings, bodyNameOverride );
}

//! @get_docstring(tabulatedEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > tabulatedEphemerisSettings(
		const std::map< double, Eigen::Vector6d >& bodyStateHistory,
		std::string frameOrigin = "SSB",
		std::string frameOrientation = "ECLIPJ2000" )
{
	return std::make_shared< TabulatedEphemerisSettings >(
			bodyStateHistory, frameOrigin, frameOrientation	);
}

inline std::shared_ptr< EphemerisSettings > tabulatedEphemerisSettings(
        const std::shared_ptr< EphemerisSettings > ephemerisSettings,
        const double startTime,
        const double endTime,
        const double timeStep,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{
    interpolatorSettings->resetBoundaryHandling( interpolators::throw_exception_at_boundary );
    return std::make_shared< AutoGeneratedTabulatedEphemerisSettings >(
            ephemerisSettings, startTime, endTime, timeStep, interpolatorSettings );
}

//! @get_docstring(constantEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > constantEphemerisSettings(
		const Eigen::Vector6d& constantState,
		const std::string& frameOrigin = "SSB",
		const std::string& frameOrientation = "ECLIPJ2000" )
{
	return std::make_shared< ConstantEphemerisSettings >(
			constantState, frameOrigin, frameOrientation );
}

//! @get_docstring(customEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > customEphemerisSettings(
		const std::function< Eigen::Vector6d( const double ) > customStateFunction,
		const std::string& frameOrigin = "SSB",
		const std::string& frameOrientation = "ECLIPJ2000" )
{
	return std::make_shared< CustomEphemerisSettings >(
			customStateFunction, frameOrigin, frameOrientation );
}

inline std::shared_ptr< EphemerisSettings > directTleEphemerisSettings(
		std::shared_ptr< ephemerides::Tle > tle,
		const std::string frameOrigin = "Earth",
		const std::string frameOrientation = "J2000" )
{
	return std::make_shared< DirectTleEphemerisSettings >( tle, frameOrigin, frameOrientation );
}

inline std::shared_ptr< EphemerisSettings > directTleEphemerisSettingsFromTleLines(
    const std::string& tleLine1, const std::string& tleLine2,
    const std::string frameOrigin = "Earth",
    const std::string frameOrientation = "J2000" )
{
    return std::make_shared< DirectTleEphemerisSettings >( std::make_shared< ephemerides::Tle >( tleLine1, tleLine2 ), frameOrigin, frameOrientation );
}

inline std::shared_ptr< EphemerisSettings > interpolatedTleEphemerisSettings(
		const double initialTime, const double finalTime,
		const double timeStep, std::shared_ptr< ephemerides::Tle > tle,
		std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
		std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ),
		const std::string& frameOrigin = "Earth", const std::string& frameOrientation = "J2000" )
{
	return std::make_shared< InterpolatedTleEphemerisSettings >(
            initialTime, finalTime, timeStep, tle, interpolatorSettings,
			frameOrigin, frameOrientation );
}

//! @get_docstring(scaledEphemerisSettings)
inline std::shared_ptr< EphemerisSettings > scaledEphemerisSettings(
        const std::shared_ptr< EphemerisSettings > baseSettings,
        const double scaling,
        const bool isScalingAbsolute )
{
    return std::make_shared< ScaledEphemerisSettings >( baseSettings, scaling, isScalingAbsolute );
}

//! @get_docstring(scaledEphemerisSettings, 1)
inline std::shared_ptr< EphemerisSettings > scaledEphemerisSettings(
        const std::shared_ptr< EphemerisSettings > baseSettings,
        const Eigen::Vector6d scaling,
        const bool isScalingAbsolute )
{
    return std::make_shared< ScaledEphemerisSettings >( baseSettings, scaling, isScalingAbsolute );
}

//! @get_docstring(scaledEphemerisSettings, 2)
inline std::shared_ptr< EphemerisSettings > scaledEphemerisSettings(
        const std::shared_ptr< EphemerisSettings > baseSettings,
        const std::function< Eigen::Vector6d( const double ) > scaling,
        const bool isScalingAbsolute )
{
    return std::make_shared< ScaledEphemerisSettings >( baseSettings, scaling, isScalingAbsolute );
}


// Function to create a ephemeris model.
/*
 *  Function to create a ephemeris model based on model-specific settings for the ephemeris.
 *  \param ephemerisSettings Settings for the ephemeris model that is to be created, defined
 *  a pointer to an object of class (derived from) EphemerisSettings.
 *  \param bodyName Name of the body for which the ephemeris model is to be created.
 *  \return Ephemeris model created according to settings in ephemerisSettings.
 */
template< typename StateScalarType = double , typename TimeType = double >
std::shared_ptr< ephemerides::Ephemeris > createBodyEphemeris(
        const std::shared_ptr< EphemerisSettings > ephemerisSettings,
        const std::string& bodyName )
{
    using namespace ephemerides;
    // Declare return object.
    std::shared_ptr< ephemerides::Ephemeris > ephemeris;

    if( ephemerisSettings->getMakeMultiArcEphemeris( ) )
    {
        std::map< double, std::shared_ptr< Ephemeris > > singleArcEphemerides;
        ephemerisSettings->resetMakeMultiArcEphemeris( false );

        singleArcEphemerides[ -std::numeric_limits< double >::lowest( ) ] = createBodyEphemeris< StateScalarType, TimeType >(
                    ephemerisSettings, bodyName );

        ephemeris = std::make_shared< MultiArcEphemeris >(
                    singleArcEphemerides, ephemerisSettings->getFrameOrigin( ),
                    ephemerisSettings->getFrameOrientation( ) );
    }
    else
    {

        // Check which type of ephemeris model is to be created.
        switch( ephemerisSettings->getEphemerisType( ) )
        {
        case direct_spice_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< DirectSpiceEphemerisSettings > directEphemerisSettings =
                    std::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
            if( directEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected direct spice ephemeris settings for body " + bodyName );
            }
            else
            {
                std::string inputName = ( directEphemerisSettings->getBodyNameOverride( ) == "" ) ?
                            bodyName : directEphemerisSettings->getBodyNameOverride( );

                // Create corresponding ephemeris object.
                ephemeris = std::make_shared< SpiceEphemeris >(
                            inputName,
                            directEphemerisSettings->getFrameOrigin( ),
                            directEphemerisSettings->getCorrectForStellarAberration( ),
                            directEphemerisSettings->getCorrectForLightTimeAberration( ),
                            directEphemerisSettings->getConvergeLighTimeAberration( ),
                            directEphemerisSettings->getFrameOrientation( ) );
            }
            break;
        }
        case interpolated_spice:
        {
            // Check consistency of type and class.
            std::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedEphemerisSettings =
                    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >(
                        ephemerisSettings );
            if( interpolatedEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected tabulated spice ephemeris settings for body " + bodyName );
            }
            else
            {
                // Since only the barycenters of planetary systems are included in the standard DE
                // ephemerides, append 'Barycenter' to body name.
                std::string inputName = ( interpolatedEphemerisSettings->getBodyNameOverride( ) == "" ) ?
                            bodyName : interpolatedEphemerisSettings->getBodyNameOverride( );
                ephemeris = createTabulatedEphemerisFromSpice< StateScalarType, TimeType >(
                            inputName,
                            interpolatedEphemerisSettings->getInitialTime( ),
                            interpolatedEphemerisSettings->getFinalTime( ),
                            interpolatedEphemerisSettings->getTimeStep( ),
                            interpolatedEphemerisSettings->getFrameOrigin( ),
                            interpolatedEphemerisSettings->getFrameOrientation( ),
                            interpolatedEphemerisSettings->getInterpolatorSettings( ) );
            }
        }
        break;
        case tabulated_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                    std::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
            if( tabulatedEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected tabulated ephemeris settings for body " + bodyName );
            }
            else
            {
                // Create corresponding ephemeris object.

                // Cast input history to required type.
                if( tabulatedEphemerisSettings->getBodyStateHistory( ).size( ) != 0 )
                {
                    std::map< double, Eigen::Vector6d > originalStateHistory =
                            tabulatedEphemerisSettings->getBodyStateHistory( );
                    std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > templatedStateHistory;

                    for( std::map< double, Eigen::Vector6d >::const_iterator stateIterator =
                         originalStateHistory.begin( ); stateIterator != originalStateHistory.end( ); stateIterator++ )
                    {
                        templatedStateHistory[ stateIterator->first ] = stateIterator->second.cast< StateScalarType>( );
                    }
                    ephemeris =
                            std::make_shared< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                                std::make_shared< interpolators::LagrangeInterpolator<
                                TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >
                                ( templatedStateHistory, 6,
                                  interpolators::huntingAlgorithm,
                                  interpolators::lagrange_cubic_spline_boundary_interpolation,
                                  interpolators::throw_exception_at_boundary ),
                                tabulatedEphemerisSettings->getFrameOrigin( ),
                                tabulatedEphemerisSettings->getFrameOrientation( ) );
                }
                else
                {
                    ephemeris = std::make_shared< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                                std::shared_ptr< interpolators::OneDimensionalInterpolator<
                                TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >( ),
                                tabulatedEphemerisSettings->getFrameOrigin( ),
                                tabulatedEphemerisSettings->getFrameOrientation( ) );
                }
            }
            break;
        }
        case auto_generated_tabulated_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< AutoGeneratedTabulatedEphemerisSettings > tabulatedEphemerisSettings =
                    std::dynamic_pointer_cast< AutoGeneratedTabulatedEphemerisSettings >( ephemerisSettings );
            if( tabulatedEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected auto-generate tabulated ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = getTabulatedEphemeris(
                            createBodyEphemeris( tabulatedEphemerisSettings->getEphemerisSettings( ), bodyName ),
                            tabulatedEphemerisSettings->getStartTime( ), tabulatedEphemerisSettings->getEndTime( ),
                            tabulatedEphemerisSettings->getTimeStep( ), tabulatedEphemerisSettings->getInterpolatorSettings( ) );
            }
            break;
        }
        case constant_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< ConstantEphemerisSettings > constantEphemerisSettings =
                    std::dynamic_pointer_cast< ConstantEphemerisSettings >( ephemerisSettings );
            if( constantEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected constant ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = std::make_shared< ConstantEphemeris >(
                            [ = ]( ){ return constantEphemerisSettings->getConstantState( ); },
                            constantEphemerisSettings->getFrameOrigin( ),
                            constantEphemerisSettings->getFrameOrientation( ) );
            }
            break;
        }
        case custom_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< CustomEphemerisSettings > customEphemerisSettings =
                    std::dynamic_pointer_cast< CustomEphemerisSettings >( ephemerisSettings );
            if( customEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected custom ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = std::make_shared< CustomEphemeris< double, double > >(
                            customEphemerisSettings->getCustomStateFunction( ),
                            customEphemerisSettings->getFrameOrigin( ),
                            customEphemerisSettings->getFrameOrientation( ) );
            }
            break;
        }

        case kepler_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
                    std::dynamic_pointer_cast< KeplerEphemerisSettings >( ephemerisSettings );
            if( keplerEphemerisSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected Kepler ephemeris settings for " + bodyName );
            }
            else
            {
                // Create ephemeris
                ephemeris = std::make_shared< KeplerEphemeris >(
                            keplerEphemerisSettings->getInitialStateInKeplerianElements( ),
                            keplerEphemerisSettings->getEpochOfInitialState( ),
                            keplerEphemerisSettings->getCentralBodyGravitationalParameter( ),
                            keplerEphemerisSettings->getFrameOrigin( ),
                            keplerEphemerisSettings->getFrameOrientation( ),
                            keplerEphemerisSettings->getRootFinderAbsoluteTolerance( ),
                            keplerEphemerisSettings->getRootFinderMaximumNumberOfIterations( ) );
            }
            break;
        }
        case approximate_planet_positions:
        {
            // Check consistency of type and class.
            std::shared_ptr< ApproximateJplEphemerisSettings > approximateEphemerisSettings =
                    std::dynamic_pointer_cast< ApproximateJplEphemerisSettings >(
                        ephemerisSettings );
            if( approximateEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected approximate ephemeris settings for body " + bodyName );
            }
            else
            {


                // Create corresponding ephemeris object.
                if( approximateEphemerisSettings->getUseCircularCoplanarApproximation( ) )
                {
                    ephemeris = std::make_shared< ApproximateJplCircularCoplanarEphemeris >(
                                approximateEphemerisSettings->getBodyName( ) );
                }
                else
                {
                    ephemeris = std::make_shared< ApproximateJplEphemeris >(
                                approximateEphemerisSettings->getBodyName( ) );
                }
            }
            break;
        }
        case direct_tle_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< DirectTleEphemerisSettings > directTleEphemerisSettings =
                    std::dynamic_pointer_cast< DirectTleEphemerisSettings >( ephemerisSettings );
            if( directTleEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error, expected direct TLE ephemeris settings for body " + bodyName );
            }
            // Check if the Earth is present in the simulation
            else
            {
                //std::string inputName = bodyName;

                // Check period of the satellite for correct SDP setting


                // Create corresponding ephemeris object.
                ephemeris = std::make_shared< TleEphemeris >(
                        directTleEphemerisSettings->getFrameOrigin(),
                        directTleEphemerisSettings->getFrameOrientation(),
                        directTleEphemerisSettings->getTle()
                        );
            }
            break;
        }
        case interpolated_tle_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< InterpolatedTleEphemerisSettings > interpolatedTleEphemerisSettings =
                    std::dynamic_pointer_cast< InterpolatedTleEphemerisSettings >(
                            ephemerisSettings );
            if( interpolatedTleEphemerisSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error, expected interpolated TLE ephemeris settings for body " + bodyName );
            }
            else
            {

                ephemeris = createTabulatedEphemerisFromTLE< double, double >(
                        bodyName,
                        interpolatedTleEphemerisSettings->getInitialTime( ),
                        interpolatedTleEphemerisSettings->getFinalTime( ),
                        interpolatedTleEphemerisSettings->getTimeStep( ),
                        interpolatedTleEphemerisSettings->getFrameOrigin( ),
                        interpolatedTleEphemerisSettings->getFrameOrientation( ),
                        interpolatedTleEphemerisSettings->getTle( ),
                        interpolatedTleEphemerisSettings->getInterpolatorSettings( )
                        );
            }
            break;
        }
        case scaled_ephemeris:
        {
            // Check consistency of type and class.
            std::shared_ptr< ScaledEphemerisSettings > scaledEphemeriSettings =
                    std::dynamic_pointer_cast< ScaledEphemerisSettings >(
                        ephemerisSettings );
            if( scaledEphemeriSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error, expected scaled ephemeris settings for body " + bodyName );
            }
            else
            {
                std::shared_ptr< Ephemeris > baseEphemeris = createBodyEphemeris(
                            scaledEphemeriSettings->getBaseSettings( ), bodyName );
                ephemeris = std::make_shared< ScaledEphemeris >(
                            baseEphemeris, scaledEphemeriSettings->getScaling( ), scaledEphemeriSettings->getIsScalingAbsolute( ) );
            }
            break;
        }
        default:
        {
            throw std::runtime_error(
                        "Error, did not recognize ephemeris model settings type " +
                        std::to_string( ephemerisSettings->getEphemerisType( ) ) );
        }
        }
    }
    return ephemeris;

}

// Function that retrieves the time interval at which an ephemeris can be safely interrogated
/*
 * Function that retrieves the time interval at which an ephemeris can be safely interrogated. For most ephemeris types,
 * this function returns the full range of double values ( lowest( ) to max( ) ). For the tabulated ephemeris, the interval
 * on which the interpolator inside this object is valid is checked and returned
 * \param ephemerisModel Ephemeris model for which the interval is to be determined.
 * \return The time interval at which the ephemeris can be safely interrogated
 */
std::pair< double, double > getSafeInterpolationInterval( const std::shared_ptr< ephemerides::Ephemeris > ephemerisModel );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEEPHEMERIS_H
