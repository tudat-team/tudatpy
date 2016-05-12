/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsBase.h"

namespace tudat
{

namespace simulation_setup
{

//! List of ephemeris models available in simulations
/*!
 *  List of ephemeris models available in simulations. Ephemeris models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum EphemerisType
{
    approximate_planet_positions,
    direct_spice_ephemeris,
    tabulated_ephemeris,
    interpolated_spice
};

//! Class for providing settings for ephemeris model.
/*!
 *  Class for providing settings for automatic ephemeris model creation. This class is a
 *  functional (base) class for settings of ephemeris models that require no information in
 *  addition to their type (and frame origin and orientation). Ephemeris model classes defining
 *  requiring additional information must be created using an object derived from this class.
 */
class EphemerisSettings
{
public:

    //! Constructor, sets type of ephemeris model.
    /*!
     *  Constructor, sets type of ephemeris model and frame origin and orientation.
     *  Settings for ephemeris models requiring additional information should be defined in a
     *  derived class.
     *  \param ephemerisType Type of ephemeris model that is to be created.
     *  \param frameOrigin Origin of frame in which ephemeris data is defined.
     *  \param frameOrientation Orientation of frame in which ephemeris data is defined.
     */
    EphemerisSettings( const EphemerisType ephemerisType,
                       const std::string& frameOrigin = "SSB",
                       const  std::string& frameOrientation = "ECLIPJ2000" ):
        ephemerisType_( ephemerisType ),
        frameOrigin_( frameOrigin ),
        frameOrientation_( frameOrientation ){ }

    //! Destructor
    virtual ~EphemerisSettings( ){ }

    //! Function to return type of gravity field model that is to be created.
    /*!
     *  Function to return type of gravity field model that is to be created.
     *  \return Type of gravity field model that is to be created.
     */
    EphemerisType getEphemerisType( ){ return ephemerisType_; }

    //! Function to return the origin of frame in which ephemeris data is defined.
    /*!
     *  Function to return the origin of frame in which ephemeris data is defined.
     *  \return Origin of frame in which ephemeris data is defined.
     */
    std::string getFrameOrigin( ){ return frameOrigin_; }

    //! Function to return the orientation of frame in which ephemeris data is defined.
    /*!
     *  Function to return the orientation of frame in which ephemeris data is defined.
     *  \return Orientation of frame in which ephemeris data is defined.
     */
    std::string getFrameOrientation( ){ return frameOrientation_;}

protected:

    //! Type of ephemeris model that is to be created.
    EphemerisType ephemerisType_;

    //! Origin of frame in which ephemeris data is defined.
    std::string frameOrigin_;

    //! Orientation of frame in which ephemeris data is defined.
    std::string frameOrientation_;
};

//! EphemerisSettings derived class for defining settings of an ephemeris linked directly to Spice.
class DirectSpiceEphemerisSettings: public EphemerisSettings
{
public:

    //! Constructor.
    /*! Constructor, sets the properties from which the Spice ephemeris is to be retrieved.
     * \param frameOrigin Name of body relative to which the ephemeris is to be calculated.
     * \param correctForStellarAbberation Boolean whether to correct for stellar Abberation in
     *          retrieved values of (observed state).
     * \param correctForLightTimeAbberation Boolean whether to correct for light time in
     *          retrieved values of (observed state).
     * \param convergeLighTimeAbberation Boolean whether to use single iteration or max. 3
     *          iterations for calculating light time.
     * \param frameOrientation Orientatioan of the reference frame in which the epehemeris is to be
     *          calculated.
     * \param ephemerisType Type of ephemeris that is to be created, always set to
     *          direct_spice_ephemeris when using this class directly. The derived class
     *          InterpolatedSpiceEphemerisSettings sets this parameter to a different value. Not
     *          to be changed by used.
     */
    DirectSpiceEphemerisSettings( const std::string frameOrigin = "SSB",
                                  const std::string frameOrientation = "ECLIPJ2000",
                                  const bool correctForStellarAbberation = 0,
                                  const bool correctForLightTimeAbberation = 0,
                                  const bool convergeLighTimeAbberation = 0,
                                  const EphemerisType ephemerisType = direct_spice_ephemeris ):
        EphemerisSettings( ephemerisType, frameOrigin, frameOrientation ),
        correctForStellarAbberation_( correctForStellarAbberation ),
        correctForLightTimeAbberation_( correctForLightTimeAbberation ),
        convergeLighTimeAbberation_( convergeLighTimeAbberation ){ }


    //! Destructor
    virtual ~DirectSpiceEphemerisSettings( ){ }

    //! Returns whether to correct for stellar abberation in retrieved values of (observed state).
    /*!
     *  Returns whether to correct for stellar abberation in retrieved values of (observed state).
     *  \return Boolean defining whether to correct for stellar abberation in retrieved
     *  values of (observed state).
     */
    bool getCorrectForStellarAbberation( ){ return correctForStellarAbberation_; }

    //! Returns whether to correct for light time in retrieved values of (observed state).
    /*!
     *  Returns whether to correct for light time in retrieved values of (observed state).
     *  \return Boolean defining whether to correct for light time in retrieved values of
     *  (observed state).
     */
    bool getCorrectForLightTimeAbberation( ){ return correctForLightTimeAbberation_; }

    //! Returns whether to use single iteration or max. 3 iterations for calculating light time.
    /*!
     *  Returns whether to use single iteration or max. 3 iterations for calculating light time.
     *  \return Boolean defining whether to use single iteration or max. 3 iterations for
     *  calculating light time.
     */
    bool getConvergeLighTimeAbberation( ){ return convergeLighTimeAbberation_; }
protected:

    //! Boolean whether to correct for stellar abberation in retrieved values of (observed state).
    bool correctForStellarAbberation_;

    //! Boolean whether to correct for light time in retrieved values of (observed state).
    bool correctForLightTimeAbberation_;

    //! Boolean whether to use single iteration or max. 3 iterations for calculating light time.
    bool convergeLighTimeAbberation_;
};

//! EphemerisSettings derived class for defining settings of a ephemeris interpolated from Spice
//! data.
/*!
 *  EphemerisSettings derived class for defining settings of a ephemeris interpolated from Spice
 *  data. Retrieving a state from Spice can be very computationally intensive. Using the
 *  settings of this class, ephemeris data for a body is pre-computed using a limited number
 *  of calls to Spice, which is then used to create an interpolator (6th order Lagrange). For
 *  many numerical integration scenarios, this approach may be faster than using
 *  DirectSpiceEphemerisSettings, with negligible influence on accuracy.
 */
class InterpolatedSpiceEphemerisSettings: public DirectSpiceEphemerisSettings
{
public:

    //! Constructor.
    /*! Constructor, sets the properties from which the tabulated spice data is to be created
     *  from which an ephemeris is to be created.
     * \param initialTime Initial time from which interpolated data from Spice should be created.
     * \param finalTime Final time from which interpolated data from Spice should be created.
     * \param timeStep Time step with which interpolated data from Spice should be created.
     * \param frameOrigin Name of body relative to which the ephemeris is to be calculated.
     * \param frameOrientation Orientatioan of the reference frame in which the epehemeris is to be
     *          calculated.
     */
    InterpolatedSpiceEphemerisSettings( double initialTime,
                                        double finalTime,
                                        double timeStep,
                                        std::string frameOrigin = "SSB",
                                        std::string frameOrientation = "ECLIPJ2000" ):
        DirectSpiceEphemerisSettings( frameOrigin, frameOrientation, 0, 0, 0,
                                      interpolated_spice ),
        initialTime_( initialTime ), finalTime_( finalTime ), timeStep_( timeStep ){ }

    //! Function to returns initial time from which interpolated data from Spice should be created.
    /*!
     *  Function to returns initial time from which interpolated data from Spice should be created.
     *  \return Initial time from which interpolated data from Spice should be created.
     */
    double getInitialTime( ){ return initialTime_; }

    //! Function to returns final time from which interpolated data from Spice should be created.
    /*!
     *  Function to returns final time from which interpolated data from Spice should be created.
     *  \return Final time from which interpolated data from Spice should be created.
     */
    double getFinalTime( ){ return finalTime_; }

    //! Function to returns time step with which interpolated data from Spice should be created.
    /*!
     *  Function to returns time step with which interpolated data from Spice should be created.
     *  \return Time step with which interpolated data from Spice should be created.
     */
    double getTimeStep( ){ return timeStep_; }

private:

    //! Initial time from which interpolated data from Spice should be created.
    double initialTime_;

    //! Final time from which interpolated data from Spice should be created.
    double finalTime_;

    //! Time step with which interpolated data from Spice should be created.
    double timeStep_;

};

//! EphemerisSettings derived class for defining settings of an approximate ephemeris for major
//! planets.
/*!
 *  EphemerisSettings derived class for defining settings of an approximate ephemeris for major
 *  planets, as inplemented in ApproximatePlanetPositions class and derived class,
 *  described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf.
 */
class ApproximatePlanetPositionSettings: public EphemerisSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param bodyIdentifier Parameter identifying for which body an ephemeris is to be created.
     *  \param useCircularCoplanarApproximation Boolean defining whether a circular, coplanar
     *  orbit of the body is to be assumed, or whether a non-zero inclination and long-period
     *  changes in the orbit are to be included.
     */
    ApproximatePlanetPositionSettings(
            const ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData
            bodyIdentifier,
            const bool useCircularCoplanarApproximation ):
        EphemerisSettings( approximate_planet_positions ),
        bodyIdentifier_( bodyIdentifier ),
        useCircularCoplanarApproximation_( useCircularCoplanarApproximation ){ }

    //! Function to return parameter identifying for which body an ephemeris is to be created.
    /*!
     *  Function to return parameter identifying for which body an ephemeris is to be created.
     *  \return Parameter identifying for which body an ephemeris is to be created.
     */
    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData getBodyIdentifier( )
    {
        return bodyIdentifier_;
    }

    //! Function to return whether a circular, coplanar orbit of the body is to be assumed.
    /*!
     *  Function to return whether a circular, coplanar orbit of the body is to be assumed.
     *  \return Boolean defining whether a circular, coplanar orbit of the body is to be assumed.
     */
    bool getUseCircularCoplanarApproximation( )
    {
        return useCircularCoplanarApproximation_;
    }

private:

    //! Parameter identifying for which body an ephemeris is to be created.
    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData bodyIdentifier_;

    //!  Boolean defining whether a circular, coplanar orbit of the body is to be assumed,
    /*!
     *  Boolean defining whether a circular, coplanar orbit of the body is to be assumed
     *  (creating an ApproximatePlanetPositionsCircularCoplanar object), or whether
     *  a non-zero inclination and long-period changes in the orbit are to be included
     *  (creating an ApproximatePlanetPositions object).
     */
     bool useCircularCoplanarApproximation_;
};

//! EphemerisSettings derived class for defining settings of an ephemeris created from tabulated
//! data.
/*!
 *  EphemerisSettings derived class for defining settings of an ephemeris created from tabulated
 *  data. Currently the use of an 6th order Lagrange interpolator is hardcoded, which is created
 *  from the data that is provided. Note that at the edges of the interpolation interval, a
 *  Cubic spline interpolator is used to suppres the influence of Runge's phenomenon.
 */
class TabulatedEphemerisSettings: public EphemerisSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param bodyStateHistory Data map (time as key, Cartesian state as values) defining data
     *  from which an interpolated ephemeris is to be created.
     * \param frameOrigin Name of body relative to which the ephemeris is to be calculated.
     * \param frameOrientation Orientatioan of the reference frame in which the epehemeris is to be
     *          calculated.
     */
    TabulatedEphemerisSettings(
            const std::map< double, basic_mathematics::Vector6d >& bodyStateHistory,
            std::string frameOrigin = "SSB",
            std::string frameOrientation = "ECLIPJ2000" ):
        EphemerisSettings( tabulated_ephemeris, frameOrigin, frameOrientation ),
        bodyStateHistory_( bodyStateHistory ){ }

    //! Function returning data map defining discrete data from which an ephemeris is to be created.
    /*!
     *  Function returning data map defining discrete data from which an ephemeris is to be created.
     *  \return Data map defining discrete data from which an ephemeris is to be created.
     */
    std::map< double, basic_mathematics::Vector6d > getBodyStateHistory( )
    { return bodyStateHistory_; }

private:

    //! Data map defining discrete data from which an ephemeris is to be created.
    /*!
     *  Data map (time as key, Cartesian state as values) defining data from which an interpolated
     *  ephemeris is to be created.
     */
    std::map< double, basic_mathematics::Vector6d > bodyStateHistory_;
};

//! Function to create a tabulated ephemeris using data from Spice.
/*!
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
 */
boost::shared_ptr< ephemerides::Ephemeris > createTabulatedEphemerisFromSpice(
        const std::string& body,
        const double initialTime,
        const double endTime,
        const double timeStep,
        const std::string& observerName,
        const std::string& referenceFrameName );

//! Function to create a ephemeris model.
/*!
 *  Function to create a ephemeris model based on model-specific settings for the ephemeris.
 *  \param ephemerisSettings Settings for the ephemeris model that is to be created, defined
 *  a pointer to an object of class (derived from) EphemerisSettings.
 *  \param bodyName Name of the body for which the ephemeris model is to be created.
 *  \return Ephemeris model created according to settings in ephemerisSettings.
 */
boost::shared_ptr< ephemerides::Ephemeris > createBodyEphemeris(
        const boost::shared_ptr< EphemerisSettings > ephemerisSettings,
        const std::string& bodyName );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEEPHEMERIS_H
