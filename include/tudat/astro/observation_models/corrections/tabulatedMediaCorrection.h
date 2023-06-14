/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          T. Moyer (2000), Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation,
 *              DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 *          J.A. Estefan and O.J. Sovers (1994), A Comparative Survey of Current and Proposed Tropospheric Refraction-Delay
 *              Models forDSN Radio Metric Data Calibration, JPL/NASA
 *          N. Jakowski, M.M. Hoque and C. Mayer (2011), A new global TEC model for estimating transionospheric radio wave
 *              propagation errors, Journal of Geodesy, 85:965–974
 *          J.A. Klobuchar (1975), A First-Order, Worldwide, Ionospheric, Time-Delay Algorithm, AIR FORCE CAMBRIDGE
 *              RESEARCH LABORATORIES
 *          820-013 TRK-2-23, Media Calibration Interface, Revision C (2008), DSN/JPL
 *          O. Olsen (2007), HELIOSAT - An orbit determination software with applications to deep space missions,
 *              University of Oslo.
 *
 */

#ifndef TUDAT_TABULATEDMEDIACORRECTION_H
#define TUDAT_TABULATEDMEDIACORRECTION_H

#include <cmath>
#include <vector>

#include "tudat/math/interpolators.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{

namespace observation_models
{

// Base class defining a tabulated media correction. The computation of the media correction should be implemented
// in a derived class.
class TabulatedMediaReferenceCorrection
{
public:

    /*!
     * Constructor. Requires specifying the initial and final time during which the media correction is valid.
     * Initial/final time should be set to TUDAT_NAN if there is no initial/final limit to the validity of there
     * corrections.
     *
     * @param startTime Start time of the validity of the media correction.
     * @param endTime End time of the validity of the media correction.
     */
    TabulatedMediaReferenceCorrection( const double startTime = TUDAT_NAN,
                                       const double endTime = TUDAT_NAN ):
       startTime_( startTime ),
       endTime_( endTime )
    { }

    // Destructor
    virtual ~TabulatedMediaReferenceCorrection( ){ }

    /*!
     * Function to compute the atmospheric correction as a function of time. Should be implemented in the derived
     * class.
     *
     * @param time Time at which to compute the correction.
     * @return Atmospheric correction value
     */
    virtual double computeReferenceCorrection( const double time ) = 0;

    // Return start time of media-correction validity.
    double getStartTime( )
    {
        return startTime_;
    }

    // Return end time of media-correction validity.
    double getEndTime( )
    {
        return endTime_;
    }

protected:

    /*!
     * Checks whether the selected time is within the bounds for which the atmospheric correction is defined.
     * Throws error if the specified time isn't valid.
     *
     * @param time Time at which to compute the correction.
     * @return Boolean indicating whether time is valid or not.
     */
    bool isTimeValid( const double time );

    // Start time of the media-correction validity.
    const double startTime_;

    // End time of the media-correction validity.
    const double endTime_;

private:

};

// Derived class implementing a constant tabulated atmospheric correction, according to TRK-2-23 (2008),
// section 3.1.7.
class ConstantReferenceCorrection: public TabulatedMediaReferenceCorrection
{
public:

    /*!
     * Constructor.

     * @param startTime Start time of the validity of the media correction.
     * @param endTime End time of the validity of the media correction.
     * @param constantCorrection Constant value of the correction.
     */
    ConstantReferenceCorrection( const double startTime,
                                 const double endTime,
                                 const double constantCorrection ):
        TabulatedMediaReferenceCorrection( startTime, endTime ),
        constantCorrection_( constantCorrection )
    { }

    /*!
     * Function to compute the atmospheric correction as a function of time. Returns the value of the constant
     * correction, according to TRK-2-23 (2008), section 3.1.7.
     *
     * @param time Time at which to compute the correction.
     * @return Atmospheric correction value.
     */
    double computeReferenceCorrection( const double time ) override
    {
        isTimeValid( time );

        return constantCorrection_;
    }

    // Return the constant correction
    double getConstantCorrection( )
    {
        return constantCorrection_;
    }

private:

    // Value of the constant correction
    const double constantCorrection_;

};

// Derived class implementing a power series tabulated atmospheric correction, according to TRK-2-23 (2008),
// section 3.1.7.
class PowerSeriesReferenceCorrection: public TabulatedMediaReferenceCorrection
{
public:

    /*!
     * Constructor.
     *
     * @param startTime Start time of the validity of the media correction.
     * @param endTime End time of the validity of the media correction.
     * @param coefficients Coefficients of the power series. The i^th coefficient applies to the power of order i.
     */
    PowerSeriesReferenceCorrection( const double startTime,
                                    const double endTime,
                                    const std::vector< double > coefficients ):
        TabulatedMediaReferenceCorrection( startTime, endTime ),
        coefficients_( coefficients )
    { }

    /*!
     * Function to compute the atmospheric correction as a function of time. Correction is computed according to
     * a power series, according to TRK-2-23 (2008), section 3.1.7.
     *
     * @param time Time at which to compute the correction.
     * @return Atmospheric correction value.
     */
    double computeReferenceCorrection( const double time ) override;

    // Returns the power series coefficients
    std::vector< double > getCoefficients( )
    {
        return coefficients_;
    }

private:

    // Power series coefficients
    const std::vector< double > coefficients_;

};

// Derived class implementing a Fourier series tabulated atmospheric correction, according to TRK-2-23 (2008),
// section 3.1.7.
class FourierSeriesReferenceCorrection: public TabulatedMediaReferenceCorrection
{
public:

   /*!
    * Constructor.
    *
    * @param startTime Start time of the validity of the media correction.
    * @param endTime End time of the validity of the media correction.
    * @param coefficients Coefficients of the Fourier series. Order of the coefficients is described by TRK-2-23 (2008),
    *      section 3.1.7.
    */
    FourierSeriesReferenceCorrection( const double startTime,
                                      const double endTime,
                                      const std::vector< double > coefficients );

    /*!
    * Function to compute the atmospheric correction as a function of time. Correction is computed according to
    * a Fourier series, according to TRK-2-23 (2008), section 3.1.7.
    *
    * @param time Time at which to compute the correction.
    * @return Atmospheric correction value.
    */
    double computeReferenceCorrection( const double time ) override;

    // Returns the sine coefficients of the Fourier series.
    std::vector< double > getSineCoefficients( )
    {
        return sineCoefficients_;
    }

    // Returns the cosine coefficients of the Fourier series.
    std::vector< double > getCosineCoefficients( )
    {
        return cosineCoefficients_;
    }

private:

    // Fundamental period of the series
    double period_;

    // Sine coefficients of the Fourier series.
    std::vector< double > sineCoefficients_;

    // Cosine coefficients of the Fourier series.
    std::vector< double > cosineCoefficients_;

};

// Class for managing multiple objects computing media corrections. Effectively serves to select the appropriate
// for computing the correction at a given time.
class TabulatedMediaReferenceCorrectionManager
{
public:

    /*!
     * Constructor.
     */
    TabulatedMediaReferenceCorrectionManager( ):
        isLookupSchemeUpdated_( false )
    { }

    /*!
     * Constructor.
     *
     * @param correctionVector Vector of objects computing media corrections.
     */
    TabulatedMediaReferenceCorrectionManager(
            std::vector< std::shared_ptr< TabulatedMediaReferenceCorrection > > correctionVector ):
        correctionVector_( correctionVector ),
        isLookupSchemeUpdated_( false )
    {
        for ( unsigned int i = 0; i < correctionVector_.size( ); ++i )
        {
            startTimes_.push_back( correctionVector_.at( i )->getStartTime( ) );
            endTimes_.push_back( correctionVector_.at( i )->getEndTime( ) );
        }
    }

    /*!
     * Saves a correction calculation object.
     *
     * @param correctionCalculator Correction calculation object.
     */
    void pushReferenceCorrectionCalculator( std::shared_ptr< TabulatedMediaReferenceCorrection > correctionCalculator )
    {
        if ( correctionVector_.empty( ) ||
            ( !correctionVector_.empty( ) && correctionCalculator->getStartTime( ) > startTimes_.back( ) &&
                correctionCalculator->getEndTime( ) > endTimes_.back( ) ) )
        {
            startTimes_.push_back( correctionCalculator->getStartTime( ) );
            endTimes_.push_back( correctionCalculator->getEndTime( ) );
            correctionVector_.push_back( correctionCalculator );
        }
        else
        {
            throw std::runtime_error("Inconsistency in times when pushing tabulated media reference correction calculator. ");
        }

        isLookupSchemeUpdated_ = false;
    }

    /*!
    * Function to compute the atmospheric correction as a function of time.
    * For the specified time, the function selects the appropriate correction calculation object and the uses it
    * to compute the time.
    *
    * @param time Time at which to compute the correction.
    * @return Atmospheric correction value.
    */
    double computeMediaCorrection( double time );

private:

    // Vector with the start validity times of the corrections.
    std::vector< double > startTimes_;

    // Vector with the end validity times of the corrections.
    std::vector< double > endTimes_;

    // Vector containing the objects that compute the media corrections
    std::vector< std::shared_ptr< TabulatedMediaReferenceCorrection > > correctionVector_;

    // Boolean indicating whether the lookup scheme is updated (it gets out of date in case the start times vector
    // is modified)
    bool isLookupSchemeUpdated_;

    //! Lookup scheme to find the nearest correction object start time for a given time
    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;
};

// Enum listing the available tropospheric mapping models.
enum TroposphericMappingModel
{
    simplified_chao,
    niell
};

// Base class defining a tropospheric elevation mapping function (used to map a zenith tropospheric correction
// to the desired elevation)
class TroposhericElevationMapping
{
public:

    /*!
     * Constructor.
     */
    TroposhericElevationMapping( )
    { }

    // Destructor
    virtual ~TroposhericElevationMapping( ){ }

    /*!
     * Computes the factor that maps the dry troposheric correction. The computation of this factor should be
     * implemented in a derived class.
     *
     * @param transmitterState State of the transmitter at transmission.
     * @param receiverState State of the receiver at reception.
     * @param transmissionTime Transmission time.
     * @param receptionTime Reception time.
     * @return Dry mapping factor.
     */
    virtual double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) = 0;

    /*!
     * Computes the factor that maps the wet troposheric correction. The computation of this factor should be
     * implemented in a derived class.
     *
     * @param transmitterState State of the transmitter at transmission.
     * @param receiverState State of the receiver at reception.
     * @param transmissionTime Transmission time.
     * @param receptionTime Reception time.
     * @return Wet mapping factor.
     */
    virtual double computeWetTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) = 0;

protected:

private:

};

// Used to map a zenith range correction to a different elevation
// Derived class that implements the simplified Chao tropospheric mapping function, according to Moyer (2000), section eqs. 10-8 to 10-10.
// Mapping function is used to map a zenith range correction to a different elevation.
class SimplifiedChaoTroposphericMapping: public TroposhericElevationMapping
{
public:

    /*!
     * Constructor
     * @param elevationFunction Function that computes the elevation as seen from the ground station, given the vector to
     *      the target and the current time.
     * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
      *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
     */
    SimplifiedChaoTroposphericMapping(
            std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
            bool isUplinkCorrection ):
        TroposhericElevationMapping( ),
        elevationFunction_( elevationFunction ),
        isUplinkCorrection_( isUplinkCorrection )
    { }

    /*!
     * Computes the factor that maps the dry troposheric correction, using the simplified Chao mapping model.
     * Implementation according to Moyer (2000), eqs. 10-8 and 10-9.
     *
     * @param transmitterState State of the transmitter at transmission.
     * @param receiverState State of the receiver at reception.
     * @param transmissionTime Transmission time.
     * @param receptionTime Reception time.
     * @return Dry mapping factor.
     */
    double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override
    {
        computeCurrentElevation( transmitterState, receiverState, transmissionTime, receptionTime );
        return troposphericSimplifiedChaoMapping( currentElevation_, true );
    }

    /*!
     * Computes the factor that maps the wet troposheric correction, using the simplified Chao mapping model.
     * Implementation according to Moyer (2000), eqs. 10-8 and 10-10.
     *
     * @param transmitterState State of the transmitter at transmission.
     * @param receiverState State of the receiver at reception.
     * @param transmissionTime Transmission time.
     * @param receptionTime Reception time.
     * @return Wet mapping factor.
     */
    double computeWetTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override
    {
        computeCurrentElevation( transmitterState, receiverState, transmissionTime, receptionTime );
        return troposphericSimplifiedChaoMapping( currentElevation_, false );
    }

private:

    /*!
     * Computes the wet/dry mapping factor, using the simplified Chao mapping model.
     * Implementation according to Moyer (2000), eqs. 10-8, 10-9, and 10-10
     *
     * @param elevation Current elevation.
     * @param dryCorrection Bool indicating whether to compute the dry mapping (if true) or the wet correction (if false).
     * @return Dry or wet mapping factor. 
     */
    double troposphericSimplifiedChaoMapping( const double elevation,
                                              const bool dryCorrection );

    /*!
     * Retrieves the current elevation.
     *
     * @param transmitterState State of the transmitter at transmission.
     * @param receiverState State of the receiver at reception.
     * @param transmissionTime Transmission time.
     * @param receptionTime Reception time.
     * @return Elevation.
     */
    void computeCurrentElevation(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime );

    // Value of the current elevation
    double currentElevation_;

    // Function that computes the elevation as seen from the ground station, given the vector to the target and the current time.
    std::function< double ( Eigen::Vector3d, double ) > elevationFunction_;

    // Boolean indicating whether the correction is for uplink or donwlink (necessary when computing the elevation)
    bool isUplinkCorrection_;

};

// Derived class that implements the Niell tropospheric mapping function, according to Moyer (2000), section 10.2.1.3.2.
// Mapping function is used to map a zenith range correction to a different elevation.
class NiellTroposphericMapping: public TroposhericElevationMapping
{
public:

    /*!
     * Constructor
     * @param elevationFunction Function that computes the elevation as seen from the ground station, given the vector to
     *      the target and the current time.
     * @param groundStationGeodeticPositionFunction Geodetic position of the ground station as a function of time
     * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
      *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
     */
    NiellTroposphericMapping(
            std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
            std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction,
            bool isUplinkCorrection );


    /*!
     * Computes the factor that maps the wet troposheric correction, using the Niell mapping model. Implementation
     * according to Moyer (2000), section 10.2.1.3.2.
     *
     * @param transmitterState State of the transmitter at transmission.
     * @param receiverState State of the receiver at reception.
     * @param transmissionTime Transmission time.
     * @param receptionTime Reception time.
     * @return Wet mapping factor.
     */
    double computeWetTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override;

   /*!
    * Computes the factor that maps the dry troposheric correction, using the Niell mapping model. Implementation
    * according to Moyer (2000), section 10.2.1.3.2.
    *
    * @param transmitterState State of the transmitter at transmission.
    * @param receiverState State of the receiver at reception.
    * @param transmissionTime Transmission time.
    * @param receptionTime Reception time.
    * @return Dry mapping factor.
    */
    double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override;

private:

    /*!
     * Computes the value of the function M, according to Moyer (2000), section 10.2.1.3.2, and Estefan and Sovers (1994),
     * eq. 51. The equation of Estefan and Sovers (1994) is wrong, the correct version is described by Moyer (2000).
     *
     * @param a Mapping coefficient a.
     * @param b Mapping coefficient b.
     * @param c Mapping coefficient c.
     * @param elevation Elevation of the spacecraft as seen from the ground station.
     * @return Value of M function.
     */
    double computeMFunction ( const double a, const double b, const double c, const double elevation );

    /*!
     * Computes a given dry mapping function coefficient as a function of time and latitude of the ground station,
     * according to Moyer (2000), section 10.2.1.3.2.
     *
     * @param averageInterpolator Interpolator for the average value of the mapping coefficient, giving the value of the
     *      average as a function of the latitude of the ground station.
     * @param amplitudeInterpolator Interpolator for the amplitude of the mapping coefficient, giving the value of the
     *      amplitude as a function of the latitude of the ground station.
     * @param time Time at which
     * @param geodeticLatitude Geodetic latitude of the ground station.
     * @return Dry mapping coefficient.
     */
    double computeDryCoefficient(
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > averageInterpolator,
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > amplitudeInterpolator,
            const double time,
            const double geodeticLatitude );

    // Function that computes the elevation as seen from the ground station, given the vector to the target and the current time.
    std::function< double ( Eigen::Vector3d, double ) > elevationFunction_;

    // Geodetic position of the ground station as a function of time
    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction_;

    // Boolean indicating whether the correction is for uplink or donwlink (necessary when computing the elevation)
    bool isUplinkCorrection_;

    // Estefan and Sovers (1994), table 4a/4b
    const std::vector< double > referenceGeodeticLatitudes_ = { 15.0 * mathematical_constants::PI / 180.0,
                                                                30.0 * mathematical_constants::PI / 180.0,
                                                                45.0 * mathematical_constants::PI / 180.0,
                                                                60.0 * mathematical_constants::PI / 180.0,
                                                                75.0 * mathematical_constants::PI / 180.0 };

    // Estefan and Sovers (1994), table 4a
    const std::vector< double > aDryAverage_ = { 1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3 };
    const std::vector< double > bDryAverage_ = { 2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3 };
    const std::vector< double > cDryAverage_ = { 62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3 };

    // Estefan and Sovers (1994), table 4a
    const std::vector< double > aDryAmplitude_ = { 0.0e-5, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5 };
    const std::vector< double > bDryAmplitude_ = { 0.0e-5, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5 };
    const std::vector< double > cDryAmplitude_ = { 0.0e-5, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5 };

    // Estefan and Sovers (1994), eqs. 55a, 55b, 55c
    const double aHt_ = 2.53e-5;
    const double bHt_ = 5.49e-3;
    const double cHt_ = 1.14e-3;

    // Estefan and Sovers (1994), table 4b
    const std::vector< double > aWet_ = { 5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4 };
    const std::vector< double > bWet_ = { 1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3 };
    const std::vector< double > cWet_ = { 4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2 };

    // Interpolators to compute the mapping coefficients for the desired latitude values
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > aDryAverageInterpolator_;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > bDryAverageInterpolator_;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > cDryAverageInterpolator_;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > aDryAmplitudeInterpolator_;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > bDryAmplitudeInterpolator_;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > cDryAmplitudeInterpolator_;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > aWetInterpolator_;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > bWetInterpolator_;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > cWetInterpolator_;
};

// Base class implementing a mapped tropospheric light-time correction, according to Moyer (2000), section 10.2.1.
// Class implements the computation of the light time. This computation requires functions that return the wet/dry zenith
// corrections, which should be implemented in derived classes.
class MappedTroposphericCorrection: public LightTimeCorrection
{
public:

    /*!
     * Constructor
     * @param lightTimeCorrectionType Type of light-time correction represented by instance of class.
     * @param elevationMapping Mapping function of range corrections from zenith to other elevations
     * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
      *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
     * @param dryZenithRangeCorrectionFunction Function computing the dry atmosphere zenith correction for a given time
     * @param wetZenithRangeCorrectionFunction Function computing the wet atmosphere zenith correction for a given time
     */
    MappedTroposphericCorrection(
            const LightTimeCorrectionType lightTimeCorrectionType,
            std::shared_ptr< TroposhericElevationMapping > elevationMapping,
            bool isUplinkCorrection,
            std::function< double ( double time ) > dryZenithRangeCorrectionFunction = [] ( double ) { return 0.0; },
            std::function< double ( double time ) > wetZenithRangeCorrectionFunction = [] ( double ) { return 0.0; } ):
        LightTimeCorrection( lightTimeCorrectionType ),
        dryZenithRangeCorrectionFunction_( dryZenithRangeCorrectionFunction ),
        wetZenithRangeCorrectionFunction_( wetZenithRangeCorrectionFunction ),
        elevationMapping_( elevationMapping ),
        isUplinkCorrection_( isUplinkCorrection )
    { }

    /*!
      * Function to computes a mapped tropospheric light-time correction. Correction is computed according to Moyer
      * (2000), section 10.2.1, using dry/wet zenith correction functions defined in the derived classes.
      *
      * @param linkEndsStates List of states at each link end during observation.
      * @param linkEndsTimes List of times at each link end during observation.
      * @param currentMultiLegTransmitterIndex Index in the linkEndsStates and linkEndsTimes of the transmitter in the current link.
      * @param ancillarySettings Observation ancillary simulation settings.
      * @return Light-time correction
      */
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. observation time. Partial is
     * currently not implemented, function returns 0.
     *
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param fixedLinkEnd Reference link end for observation
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
     * \return Partial of light-time correction w.r.t. observation time
     */
    double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType fixedLinkEnd,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        return 0.0;
    }

    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. link end position. Partial is
     * currently not implemented, function returns 0.
     *
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
     * \return Partial of ight-time correction w.r.t. link end position
     */
    Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        return Eigen::Vector3d::Zero( );
    }

    // Returns the function that computes the dry zenith correction as a function of time.
    std::function< double ( double time ) > getDryZenithRangeCorrectionFunction( )
    {
        return dryZenithRangeCorrectionFunction_;
    }

    // Returns the function that computes the wet zenith correction as a function of time.
    std::function< double ( double time ) > getWetZenithRangeCorrectionFunction( )
    {
        return wetZenithRangeCorrectionFunction_;
    }

protected:

    // Dry atmosphere zenith range correction (in meters)
    std::function< double ( double time ) > dryZenithRangeCorrectionFunction_;

    // Wet atmosphere zenith range correction (in meters)
    std::function< double ( double time ) > wetZenithRangeCorrectionFunction_;

    // Object that maps the zenith tropospheric correction to the desired elevation.
    std::shared_ptr< TroposhericElevationMapping > elevationMapping_;

    // Boolean indicating whether the correction is for uplink or donwlink (necessary when computing the elevation)
    bool isUplinkCorrection_;
};

// Class to compute the tabulated tropospheric corrections using DSN data, according to Moyer (2000), section 10.2.1.
// Derived class from MappedTroposphericCorrection.
class TabulatedTroposphericCorrection: public MappedTroposphericCorrection
{
public:

    /*!
     * Constructor
     * @param seasonalModelDryZenithCorrectionCalculator, seasonalModelWetZenithCorrectionCalculator Correction calculators
     *      based on seasonal model. These corrections should either correspond to Figure 3a or 3b of Estefan and Sovers (1994).
     * @param dryZenithCorrectionAdjustmentCorrectionCalculator, wetZenithCorrectionAdjustmentCorrectionCalculator Corrections
     *      to the seasonal model based on real time data. Should be read from DSN TRK-2-23 files.
     * @param elevationMapping Mapping function of range corrections from zenith to other elevations
     * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
      *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
     */
    TabulatedTroposphericCorrection(
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > seasonalModelDryZenithCorrectionCalculator,
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > seasonalModelWetZenithCorrectionCalculator,
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > dryZenithCorrectionAdjustmentCalculator,
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > wetZenithCorrectionAdjustmentCalculator,
            std::shared_ptr< TroposhericElevationMapping > elevationMapping,
            bool isUplinkCorrection ):
        MappedTroposphericCorrection(
                tabulated_tropospheric,
                elevationMapping,
                isUplinkCorrection ),
        seasonalModelDryZenithCorrectionCalculator_( seasonalModelDryZenithCorrectionCalculator ),
        seasonalModelWetZenithCorrectionCalculator_( seasonalModelWetZenithCorrectionCalculator ),
        dryZenithCorrectionAdjustmentCalculator_( dryZenithCorrectionAdjustmentCalculator ),
        wetZenithCorrectionAdjustmentCalculator_( wetZenithCorrectionAdjustmentCalculator )
    {
        // Override default values for dry and wet zenith corrections
        // TRK-2-23 (2008), section 3.2.2
        dryZenithRangeCorrectionFunction_ = [=] ( double time ) {
            return seasonalModelDryZenithCorrectionCalculator_->computeMediaCorrection( time ) +
                dryZenithCorrectionAdjustmentCalculator_->computeMediaCorrection( time ); };
        wetZenithRangeCorrectionFunction_ = [=] ( double time ) {
            return seasonalModelWetZenithCorrectionCalculator_->computeMediaCorrection( time ) +
                wetZenithCorrectionAdjustmentCalculator_->computeMediaCorrection( time ); };
    }

private:

    // Dry atmosphere zenith range correction (in meters): based on seasonal model
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > seasonalModelDryZenithCorrectionCalculator_;

    // Wet atmosphere zenith range correction (in meters): based on seasonal model
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > seasonalModelWetZenithCorrectionCalculator_;

    // Dry atmosphere zenith range correction (in meters): tabulated correction to seasonal model
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > dryZenithCorrectionAdjustmentCalculator_;

    // Wet atmosphere zenith range correction (in meters): tabulated correction to seasonal model
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > wetZenithCorrectionAdjustmentCalculator_;
};

//! Enum defining different types of water vapor partial pressure models.
enum WaterVaporPartialPressureModel
{
    tabulated,
    bean_and_dutton
};

/*! Calculate the partial vapor pressure.
 *
 * Calculate the partial vapor pressure according to the Bean and Dutton (1966) model, as described by Estefan and Sovers
 * (1994), Eq. 16.
 *
 * @param relativeHumidity Relative humidity, defined in [0,1]
 * @param temperature Temperature in Kelvin
 * @return Partial vapor pressure in Pa
 */
double calculateBeanAndDuttonWaterVaporPartialPressure( double relativeHumidity,
                                                        double temperature );

/*!
 * Returns a function that computes the water vapor partial pressure as a function of time, according to the Bean and
 * Dutton (1966) model.
 *
 * @param relativeHumidity Relative humidity as a function of time.
 * @param temperature Temperature as a function of time.
 * @return Water vapor partial pressure as a function of time.
 */
std::function< double ( const double time ) > getBeanAndDuttonWaterVaporPartialPressureFunction(
        std::function< double ( const double time ) > relativeHumidity,
        std::function< double ( const double time ) > temperature );

// Class to compute the Saastamaoinen tropospheric corrections, according to Estefan and Sovers (1994). Derived class
// from MappedTroposphericCorrection.
class SaastamoinenTroposphericCorrection: public MappedTroposphericCorrection
{
public:

    /*!
     * Constructor
     * @param groundStationGeodeticPositionFunction  Geodetic position of the ground station as a function of time
     * @param pressureFunction Pressure at the ground station as a function of time
     * @param temperatureFunction Temperature at the ground station as a function of time
     * @param waterVaporPartialPressureFunction Water vapor partial pressure at the ground station as a function of time
     * @param elevationMapping Mapping function of range corrections from zenith to other elevations
     * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
     *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
     */
    SaastamoinenTroposphericCorrection(
            std::function< Eigen::Vector3d ( const double time ) > groundStationGeodeticPositionFunction,
            std::function< double ( const double time ) > pressureFunction,
            std::function< double ( const double time ) > temperatureFunction,
            std::function< double ( const double time ) > waterVaporPartialPressureFunction,
            std::shared_ptr< TroposhericElevationMapping > elevationMapping,
            bool isUplinkCorrection ):
        MappedTroposphericCorrection(
                saastamoinen_tropospheric,
                elevationMapping,
                isUplinkCorrection ),
        groundStationGeodeticPositionFunction_( groundStationGeodeticPositionFunction ),
        pressureFunction_( pressureFunction ),
        temperatureFunction_( temperatureFunction ),
        waterVaporPartialPressureFunction_( waterVaporPartialPressureFunction )
    {
        // Override default values for dry and wet zenith corrections
        dryZenithRangeCorrectionFunction_ = std::bind(
                &SaastamoinenTroposphericCorrection::computeDryZenithRangeCorrection, this, std::placeholders::_1 );
        wetZenithRangeCorrectionFunction_ = std::bind(
                &SaastamoinenTroposphericCorrection::computeWetZenithRangeCorrection, this, std::placeholders::_1 );
    }

private:

    // Computes the dry atmosphere zenith range correction (in meters)
    double computeDryZenithRangeCorrection( const double stationTime );

    // Computes the wet atmosphere zenith range correction (in meters)
    double computeWetZenithRangeCorrection( const double stationTime );

    // Geodetic position of the ground station as a function of time
    std::function< Eigen::Vector3d ( const double ) > groundStationGeodeticPositionFunction_;

    // Pressure at the ground station as a function of time
    std::function< double ( const double ) > pressureFunction_;

    // Temperature at the ground station as a function of time
    std::function< double ( const double ) > temperatureFunction_;

    // Water vapor partial pressure at the ground station as a function of time
    std::function< double ( const double ) > waterVaporPartialPressureFunction_;
};

// Tabulated ionospheric corrections using DSN data, according to Moyer (2000), section 10.2.2
class TabulatedIonosphericCorrection: public LightTimeCorrection
{
public:

     /*!
      * Constructor
      * @param referenceCorrectionCalculator Range correction calculator. Corrections based on real time data, which should
      *      read from DSN TRK-2-23 files.
      * @param transmittedFrequencyFunction Function calculating the frequency at the current link given a vector with
      *     the frequency bands in each link of the model and the transmission time.
      * @param baseObservableType Observable type associated with the correction.
      * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
      *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
      * @param referenceFrequency Frequency for which the reference corrections are given.
      */
    TabulatedIonosphericCorrection(
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > referenceCorrectionCalculator,
            std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction,
            ObservableType baseObservableType,
            bool isUplinkCorrection,
            double referenceFrequency = 2295e6 );

    /*!
      * Function to compute the ionospheric light-time correction, using tabulated DSN data, according to Moyer (2000),
      * section 10.2.2.
      *
      * @param linkEndsStates List of states at each link end during observation.
      * @param linkEndsTimes List of times at each link end during observation.
      * @param currentMultiLegTransmitterIndex Index in the linkEndsStates and linkEndsTimes of the transmitter in the current link.
      * @param ancillarySettings Observation ancillary simulation settings.
      * @return Light-time correction
      */
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. observation time. Partial is
     * currently not implemented, function returns 0.
     *
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param fixedLinkEnd Reference link end for observation
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
     * \return Partial of light-time correction w.r.t. observation time
     */
    double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType fixedLinkEnd,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        // TODO: Add computation of partial
        return 0.0;
    }

    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. link end position. Partial is
     * currently not implemented, function returns 0.
     *
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
     * \return Partial of ight-time correction w.r.t. link end position
     */
    Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        // TODO: Add computation of partial
        return Eigen::Vector3d::Zero( );
    }

private:

    // Range correction calculator. Correction determined for referenceFrequency_
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > referenceCorrectionCalculator_;

    // Frequency at the link as a function of the frequency bands per link, and of the current time
    std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction_;

    // Reference frequency for which the reference corrections where calculated
    double referenceFrequency_;

    // Sign of the correction (+1 or -1)
    int sign_;

    // Boolean indicating whether the correction is for uplink or downlink
    bool isUplinkCorrection_;

};

// Base class for calculating the vertical total electron content (VTEC) of the ionosphere.
class VtecCalculator
{
public:

    /*!
     * Constructor.
     *
     * @param referenceIonosphereHeight Reference height of the ionosphere layer used when computing VTEC.
     */
    VtecCalculator( const double referenceIonosphereHeight ):
        referenceIonosphereHeight_( referenceIonosphereHeight )
    { }

    // Destructor
    virtual ~VtecCalculator( ){ }

    /*!
     * Function to calculate the VTEC in [m^-2]. (m^-2 = 1e16 TECU). Pure virtual function, should be implemented in
     * derived class.
     *
     * @param time Time at which to compute the VTEC
     * @param subIonosphericPointGeodeticPosition Geodetic position of the sub-ionospheric point where to compute the VTEC.
     * @return VTEC
     */
    virtual double calculateVtec( const double time,
                                  const Eigen::Vector3d subIonosphericPointGeodeticPosition ) = 0;

    // Returns the reference ionoshere height
    double getReferenceIonosphereHeight( )
    {
        return referenceIonosphereHeight_;
    }

private:

    // Reference height of the ionosphere layer used when computing VTEC.
    const double referenceIonosphereHeight_;

};

// Derived class implementing the computation of the VTEC using the model from Jakowski et al. (2011)
// The geomagnetic latitude is computed according to Klobuchar (1975)
// The local time is computed according to Moyer (2000)
class JakowskiVtecCalculator: public VtecCalculator
{
public:

    /*!
     * Constructor
     *
     * Default values copied from GODOT
     * - Geomagnetic pole latitude and longitude:
     *      "Geomagnetic North pole latitude and longitude, values for 2025 from http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
     *      Variation of Geomagnetic poles is historically slow with very minor effect on model results
     *      e.g. (1.5,2.7) deg difference between 1960 and 2010, resulting in <3 cm effect on range"
     * - Reference ionosphere height:
     *      "400km confirmed as good choice by Jakowski in email to G.Bellei, 11.11.2014"
     *
     * @param sunDeclinationFunction Declination of the Sun as seen from the ground station as a function of time
     * @param observedSolarRadioFlux107Function Observed F10.7 flux as a function of time
     * @param useUtcTimeForLocalTime Boolean indicating whether to use UTC (if true) or TDB (if false) time when computing
     *      the local time. According to Jakowski, UT1 should be used, but UTC is the closest time scale that doesn't
     *      require loading file data. Anyway, TDB is most likely accurate enough, hence the default is here selected to
     *      use TDB.
     * @param geomagneticPoleLatitude Latitude of the geomagnetic pole
     * @param geomagneticPoleLongitude Longitude of the geomagnetic pole
     * @param referenceIonosphereHeight Ionosphere height
     */
    JakowskiVtecCalculator(
            std::function< double ( const double time ) > sunDeclinationFunction,
            std::function< double ( const double time ) > observedSolarRadioFlux107Function,
            bool useUtcTimeForLocalTime = false,
            double geomagneticPoleLatitude = unit_conversions::convertDegreesToRadians( 80.9 ),
            double geomagneticPoleLongitude = unit_conversions::convertDegreesToRadians( -72.6 ),
            double referenceIonosphereHeight = 400.0e3 ):
        VtecCalculator( referenceIonosphereHeight ),
        sunDeclinationFunction_( sunDeclinationFunction ),
        observedSolarRadioFlux107Function_( observedSolarRadioFlux107Function ),
        geomagneticPoleLatitude_( geomagneticPoleLatitude ),
        geomagneticPoleLongitude_( geomagneticPoleLongitude )
    {
        if ( useUtcTimeForLocalTime )
        {
            timeScaleConverter_ = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( );
        }
        else
        {
            timeScaleConverter_ = nullptr;
        }
    }

    /*!
     * Function to calculate the VTEC in [m^-2]. (m^-2 = 1e16 TECU).
     *
     * @param time Time at which to compute the VTEC
     * @param subIonosphericPointGeodeticPosition Geodetic position of the sub-ionospheric point where to compute the VTEC.
     * @return VTEC
     */
    double calculateVtec( const double time,
                          const Eigen::Vector3d subIonosphericPointGeodeticPosition ) override;

private:

    /*!
     * Gets the time used to compute the local time. If a time scale converter was provided, this time is obtained by
     * converting from TDB from UTC. If no time scale converter is provided, the function simply returns the TDB time.
     *
     * @param time TDB time at ground station.
     * @return Time (TDB or UTC).
     */
    double getUtcTime( double time )
    {
        if ( timeScaleConverter_ == nullptr )
        {
            return time;
        }
        else
        {
            return timeScaleConverter_->getCurrentTime(
                    basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time );
        }
    }

    // Coefficients of Jakowski model, Jakowski et al. (2011), tab. 1
    const std::vector< double > jakowskiCoefficients_ = { 0.89656, 0.16984, -0.02166, 0.05928, 0.00738, 0.13912,
                                                          -0.17593, -0.34545, 1.1167, 1.1573, -4.3356, 0.17775 };

    // Declination of the Sun as seen from the ground station as a function of time
    const std::function< double ( const double time ) > sunDeclinationFunction_;

    //! Observed (unadjusted) value of F10.7. Expressed in units of 10-22 W/m2/Hz.
    const std::function< double ( const double time ) > observedSolarRadioFlux107Function_;

    // Latitude of the geomagnetic pole
    const double geomagneticPoleLatitude_;

    // Longitude of the geomagnetic pole
    const double geomagneticPoleLongitude_;

    // Time scale converter. Used to convert the time from TDB to UTC.
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter_;

};

// Computes the ionospheric delay by mapping the vertical TEC to slant TEC using a very simple mapping function, following
// Moyer (2000), section 10.3.1.
class MappedVtecIonosphericCorrection: public LightTimeCorrection
{
public:

    /*!
     * Constructor.
     * @param vtecCalculator Class to calculate the vertical total electron content (VTEC)
     * @param transmittedFrequencyFunction Function calculating the frequency at the current link given a vector with
     *     the frequency bands in each link of the model and the transmission time.
     * @param elevationFunction Function that computes the elevation as seen from the ground station, given the vector to
     *      the target and the current time.
     * @param azimuthFunction Function that computes the azimuth as seen from the ground station, given the vector to
     *      the target and the current time.
     * @param groundStationGeodeticPositionFunction Geodetic position of the ground station as a function of time
     * @param baseObservableType Observable type associated with the correction.
     * @param isUplinkCorrection Boolean indicating whether correction is for uplink (i.e. transmitting station on planet,
      *      reception on spacecraft) or downlink (i.e. transmission from spacecraft, reception at ground station)
     * @param bodyWithAtmosphereMeanEquatorialRadius Mean equatorial radius of the body with the ionosphere.
     * @param firstOrderDelayCoefficient Value of the 1st order delay coefficient. Default value from Jakowski (2011);
     *      also see IERS conventions 2010, eq. 9.24.
     */
    MappedVtecIonosphericCorrection(
            std::shared_ptr< VtecCalculator > vtecCalculator,
            std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction,
            std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
            std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > azimuthFunction,
            std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction,
            ObservableType baseObservableType,
            bool isUplinkCorrection,
            double bodyWithAtmosphereMeanEquatorialRadius,
            double firstOrderDelayCoefficient = 40.3 );

    /*!
      * Function to compute the ionospheric light-time correction, by mapping the vertical TEC to slant TEC using a very
      * simple mapping function, following Moyer (2000), section 10.3.1.
      *
      * Note 1: Function uses a very simple mapping function, following Moyer (2000). At some point, it might be worth
      *     using other mapping functions, in which case the part of the function where the mapping is executed should
      *     be moved to a new class.
      * Note 2: Currently the STEC is calculated from the VTEC using the algorithm for computing the sub-ionospheric point
      *     described by Moyer (2000). This algorithm fails if the line of sight passes near the poles. For an alternative
      *     method see Olsen (2007), section 5.7.
      *
      * @param linkEndsStates List of states at each link end during observation.
      * @param linkEndsTimes List of times at each link end during observation.
      * @param currentMultiLegTransmitterIndex Index in the linkEndsStates and linkEndsTimes of the transmitter in the current link.
      * @param ancillarySettings Observation ancillary simulation settings.
      * @return Light-time correction
      */
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. observation time. Partial is
     * currently not implemented, function returns 0.
     *
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param fixedLinkEnd Reference link end for observation
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
     * \return Partial of light-time correction w.r.t. observation time
     */
    double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType fixedLinkEnd,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        // TODO: Add computation of partial
        return 0.0;
    }

    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. link end position. Partial is
     * currently not implemented, function returns 0.
     *
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
     * \return Partial of ight-time correction w.r.t. link end position
     */
    Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        // TODO: Add computation of partial
        return Eigen::Vector3d::Zero( );
    }

private:

    // Class to calculate the vertical total electron content (VTEC)
    std::shared_ptr< VtecCalculator > vtecCalculator_;

    // Frequency at the link as a function of the frequency bands per link, and of the current time
    std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction_;

    // Function that computes the elevation as seen from the ground station, given the vector to the target and the current time.
    std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction_;

    // Function that computes the azimuth as seen from the ground station, given the vector to the target and the current time.
    std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > azimuthFunction_;

    // Geodetic position of the ground station as a function of time
    std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction_;

    // Mean equatorial radius of the body with the ionosphere
    const double bodyWithAtmosphereMeanEquatorialRadius_;

    // Value of the 1st order delay coefficient
    const double firstOrderDelayCoefficient_;

    // Sign of the correction (+1 or -1)
    int sign_;

    // Boolean indicating whether the correction is for uplink or downlink
    bool isUplinkCorrection_;

};

} // namespace observation_models

} // namespace tudat

#endif //TUDATBUNDLE_TABULATEDMEDIACORRECTION_H
