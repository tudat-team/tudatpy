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

class TabulatedMediaReferenceCorrection
{
public:

    TabulatedMediaReferenceCorrection( const double startTime = TUDAT_NAN,
                                       const double endTime = TUDAT_NAN ):
       startTime_( startTime ),
       endTime_( endTime )
    { }

    virtual ~TabulatedMediaReferenceCorrection( ){ }

    virtual double computeReferenceCorrection( const double time ) = 0;

    double getStartTime( )
    {
        return startTime_;
    }

    double getEndTime( )
    {
        return endTime_;
    }

protected:

    bool isTimeValid( const double time );

    const double startTime_;
    const double endTime_;

private:

};

// TRK-2-23 (2008), section 3.1.7
class ConstantReferenceCorrection: public TabulatedMediaReferenceCorrection
{
public:

    ConstantReferenceCorrection( const double startTime,
                                 const double endTime,
                                 const double constantCorrection ):
        TabulatedMediaReferenceCorrection( startTime, endTime ),
        constantCorrection_( constantCorrection )
    { }

    double computeReferenceCorrection( const double time ) override
    {
        isTimeValid( time );
        
        return constantCorrection_;
    }

    double getConstantCorrection( )
    {
        return constantCorrection_;
    }

private:

    const double constantCorrection_;

};

// TRK-2-23 (2008), section 3.1.7
class PowerSeriesReferenceCorrection: public TabulatedMediaReferenceCorrection
{
public:

    PowerSeriesReferenceCorrection( const double startTime,
                                    const double endTime,
                                    const std::vector< double > coefficients ):
        TabulatedMediaReferenceCorrection( startTime, endTime ),
        coefficients_( coefficients )
    { }

    double computeReferenceCorrection( const double time ) override;

    std::vector< double > getCoefficients( )
    {
        return coefficients_;
    }

private:

    const std::vector< double > coefficients_;

};

// TRK-2-23 (2008), section 3.1.7
class FourierSeriesReferenceCorrection: public TabulatedMediaReferenceCorrection
{
public:

    FourierSeriesReferenceCorrection( const double startTime,
                                      const double endTime,
                                      const std::vector< double > coefficients );

    double computeReferenceCorrection( const double time ) override;

    std::vector< double > getSineCoefficients( )
    {
        return sineCoefficients_;
    }

    std::vector< double > getCosineCoefficients( )
    {
        return cosineCoefficients_;
    }

private:

    double period_;

    std::vector< double > sineCoefficients_;
    std::vector< double > cosineCoefficients_;

};

class TabulatedMediaReferenceCorrectionManager
{
public:

    TabulatedMediaReferenceCorrectionManager( ):
        isLookupSchemeUpdated_( false )
    { }

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

    double computeMediaCorrection( double time );

private:

    std::vector< double > startTimes_;

    std::vector< double > endTimes_;

    std::vector< std::shared_ptr< TabulatedMediaReferenceCorrection > > correctionVector_;

    bool isLookupSchemeUpdated_;

    //! Lookup scheme to find the nearest correction object start time for a given time
    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;
};

enum TroposphericMappingModel
{
    simplified_chao,
    niell
};

class TroposhericElevationMapping
{
public:

    TroposhericElevationMapping( )
    { }

    virtual ~TroposhericElevationMapping( ){ }

    virtual double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) = 0;

    virtual double computeWetTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) = 0;

protected:

private:

};

// Chao tropospheric mapping, according to Moyer (2000), Eq. 10-8 to 10-10
// Used to map a zenith range correction to a different elevation
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

    double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override
    {
        computeCurrentElevation( transmitterState, receiverState, transmissionTime, receptionTime );
        return troposphericSimplifiedChaoMapping( currentElevation_, true );
    }

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

    double troposphericSimplifiedChaoMapping( const double elevation,
                                              const bool dryCorrection );

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

// Niell tropospheric mapping, according to Moyer (2000), section 10.2.1.3.2
// Used to map a zenith range correction to a different elevation
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

    // Moyer (2000), section 10.2.1.3.2
    double computeWetTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override;

    // Moyer (2000), section 10.2.1.3.2
    double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override;

private:

    // Moyer (2000), section 10.2.1.3.2
    double computeMFunction ( const double a, const double b, const double c, const double elevation );

    // Moyer (2000), section 10.2.1.3.2
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

// Moyer (2000), section 10.2.1
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

    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

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

    Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
    {
        return Eigen::Vector3d::Zero( );
    }

    std::function< double ( double time ) > getDryZenithRangeCorrectionFunction( )
    {
        return dryZenithRangeCorrectionFunction_;
    }

    std::function< double ( double time ) > getWetZenithRangeCorrectionFunction( )
    {
        return wetZenithRangeCorrectionFunction_;
    }

protected:

    // Dry atmosphere zenith range correction (in meters)
    std::function< double ( double time ) > dryZenithRangeCorrectionFunction_;

    // Wet atmosphere zenith range correction (in meters)
    std::function< double ( double time ) > wetZenithRangeCorrectionFunction_;

    std::shared_ptr< TroposhericElevationMapping > elevationMapping_;

    // Boolean indicating whether the correction is for uplink or donwlink (necessary when computing the elevation)
    bool isUplinkCorrection_;
};

// Tabulated tropospheric corrections using DSN data, according to Moyer (2000), section 10.2.1
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

std::function< double ( const double ) > getBeanAndDuttonWaterVaporPartialPressureFunction(
        std::function< double ( const double time ) > relativeHumidity,
        std::function< double ( const double time ) > temperature );

// Saastamaoinen model for tropospheric corrections, according to Estefan and Sovers (1994)
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

    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

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

// Calculator of the vertical total electron content (VTEC) of the ionosphere
class VtecCalculator
{
public:

    VtecCalculator( const double referenceIonosphereHeight ):
        referenceIonosphereHeight_( referenceIonosphereHeight )
    { }

    virtual ~VtecCalculator( ){ }

    // Calculates the VTEC in [m^-2]. (m^-2 = 1e16 TECU)
    virtual double calculateVtec( const double time,
                                  const Eigen::Vector3d subIonosphericPointGeodeticPosition ) = 0;

    double getReferenceIonosphereHeight( )
    {
        return referenceIonosphereHeight_;
    }

private:

    const double referenceIonosphereHeight_;

};

// Computation of the VTEC using the model from Jakowski et al. (2011)
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

    double calculateVtec( const double time,
                          const Eigen::Vector3d subIonosphericPointGeodeticPosition ) override;

private:

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

    // Note 1: Function uses a very simple mapping function, following Moyer (2000). At some point, it might be worth
    //      using other mapping functions, in which case the part of the function where the mapping is executed should
    //      be moved to a new class.
    // Note 2: Currently the STEC is calculated from the VTEC using the algorithm for computing the sub-ionospheric point
    //      described by Moyer (2000). This algorithm fails if the line of sight passes near the poles. For an alternative
    //      method see Olsen (2007), section 5.7.
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

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
