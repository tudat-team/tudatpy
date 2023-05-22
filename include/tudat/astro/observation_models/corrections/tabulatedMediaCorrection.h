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

typedef std::map< std::pair< std::string, std::string >, std::map< observation_models::ObservableType,
    std::shared_ptr< observation_models::TabulatedMediaReferenceCorrectionManager > > > AtmosphericCorrectionPerStationAndSpacecraftType;

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

// Moyer (2000), Eq. 10-8 to 10-10
class SimplifiedChaoTroposphericMapping: public TroposhericElevationMapping
{
public:
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

    double currentElevation_;

    std::function< double ( Eigen::Vector3d, double ) > elevationFunction_;

    // Boolean indicating whether the correction is for uplink or donwlink (necessary when computing the elevation)
    bool isUplinkCorrection_;

};

// Moyer (2000), section 10.2.1.3.2
class NiellTroposphericMapping: public TroposhericElevationMapping
{
public:
    NiellTroposphericMapping(
            std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
            std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction,
            bool isUplinkCorrection );

    double computeWetTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override;

    double computeDryTroposphericMapping(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) override;

private:

    double computeMFunction ( const double a, const double b, const double c, const double elevation );

    double computeDryCoefficient(
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > averageInterpolator,
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > amplitudeInterpolator,
            const double time,
            const double geodeticLatitude );

    std::function< double ( Eigen::Vector3d, double ) > elevationFunction_;

    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction_;

    // Boolean indicating whether the correction is for uplink or donwlink (necessary when computing the elevation)
    bool isUplinkCorrection_;

    const std::vector< double > referenceGeodeticLatitudes_ = { 15.0 * mathematical_constants::PI / 180.0,
                                                                30.0 * mathematical_constants::PI / 180.0,
                                                                45.0 * mathematical_constants::PI / 180.0,
                                                                60.0 * mathematical_constants::PI / 180.0,
                                                                75.0 * mathematical_constants::PI / 180.0 };

    const std::vector< double > aDryAverage_ = { 1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3 };
    const std::vector< double > bDryAverage_ = { 2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3 };
    const std::vector< double > cDryAverage_ = { 62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3 };

    const std::vector< double > aDryAmplitude_ = { 0.0e-5, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5 };
    const std::vector< double > bDryAmplitude_ = { 0.0e-5, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5 };
    const std::vector< double > cDryAmplitude_ = { 0.0e-5, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5 };

    const double aHt_ = 2.53e-5;
    const double bHt_ = 5.49e-3;
    const double cHt_ = 1.14e-3;

    const std::vector< double > aWet_ = { 5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4 };
    const std::vector< double > bWet_ = { 1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3 };
    const std::vector< double > cWet_ = { 4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2 };

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

// Moyer (2000), section 10.2.1
class TabulatedTroposphericCorrection: public MappedTroposphericCorrection
{
public:

    /*!
     *
     * @param seasonalModelDryZenithCorrectionCalculator, seasonalModelWetZenithCorrectionCalculator Correction calculators
     *      based on seasonal model. These corrections should either correspond to Figure 3a or 3b of Estefan and Sovers (1994).
     * @param dryZenithCorrectionAdjustmentCorrectionCalculator, wetZenithCorrectionAdjustmentCorrectionCalculator Corrections
     *      to the seasonal model based on real time data. Should be read from DSN TRK-2-23 files.
     * @param elevationMapping
     * @param isUplinkCorrection
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
 * Calculate the partial vapor pressure according to the Bean and Dutton (1966) model, as described by Estefan (1994),
 * Eq. 16.
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

// Estefan (1994)
class SaastamoinenTroposphericCorrection: public MappedTroposphericCorrection
{
public:

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

    std::function< Eigen::Vector3d ( const double ) > groundStationGeodeticPositionFunction_;

    std::function< double ( const double ) > pressureFunction_;

    std::function< double ( const double ) > temperatureFunction_;

    std::function< double ( const double ) > waterVaporPartialPressureFunction_;
};

// Moyer (2000), section 10.2.2
class TabulatedIonosphericCorrection: public LightTimeCorrection
{
public:

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

    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > referenceCorrectionCalculator_;

    std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction_;

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

    ~VtecCalculator( ){ }

    virtual double calculateVtec( const double time,
                                  const Eigen::Vector3d subIonosphericPointGeodeticPosition ) = 0;

    double getReferenceIonosphereHeight( )
    {
        return referenceIonosphereHeight_;
    }

private:

    const double referenceIonosphereHeight_;

};

// Model from Jakowski et al. (2011)
// Computation of geomagnetic latitude from Klobuchar (1975)
class JakowskiVtecCalculator: public VtecCalculator
{
public:

    // Default values copied from GODOT
    // - Geomagnetic pole latitude and longitude:
    //      "Geomagnetic North pole latitude and longitude, values for 2025 from http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
    //      Variation of Geomagnetic poles is historically slow with very minor effect on model results
    //      e.g. (1.5,2.7) deg difference between 1960 and 2010, resulting in <3 cm effect on range"
    // - Reference ionosphere height:
    //      "400km confirmed as good choice by Jakowski in email to G.Bellei, 11.11.2014"
    JakowskiVtecCalculator(
            std::function< double ( const double time ) > sunDeclinationFunction,
            std::function< double ( const double time ) > observedSolarRadioFlux107Function,
            bool useUtcTime = false,
            double geomagneticPoleLatitude = unit_conversions::convertDegreesToRadians( 80.9 ),
            double geomagneticPoleLongitude = unit_conversions::convertDegreesToRadians( -72.6 ),
            double referenceIonosphereHeight = 400.0e3 ):
        VtecCalculator( referenceIonosphereHeight ),
        sunDeclinationFunction_( sunDeclinationFunction ),
        observedSolarRadioFlux107Function_( observedSolarRadioFlux107Function ),
        geomagneticPoleLatitude_( geomagneticPoleLatitude ),
        geomagneticPoleLongitude_( geomagneticPoleLongitude )
    {
        if ( useUtcTime )
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

    // Coefficients of Jakowski model
    const std::vector< double > jakowskiCoefficients_ = { 0.89656, 0.16984, -0.02166, 0.05928, 0.00738, 0.13912,
                                                          -0.17593, -0.34545, 1.1167, 1.1573, -4.3356, 0.17775 };

    const std::function< double ( const double time ) > sunDeclinationFunction_;

    //! Observed (unadjusted) value of F10.7. Expressed in units of 10-22 W/m2/Hz.
    const std::function< double ( const double time ) > observedSolarRadioFlux107Function_;

    const double geomagneticPoleLatitude_;

    const double geomagneticPoleLongitude_;

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter_;

};

// Moyer (2000), section 10.3.1
class MappedVtecIonosphericCorrection: public LightTimeCorrection
{
public:

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

    std::shared_ptr< VtecCalculator > vtecCalculator_;

    std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction_;

    std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction_;

    std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > azimuthFunction_;

    std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction_;

    const double bodyWithAtmosphereMeanEquatorialRadius_;

    const double firstOrderDelayCoefficient_;

    // Sign of the correction (+1 or -1)
    int sign_;

    // Boolean indicating whether the correction is for uplink or downlink
    bool isUplinkCorrection_;

};

} // namespace observation_models

} // namespace tudat

#endif //TUDATBUNDLE_TABULATEDMEDIACORRECTION_H
