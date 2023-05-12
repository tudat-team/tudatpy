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
 *
 */

#ifndef TUDAT_TABULATEDMEDIACORRECTION_H
#define TUDAT_TABULATEDMEDIACORRECTION_H

#include <cmath>
#include <vector>

#include "tudat/math/interpolators.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

class TabulatedMediaReferenceCorrection
{
public:

    TabulatedMediaReferenceCorrection( const double startTime,
                                       const double endTime ):
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

    bool isTimeValid( const double time )
    {
        if ( time < startTime_ || time > endTime_ )
        {
            throw std::runtime_error(
                    "Error when computing tabulated media reference correction: selected time (" + std::to_string( time ) +
                    ") is outside validity interval (" + std::to_string( startTime_ ) + " to " + std::to_string( endTime_ ) +
                    ")." );
        }

        return true;
    }

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

    TabulatedTroposphericCorrection(
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > dryReferenceCorrectionCalculator,
            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > wetReferenceCorrectionCalculator,
            std::shared_ptr< TroposhericElevationMapping > elevationMapping,
            bool isUplinkCorrection ):
        MappedTroposphericCorrection(
                tabulated_tropospheric,
                elevationMapping,
                isUplinkCorrection,
                std::bind( &TabulatedMediaReferenceCorrectionManager::computeMediaCorrection,
                           dryReferenceCorrectionCalculator, std::placeholders::_1 ),
                std::bind( &TabulatedMediaReferenceCorrectionManager::computeMediaCorrection,
                           wetReferenceCorrectionCalculator, std::placeholders::_1 ) ),
        dryReferenceCorrectionCalculator_( dryReferenceCorrectionCalculator ),
        wetReferenceCorrectionCalculator_( wetReferenceCorrectionCalculator )
    { }

private:

    // Dry atmosphere zenith range correction (in meters)
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > dryReferenceCorrectionCalculator_;

    // Wet atmosphere zenith range correction (in meters)
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > wetReferenceCorrectionCalculator_;
};

// Estefan (1994)
class SaastamoinenTroposphericCorrection: public MappedTroposphericCorrection
{
public:

    SaastamoinenTroposphericCorrection(
            std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction,
            std::function< double ( double time ) > pressureFunction,
            std::function< double ( double time ) > temperatureFunction,
            std::function< double ( double time ) > waterVaporPartialPressureFunction,
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

    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction_;

    std::function< double ( double ) > pressureFunction_;

    std::function< double ( double ) > temperatureFunction_;

    std::function< double ( double ) > waterVaporPartialPressureFunction_;
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

} // namespace observation_models

} // namespace tudat

#endif //TUDATBUNDLE_TABULATEDMEDIACORRECTION_H
