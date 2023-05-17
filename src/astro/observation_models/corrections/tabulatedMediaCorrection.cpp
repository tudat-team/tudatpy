/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/tabulatedMediaCorrection.h"

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/utilities.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"



namespace tudat
{

namespace observation_models
{

double PowerSeriesReferenceCorrection::computeReferenceCorrection( const double time )
{
    isTimeValid( time );

    const double normalizedTime = 2.0 * ( ( time - startTime_ ) / ( endTime_ - startTime_ ) ) - 1.0;

    double correction = 0;
    for ( unsigned int i = 0; i < coefficients_.size( ); ++i )
    {
        correction += coefficients_.at( i ) * std::pow( normalizedTime, i );
    }

    return correction;
}

FourierSeriesReferenceCorrection::FourierSeriesReferenceCorrection(
        const double startTime,
        const double endTime,
        const std::vector< double > coefficients ):
    TabulatedMediaReferenceCorrection( startTime, endTime )
{
    if ( coefficients.size( ) < 2 || coefficients.size( ) % 2 != 0 )
    {
        throw std::runtime_error(
                "Error when computing Fourier series tabulated media reference correction: size of specified coefficients ("
                + std::to_string( coefficients.size( ) ) + ") is invalid." );
    }

    period_ = coefficients.at( 0 );

    cosineCoefficients_.push_back( coefficients.at( 1 ) );
    sineCoefficients_.push_back( 0.0 );

    for ( unsigned int i = 2; i < coefficients.size( ); i = i + 2 )
    {
        cosineCoefficients_.push_back( coefficients.at( i ) );
        sineCoefficients_.push_back( coefficients.at( i + 1 ) );
    }
}

double FourierSeriesReferenceCorrection::computeReferenceCorrection( const double time )
{
    isTimeValid( time );

    const double normalizedTime = 2.0 * mathematical_constants::PI * ( time - startTime_ ) / period_;

    double correction = 0;
    for ( unsigned int i = 0; i < sineCoefficients_.size( ); ++i )
    {
        correction += cosineCoefficients_.at( i ) * std::cos( i * normalizedTime );
        correction += sineCoefficients_.at( i ) * std::sin( i * normalizedTime );
    }

    return correction;
}

double TabulatedMediaReferenceCorrectionManager::computeMediaCorrection( double time )
{

    if ( correctionVector_.empty( ) )
    {
        throw std::runtime_error("Error when computing reference media correction: no correction object provided. ");
    }

    if ( !isLookupSchemeUpdated_ )
    {
        startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( startTimes_ );
        isLookupSchemeUpdated_ = true;
    }

    int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( time );

    return correctionVector_.at( lowerNearestNeighbour )->computeReferenceCorrection( time );
}

double SimplifiedChaoTroposphericMapping::troposphericSimplifiedChaoMapping( const double elevation,
                                                                             const bool dryCorrection )
{
    double a, b;

    if ( dryCorrection )
    {
        a = 0.00143;
        b = 0.0445;
    }
    else
    {
        a = 0.00035;
        b = 0.017;
    }

    return 1.0 / ( std::sin( elevation ) + a / ( std::tan( elevation ) + b ) );
}

void SimplifiedChaoTroposphericMapping::computeCurrentElevation(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime )
{
    Eigen::Vector6d groundStationState, spacecraftState;
    double groundStationTime;
    if ( isUplinkCorrection_ )
    {
        groundStationState = transmitterState;
        groundStationTime = transmissionTime;
        spacecraftState = receiverState;
    }
    else
    {
        groundStationState = receiverState;
        groundStationTime = receptionTime;
        spacecraftState = transmitterState;
    }

    currentElevation_ = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ),
                                            groundStationTime );
}

NiellTroposphericMapping::NiellTroposphericMapping(
        std::function< double ( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
        std::function< Eigen::Vector3d ( double time ) > groundStationGeodeticPositionFunction,
        bool isUplinkCorrection ):
    TroposhericElevationMapping( ),
    elevationFunction_( elevationFunction ),
    groundStationGeodeticPositionFunction_( groundStationGeodeticPositionFunction ),
    isUplinkCorrection_( isUplinkCorrection )
{
    aDryAverageInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, aDryAverage_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
    bDryAverageInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, bDryAverage_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
    cDryAverageInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, cDryAverage_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );

    aDryAmplitudeInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, aDryAmplitude_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
    bDryAmplitudeInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, bDryAmplitude_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
    cDryAmplitudeInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, cDryAmplitude_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );

    aWetInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, aWet_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
    bWetInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, bWet_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
    cWetInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, cWet_ ),
            interpolators::huntingAlgorithm, interpolators::use_boundary_value );
}

double NiellTroposphericMapping::computeWetTroposphericMapping(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    Eigen::Vector6d groundStationState, spacecraftState;
    double groundStationTime;
    if ( isUplinkCorrection_ )
    {
        groundStationState = transmitterState;
        groundStationTime = transmissionTime;
        spacecraftState = receiverState;
    }
    else
    {
        groundStationState = receiverState;
        groundStationTime = receptionTime;
        spacecraftState = transmitterState;
    }

    double elevation = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ),
                                           groundStationTime );
    Eigen::Vector3d groundStationGeodeticPosition = groundStationGeodeticPositionFunction_( groundStationTime );
    double geodeticLatitude = groundStationGeodeticPosition( 1 );

    double aWetInterpolated = aWetInterpolator_->interpolate( std::abs( geodeticLatitude ) );
    double bWetInterpolated = bWetInterpolator_->interpolate( std::abs( geodeticLatitude ) );
    double cWetInterpolated = cWetInterpolator_->interpolate( std::abs( geodeticLatitude ) );

    return computeMFunction( aWetInterpolated, bWetInterpolated, cWetInterpolated, elevation );
}

double NiellTroposphericMapping::computeDryTroposphericMapping(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    Eigen::Vector6d groundStationState, spacecraftState;
    double groundStationTime;
    if ( isUplinkCorrection_ )
    {
        groundStationState = transmitterState;
        groundStationTime = transmissionTime;
        spacecraftState = receiverState;
    }
    else
    {
        groundStationState = receiverState;
        groundStationTime = receptionTime;
        spacecraftState = transmitterState;
    }

    double elevation = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ),
                                           groundStationTime );
    Eigen::Vector3d groundStationGeodeticPosition = groundStationGeodeticPositionFunction_( groundStationTime );
    double altitude = groundStationGeodeticPosition( 0 );
    double geodeticLatitude = groundStationGeodeticPosition( 1 );

    double aDry = computeDryCoefficient( aDryAverageInterpolator_, aDryAmplitudeInterpolator_, groundStationTime, geodeticLatitude );
    double bDry = computeDryCoefficient( bDryAverageInterpolator_, bDryAmplitudeInterpolator_, groundStationTime, geodeticLatitude );
    double cDry = computeDryCoefficient( cDryAverageInterpolator_, cDryAmplitudeInterpolator_, groundStationTime, geodeticLatitude );

    double altitudeCorrection = 0.0;
    double sinElevation = std::sin( elevation );
    if ( std::abs( sinElevation ) > 1e-12 )
    {
        // Altitude in equation should be in km
        altitudeCorrection = ( 1.0 / sinElevation - computeMFunction( aHt_, bHt_, cHt_, elevation ) ) * altitude * 1e-3;
    }

    return computeMFunction( aDry, bDry, cDry, elevation ) + altitudeCorrection;
}

double NiellTroposphericMapping::computeMFunction ( const double a, const double b, const double c, const double elevation )
{
    double numerator = 1.0 + a / ( 1.0 + b / ( 1.0 + c ) );
    double sinEl = std::sin( elevation );
    double denominator = sinEl + a / ( sinEl + b / ( sinEl + c ) );

    return numerator / denominator;
}

double NiellTroposphericMapping::computeDryCoefficient(
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > averageInterpolator,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > amplitudeInterpolator,
        const double time,
        const double geodeticLatitude )
{
    double coefficientAverage = averageInterpolator->interpolate( std::abs( geodeticLatitude ) );
    double coefficientAmplitude = amplitudeInterpolator->interpolate( std::abs( geodeticLatitude ) );

    // Get calendar date
    int year, month, day;
    double fractionOfDay;
    iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000,
               time / physical_constants::JULIAN_DAY, &year, &month, &day, &fractionOfDay );

    // Get time since start of the calendar year
    double startOfYearTime = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            year, 1, 1, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    double normalizedTime = ( time - startOfYearTime - 28 * physical_constants::JULIAN_DAY ) / ( 365.25 * physical_constants::JULIAN_DAY );

    double dryCoefficient;
    if ( geodeticLatitude >= 0 )
    {
        dryCoefficient = coefficientAverage - coefficientAmplitude * std::cos(
                2.0 * mathematical_constants::PI * normalizedTime );
    }
    else
    {
        dryCoefficient = coefficientAverage - coefficientAmplitude * std::cos(
                2.0 * mathematical_constants::PI * ( normalizedTime + 0.5 ) );
    }

    return dryCoefficient;
}

double MappedTroposphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
{
    // Retrieve state and time of receiver and transmitter
    Eigen::Vector6d transmitterState, receiverState;
    double transmissionTime, receptionTime;
    getTransmissionReceptionTimesAndStates(
            linkEndsStates, linkEndsTimes, currentMultiLegTransmitterIndex, transmitterState, receiverState,
            transmissionTime, receptionTime );

    double stationTime;
    if ( isUplinkCorrection_ )
    {
        stationTime = transmissionTime;
    }
    else
    {
        stationTime = receptionTime;
    }

    return ( dryZenithRangeCorrectionFunction_( stationTime ) *
             elevationMapping_->computeDryTroposphericMapping(
                transmitterState, receiverState, transmissionTime, receptionTime ) +
            wetZenithRangeCorrectionFunction_( stationTime ) *
            elevationMapping_->computeWetTroposphericMapping(
                transmitterState, receiverState, transmissionTime, receptionTime ) ) /
                physical_constants::getSpeedOfLight< double >( );
}

double calculateBeanAndDuttonWaterVaporPartialPressure( double relativeHumidity,
                                                        double temperature )
{
    if ( relativeHumidity < 0 || relativeHumidity > 1 )
    {
        throw std::runtime_error( "Error when computing partial vapor pressure: invalid relative humidity value (" +
            std::to_string( relativeHumidity ) + "), should be in [0,1]." );
    }

    // 1e2 factor is conversion from mbar to Pa
    return 6.11 * relativeHumidity * std::pow( 10.0, 7.5 * ( ( temperature - 273.15 ) / ( temperature - 35.85 ) ) ) * 1e2;
}

std::function< double ( const double ) > getBeanAndDuttonWaterVaporPartialPressureFunction(
        std::function< double ( const double time ) > relativeHumidity,
        std::function< double ( const double time ) > temperature )
{
    return [=] ( double time ) { return calculateBeanAndDuttonWaterVaporPartialPressure(
            relativeHumidity( time ), temperature( time ) ); };
}

double SaastamoinenTroposphericCorrection::computeDryZenithRangeCorrection( const double stationTime )
{
    Eigen::Vector3d stationGeodeticPosition = groundStationGeodeticPositionFunction_( stationTime );
    double altitude = stationGeodeticPosition( 0 );
    double geodeticLatitude = stationGeodeticPosition( 1 );

    double gravitationalAccelerationFactor = 1.0 - 0.00266 * std::cos( 2.0 * geodeticLatitude ) - 2.8e-7 * altitude;

    // 1e2 factor is conversion of pressure from Pa to mBar
    return 0.0022768 * pressureFunction_( stationTime ) * 1e2 / gravitationalAccelerationFactor;
}

double SaastamoinenTroposphericCorrection::computeWetZenithRangeCorrection( const double stationTime )
{
    // 1e2 factor is conversion of pressure from Pa to mBar
    return 0.002277 * waterVaporPartialPressureFunction_( stationTime ) * 1e2 * (
            1255.0 / temperatureFunction_( stationTime ) + 0.05 );
}

TabulatedIonosphericCorrection::TabulatedIonosphericCorrection(
        std::shared_ptr< TabulatedMediaReferenceCorrectionManager > referenceCorrectionCalculator,
        std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction,
        ObservableType observableType,
        bool isUplinkCorrection,
        double referenceFrequency ):
    LightTimeCorrection( tabulated_ionospheric ),
    referenceCorrectionCalculator_( referenceCorrectionCalculator ),
    transmittedFrequencyFunction_( transmittedFrequencyFunction ),
    referenceFrequency_( referenceFrequency ),
    isUplinkCorrection_( isUplinkCorrection )
{
    if ( isRadiometricObservableType( observableType ) )
    {
        if ( isGroupVelocityBasedObservableType( observableType ) )
        {
            sign_ = 1;
        }
        else if ( isPhaseVelocityBasedObservableType( observableType ) )
        {
            sign_ = -1;
        }
        else
        {
            throw std::runtime_error( "Error when creating tabulated ionospheric correction: radiometric correction not "
                                      "recognized." );
        }
    }
    else
    {
        throw std::runtime_error( "Error when creating tabulated ionospheric correction: correction is only valid for "
                                  "radiometric types." );
    }

}

double TabulatedIonosphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
{
    // Retrieve state and time of receiver and transmitter
    Eigen::Vector6d legTransmitterState, legReceiverState;
    double legTransmissionTime, legReceptionTime;
    getTransmissionReceptionTimesAndStates(
            linkEndsStates, linkEndsTimes, currentMultiLegTransmitterIndex, legTransmitterState, legReceiverState,
            legTransmissionTime, legReceptionTime );

    double stationTime;
    if ( isUplinkCorrection_ )
    {
        stationTime = legTransmissionTime;
    }
    else
    {
        stationTime = legReceptionTime;
    }

    double firstLegTransmissionTime = linkEndsTimes.front( );

    // Retrieve frequency bands
    std::vector< FrequencyBands > frequencyBands;
    if( ancillarySettings == nullptr )
    {
        throw std::runtime_error(
                "Error when computing tabulated ionospheric corrections: no ancillary settings found. " );
    }
    try
    {
        frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error(
                "Error when retrieving frequency bands for tabulated ionospheric corrections: " +
                std::string( caughtException.what( ) ) );
    }

    return sign_ * referenceCorrectionCalculator_->computeMediaCorrection( stationTime ) *
        std::pow( referenceFrequency_ / transmittedFrequencyFunction_( frequencyBands, firstLegTransmissionTime ), 2.0 );
}

} // namespace observation_models

} // namespace tudat