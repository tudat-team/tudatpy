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

    if ( groundStationTime != currentStationTime_ )
    {
        currentElevation_ = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ),
                                                groundStationTime );
        currentStationTime_ = groundStationTime;
    }

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

    return computeMFunction( aDry, bDry, cDry, elevation ) + ( 1.0 / std::sin( elevation ) - computeMFunction(
            aHt, bHt, cHt, elevation ) ) * altitude;

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

    double normalizedTime = ( time - startOfYearTime - 28 * physical_constants::JULIAN_DAY ) / 365.25 * physical_constants::JULIAN_DAY;

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

double TabulatedTroposphericCorrection::calculateLightTimeCorrection(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    double stationTime;
    if ( isUplinkCorrection_ )
    {
        stationTime = transmissionTime;
    }
    else
    {
        stationTime = receptionTime;
    }

    return dryReferenceCorrectionCalculator_->computeMediaCorrection( stationTime ) *
        elevationMapping_->computeDryTroposphericMapping(
                transmitterState, receiverState, transmissionTime, receptionTime ) +
        wetReferenceCorrectionCalculator_->computeMediaCorrection( stationTime ) *
        elevationMapping_->computeWetTroposphericMapping(
                transmitterState, receiverState, transmissionTime, receptionTime );
}

} // namespace observation_models

} // namespace tudat