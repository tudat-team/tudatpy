/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/atmosphereCorrection.h"

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/utilities.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace observation_models
{

bool TabulatedMediaReferenceCorrection::isTimeValid( const double time )
{
    /*
    if( !std::isnan( startTime_ ) && time < startTime_ )
    {
        throw std::runtime_error( "Error when computing tabulated media reference correction: selected time (" + std::to_string( time ) +
                                  ") is below start time (" + std::to_string( startTime_ ) + ")." );
    }

    if( !std::isnan( endTime_ ) && time > endTime_ )
    {
        throw std::runtime_error( "Error when computing tabulated media reference correction: selected time (" + std::to_string( time ) +
                                  ") is over end time (" + std::to_string( endTime_ ) + ")." );
    }
    */

    ///// Warnings instead of error:
    if( !std::isnan( startTime_ ) && time < startTime_ )

    {
        std::cerr << "Warning when computing tabulated media reference correction: selected time (" + std::to_string( time ) +
                        ") is below start time (" + std::to_string( startTime_ ) + "). Applying 0 correction."
                  << std::endl;

        return false;
    }

    if( !std::isnan( endTime_ ) && time > endTime_ )

    {
        std::cerr << "Warning when computing tabulated media reference correction: selected time (" + std::to_string( time ) +
                        ") is over end time (" + std::to_string( endTime_ ) + "). Applying 0 correction."
                  << std::endl;

        return false;
    }

    return true;
}

double PowerSeriesReferenceCorrection::computeReferenceCorrection( const double time )
{
    bool timeCover = isTimeValid( time );
    const double normalizedTime = 2.0 * ( ( time - startTime_ ) / ( endTime_ - startTime_ ) ) - 1.0;
    double correction = 0;

    if( timeCover == true )

    {
        for( unsigned int i = 0; i < coefficients_.size( ); ++i )
        {
            correction += coefficients_.at( i ) * std::pow( normalizedTime, i );
        }
    }
    // else correction 0

    return correction;
}

FourierSeriesReferenceCorrection::FourierSeriesReferenceCorrection( const double startTime,
                                                                    const double endTime,
                                                                    const std::vector< double > coefficients ):
    TabulatedMediaReferenceCorrection( startTime, endTime )
{
    if( coefficients.size( ) < 2 || coefficients.size( ) % 2 != 0 )
    {
        throw std::runtime_error(
                "Error when computing Fourier series tabulated media reference correction: size of specified coefficients (" +
                std::to_string( coefficients.size( ) ) + ") is invalid." );
    }

    period_ = coefficients.at( 0 );

    cosineCoefficients_.push_back( coefficients.at( 1 ) );
    sineCoefficients_.push_back( 0.0 );

    for( unsigned int i = 2; i < coefficients.size( ); i = i + 2 )
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
    for( unsigned int i = 0; i < sineCoefficients_.size( ); ++i )
    {
        correction += cosineCoefficients_.at( i ) * std::cos( i * normalizedTime );
        correction += sineCoefficients_.at( i ) * std::sin( i * normalizedTime );
    }

    return correction;
}

double TabulatedMediaReferenceCorrectionManager::computeMediaCorrection( double time )
{
    if( correctionVector_.empty( ) )
    {
        throw std::runtime_error( "Error when computing reference media correction: no correction object provided. " );
    }

    int lowerNearestNeighbour = 0;
    if( startTimes_.size( ) > 1 )
    {
        if( !isLookupSchemeUpdated_ )
        {
            startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( startTimes_ );
            isLookupSchemeUpdated_ = true;
        }
        lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( time );
    }

    return correctionVector_.at( lowerNearestNeighbour )->computeReferenceCorrection( time );
}

double SimplifiedChaoTroposphericMapping::troposphericSimplifiedChaoMapping( const double elevation, const bool dryCorrection )
{
    double a, b;

    // Moyer (2000), eq. 10-9
    if( dryCorrection )
    {
        a = 0.00143;
        b = 0.0445;
    }
    // Moyer (2000), eq. 10-10
    else
    {
        a = 0.00035;
        b = 0.017;
    }

    // Moyer (2000), eq. 10-8
    return 1.0 / ( std::sin( elevation ) + a / ( std::tan( elevation ) + b ) );
}

void SimplifiedChaoTroposphericMapping::computeCurrentElevation( const Eigen::Vector6d& transmitterState,
                                                                 const Eigen::Vector6d& receiverState,
                                                                 const double transmissionTime,
                                                                 const double receptionTime )
{
    Eigen::Vector6d groundStationState, spacecraftState;
    double groundStationTime;
    if( isUplinkCorrection_ )
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

    currentElevation_ = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ), groundStationTime );
}

NiellTroposphericMapping::NiellTroposphericMapping(
        std::function< double( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
        std::function< Eigen::Vector3d( double time ) > groundStationGeodeticPositionFunction,
        bool isUplinkCorrection ):
    TroposhericElevationMapping( ), elevationFunction_( elevationFunction ),
    groundStationGeodeticPositionFunction_( groundStationGeodeticPositionFunction ), isUplinkCorrection_( isUplinkCorrection )
{
    aDryAverageInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, aDryAverage_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
    bDryAverageInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, bDryAverage_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
    cDryAverageInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, cDryAverage_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );

    aDryAmplitudeInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, aDryAmplitude_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
    bDryAmplitudeInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, bDryAmplitude_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
    cDryAmplitudeInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, cDryAmplitude_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );

    aWetInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, aWet_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
    bWetInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, bWet_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
    cWetInterpolator_ = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            utilities::createMapFromVectors( referenceGeodeticLatitudes_, cWet_ ),
            interpolators::huntingAlgorithm,
            interpolators::use_boundary_value );
}

double NiellTroposphericMapping::computeWetTroposphericMapping( const Eigen::Vector6d& transmitterState,
                                                                const Eigen::Vector6d& receiverState,
                                                                const double transmissionTime,
                                                                const double receptionTime )
{
    Eigen::Vector6d groundStationState, spacecraftState;
    double groundStationTime;
    if( isUplinkCorrection_ )
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

    double elevation = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ), groundStationTime );
    Eigen::Vector3d groundStationGeodeticPosition = groundStationGeodeticPositionFunction_( groundStationTime );
    double geodeticLatitude = groundStationGeodeticPosition( 1 );

    double aWetInterpolated = aWetInterpolator_->interpolate( std::abs( geodeticLatitude ) );
    double bWetInterpolated = bWetInterpolator_->interpolate( std::abs( geodeticLatitude ) );
    double cWetInterpolated = cWetInterpolator_->interpolate( std::abs( geodeticLatitude ) );

    return computeMFunction( aWetInterpolated, bWetInterpolated, cWetInterpolated, elevation );
}

double NiellTroposphericMapping::computeDryTroposphericMapping( const Eigen::Vector6d& transmitterState,
                                                                const Eigen::Vector6d& receiverState,
                                                                const double transmissionTime,
                                                                const double receptionTime )
{
    Eigen::Vector6d groundStationState, spacecraftState;
    double groundStationTime;
    if( isUplinkCorrection_ )
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

    double elevation = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ), groundStationTime );
    Eigen::Vector3d groundStationGeodeticPosition = groundStationGeodeticPositionFunction_( groundStationTime );
    double altitude = groundStationGeodeticPosition( 0 );
    double geodeticLatitude = groundStationGeodeticPosition( 1 );

    double aDry = computeDryCoefficient( aDryAverageInterpolator_, aDryAmplitudeInterpolator_, groundStationTime, geodeticLatitude );
    double bDry = computeDryCoefficient( bDryAverageInterpolator_, bDryAmplitudeInterpolator_, groundStationTime, geodeticLatitude );
    double cDry = computeDryCoefficient( cDryAverageInterpolator_, cDryAmplitudeInterpolator_, groundStationTime, geodeticLatitude );

    double altitudeCorrection = 0.0;
    double sinElevation = std::sin( elevation );
    if( std::abs( sinElevation ) > 1e-12 )
    {
        // Altitude in equation should be in km
        altitudeCorrection = ( 1.0 / sinElevation - computeMFunction( aHt_, bHt_, cHt_, elevation ) ) * altitude * 1e-3;
    }

    return computeMFunction( aDry, bDry, cDry, elevation ) + altitudeCorrection;
}

double NiellTroposphericMapping::computeMFunction( const double a, const double b, const double c, const double elevation )
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

    double normalizedTime = ( sofa_interface::convertSecondsSinceEpochToSecondsOfYear( time ) - 28.0 * physical_constants::JULIAN_DAY ) /
            ( 365.25 * physical_constants::JULIAN_DAY );

    double dryCoefficient;
    if( geodeticLatitude >= 0 )
    {
        dryCoefficient = coefficientAverage - coefficientAmplitude * std::cos( 2.0 * mathematical_constants::PI * normalizedTime );
    }
    else
    {
        dryCoefficient =
                coefficientAverage - coefficientAmplitude * std::cos( 2.0 * mathematical_constants::PI * ( normalizedTime + 0.5 ) );
    }

    return dryCoefficient;
}

MappedTroposphericCorrection::MappedTroposphericCorrection( const LightTimeCorrectionType lightTimeCorrectionType,
                                                            std::shared_ptr< TroposhericElevationMapping > elevationMapping,
                                                            bool isUplinkCorrection,
                                                            std::function< double( double time ) > dryZenithRangeCorrectionFunction,
                                                            std::function< double( double time ) > wetZenithRangeCorrectionFunction ):
    LightTimeCorrection( lightTimeCorrectionType ), dryZenithRangeCorrectionFunction_( dryZenithRangeCorrectionFunction ),
    wetZenithRangeCorrectionFunction_( wetZenithRangeCorrectionFunction ), elevationMapping_( elevationMapping ),
    isUplinkCorrection_( isUplinkCorrection )
{ }

double MappedTroposphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
{
    // Retrieve state and time of receiver and transmitter
    Eigen::Vector6d transmitterState, receiverState;
    double transmissionTime, receptionTime;
    getTransmissionReceptionTimesAndStates( linkEndsStates,
                                            linkEndsTimes,
                                            currentMultiLegTransmitterIndex,
                                            transmitterState,
                                            receiverState,
                                            transmissionTime,
                                            receptionTime );

    double stationTime;
    if( isUplinkCorrection_ )
    {
        stationTime = transmissionTime;
    }
    else
    {
        stationTime = receptionTime;
    }

    // Moyer (2000), eq. 10-1
    double delay =
            ( dryZenithRangeCorrectionFunction_( stationTime ) *
                      elevationMapping_->computeDryTroposphericMapping( transmitterState, receiverState, transmissionTime, receptionTime ) +
              wetZenithRangeCorrectionFunction_( stationTime ) *
                      elevationMapping_->computeWetTroposphericMapping(
                              transmitterState, receiverState, transmissionTime, receptionTime ) ) /
            physical_constants::getSpeedOfLight< double >( );
    return delay;
}

// double VMF1TroposphericCorrection::getDryMappingFunctions(
//     const double dryACoefficient,
//     const double wetACoefficient,
//     const double elevationAngle,
//     const double stationLatitude,
//     const double currentModifiedJulianDay )
//{
//
//     double bh = 0.0029;
//     double c0h = 0.062;
//     double c11h, c10h, quadrant;
//
//     if (stationLatitude < 0 )
//     {
//         quadrant = mathematical_constants::PI;
//         c11h = 0.007;
//         c10h = 0.002;
//     }
//     else
//     {
//         quadrant = 0.0;
//         c11h = 0.005;
//         c10h = 0.001;
//     }
//
//     double dayOfYear = currentModifiedJulianDay  - 44239.0 + 1.0 - 28.0;
//     double ch = c0h + ( ( std::cos( dayOfYear/365.25 * 2.0 * mathematical_constants::PI + quadrant ) + 1.0 )* c11h / 2.0 + c10h )*
//         (1.0 - std::cos( stationLatitude ) );
//
//     double sineElevation   = std::sin( elevationAngle );
//     double beta   = bh / ( sineElevation + ch  );
//     double gamma  = ah/( sineElevation + beta);
//     double topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)));
//     return topcon / ( sineElevation+gamma );
// }

// double VMF1TroposphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
//     const std::vector< Eigen::Vector6d >& linkEndsStates,
//     const std::vector< double >& linkEndsTimes,
//     const unsigned int currentMultiLegTransmitterIndex,
//     const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
//{
//
//     // Retrieve state and time of receiver and transmitter
//     Eigen::Vector6d transmitterState, receiverState;
//     double transmissionTime, receptionTime;
//     getTransmissionReceptionTimesAndStates( linkEndsStates,
//                                             linkEndsTimes,
//                                             currentMultiLegTransmitterIndex,
//                                             transmitterState,
//                                             receiverState,
//                                             transmissionTime,
//                                             receptionTime );
//
//     double stationTime;
//     if( isUplinkCorrection_ )
//     {
//         stationTime = transmissionTime;
//     }
//     else
//     {
//         stationTime = receptionTime;
//     }
//
//     Eigen::Vector2d zenithDelays = troposphereData_->getZenithDelay( stationTime );
//     Eigen::Vector2d mappingFunctionParameters = troposphereData_->getMappingFunction( stationTime );
//     Eigen::Vector4d gradientParameters = troposphereData_->getGradient( stationTime );
// }


double VMF3TroposphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
    const std::vector< Eigen::Vector6d >& linkEndsStates,
    const std::vector< double >& linkEndsTimes,
    const unsigned int currentMultiLegTransmitterIndex,
    const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
{
    // Extract transmitter/receiver states and times
    Eigen::Vector6d transmitterState, receiverState;
    double transmissionTime, receptionTime;
    getTransmissionReceptionTimesAndStates( linkEndsStates, linkEndsTimes,
                                            currentMultiLegTransmitterIndex,
                                            transmitterState, receiverState,
                                            transmissionTime, receptionTime );

    // Get station time (uplink or downlink)
    double stationTime = ( isUplinkCorrection_ ? transmissionTime : receptionTime );

    // Extract mapping coefficients from station data
    Eigen::Vector2d mappingCoefficients = troposphereData_->getMappingFunction( stationTime );

    // Downcast to VMF3MappingModel and update coefficients
    std::shared_ptr< VMF3MappingModel > vmf3Mapping = std::dynamic_pointer_cast< VMF3MappingModel >( elevationMapping_ );
    if( vmf3Mapping == nullptr )
    {
        throw std::runtime_error( "Error: VMF3MappingModel dynamic cast failed in calculateLightTimeCorrectionWithMultiLegLinkEndStates." );
    }
    vmf3Mapping->updateMappingCoefficients( mappingCoefficients );

    // Compute mapping contributions
    double dryMapping = vmf3Mapping->computeDryTroposphericMapping(
        transmitterState, receiverState, transmissionTime, receptionTime );

    double wetMapping = vmf3Mapping->computeWetTroposphericMapping(
        transmitterState, receiverState, transmissionTime, receptionTime );

    // Get zenith delays
    double dryZenith = computeDryZenithRangeCorrection( stationTime );
    double wetZenith = computeWetZenithRangeCorrection( stationTime );

    // Base correction
    double correction = dryZenith * dryMapping + wetZenith * wetMapping;

    // Gradient correction (optional)
    if ( useGradient_ )
    {
        try
        {
            Eigen::Vector4d gradients = troposphereData_->getGradient( stationTime );
            double gradientContribution = std::dynamic_pointer_cast< VMF3MappingModel >( elevationMapping_ )
                ->computeGradientContribution( transmitterState, receiverState, transmissionTime, receptionTime, gradients );
            correction += gradientContribution;
        }
        catch ( const std::exception& )
        {
            std::cerr << "Warning: VMF3 gradient correction requested but not available. Ignoring gradients." << std::endl;
        }
    }

    return correction / physical_constants::getSpeedOfLight< double >( );
}

double VMF3MappingModel::computeDryTroposphericMapping(
    const Eigen::Vector6d& transmitterState,
    const Eigen::Vector6d& receiverState,
    const double transmissionTime,
    const double receptionTime )
{
    computeCurrentVMFdata( transmitterState, receiverState, transmissionTime, receptionTime );
    return computeMappingFunction( currentDryMappingCoefficient_, true );
}

double VMF3MappingModel::computeWetTroposphericMapping(
    const Eigen::Vector6d& transmitterState,
    const Eigen::Vector6d& receiverState,
    const double transmissionTime,
    const double receptionTime )
{
    computeCurrentVMFdata( transmitterState, receiverState, transmissionTime, receptionTime );
    return computeMappingFunction( currentWetMappingCoefficient_, false );
}


void VMF3MappingModel::computeCurrentVMFdata(
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

    Eigen::Vector3d relativeVector = spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 );
    currentElevation_ = elevationFunction_( relativeVector, groundStationTime );
    currentAzimuth_ = azimuthFunction_( relativeVector, groundStationTime );

    Eigen::Vector3d geodetic = groundStationGeodeticPositionFunction_( groundStationTime );
    currentStationLatitude_ = geodetic( 1 );
    currentStationLongitude_ = geodetic( 2 );

    // Seasonal variation: F2
    double utcTime = timeScaleConverter_->getCurrentTime( basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, groundStationTime );
    currentDayOfYear_ = sofa_interface::convertSecondsSinceEpochToSecondsOfYear( utcTime ) / physical_constants::JULIAN_DAY;
}

double VMF3MappingModel::computeMappingFunction(
    const double mappingCoefficient,
    const bool isHydrostatic ) const
{
    const double a = mappingCoefficient;
    const double bh = evaluateSeasonalCoefficient(
        isHydrostatic ? anm_bh_.A0 : anm_bw_.A0,
        isHydrostatic ? anm_bh_.A1 : anm_bw_.A1,
        isHydrostatic ? anm_bh_.B1 : anm_bw_.B1,
        isHydrostatic ? anm_bh_.A2 : anm_bw_.A2,
        isHydrostatic ? anm_bh_.B2 : anm_bw_.B2,
        computeVnmWnmMatrix( 12, currentStationLatitude_, currentStationLongitude_ ).V,
        computeVnmWnmMatrix( 12, currentStationLatitude_, currentStationLongitude_ ).W,
        currentDayOfYear_ );

    const double ch = evaluateSeasonalCoefficient(
        isHydrostatic ? anm_ch_.A0 : anm_cw_.A0,
        isHydrostatic ? anm_ch_.A1 : anm_cw_.A1,
        isHydrostatic ? anm_ch_.B1 : anm_cw_.B1,
        isHydrostatic ? anm_ch_.A2 : anm_cw_.A2,
        isHydrostatic ? anm_ch_.B2 : anm_cw_.B2,
        computeVnmWnmMatrix( 12, currentStationLatitude_, currentStationLongitude_ ).V,
        computeVnmWnmMatrix( 12, currentStationLatitude_, currentStationLongitude_ ).W,
        currentDayOfYear_ );

    return ( 1.0 + a / ( 1.0 + bh / ( 1.0 + ch ) ) ) /
           ( std::sin( currentElevation_ ) + a / ( std::sin( currentElevation_ ) + bh / ( std::sin( currentElevation_ ) + ch ) ) );
}

double VMF3MappingModel::computeGradientContribution(
    const Eigen::Vector6d& transmitterState,
    const Eigen::Vector6d& receiverState,
    const double transmissionTime,
    const double receptionTime,
    const Eigen::Vector4d& gradients )
{
    computeCurrentVMFdata( transmitterState, receiverState, transmissionTime, receptionTime );
    double sinEl = std::sin( currentElevation_ );
    double tanEl = std::tan( currentElevation_ );
    double sintan = sinEl * tanEl;

    double mgh = 1.0 / ( sintan + 0.0031 );
    double mgw = 1.0 / ( sintan + 0.0007 );

    const double gnh = gradients( 0 );
    const double geh = gradients( 1 );
    const double gnw = gradients( 2 );
    const double gew = gradients( 3 );

    return mgh * ( gnh * std::cos( currentAzimuth_ ) + geh * std::sin( currentAzimuth_ ) ) +
           mgw * ( gnw * std::cos( currentAzimuth_ ) + gew * std::sin( currentAzimuth_ ) );
}



VMF3MappingModel::Vmf3SphericalHarmonicComponentSet VMF3MappingModel::loadCoefficientSet( const std::string& filePath )
{
    Eigen::MatrixXd raw = input_output::readMatrixFromFile( filePath, " " );

    if ( raw.cols( ) < 5 )
    {
        throw std::runtime_error( "VMF3 file must have 5 columns: A0, A1, B1, A2, B2" );
    }

    Vmf3SphericalHarmonicComponentSet set;
    set.A0 = raw.col( 0 );
    set.A1 = raw.col( 1 );
    set.B1 = raw.col( 2 );
    set.A2 = raw.col( 3 );
    set.B2 = raw.col( 4 );
    return set;
}

void VMF3MappingModel::loadLegendreCoefficientTables( )
{
    std::string path = paths::getAtmosphereTablesPath() + "/vmf3/";
    anm_bh_ = loadCoefficientSet( path + "anm_bh.txt" );
    bnm_bh_ = loadCoefficientSet( path + "bnm_bh.txt" );
    anm_bw_ = loadCoefficientSet( path + "anm_bw.txt" );
    bnm_bw_ = loadCoefficientSet( path + "bnm_bw.txt" );
    anm_ch_ = loadCoefficientSet( path + "anm_ch.txt" );
    bnm_ch_ = loadCoefficientSet( path + "bnm_ch.txt" );
    anm_cw_ = loadCoefficientSet( path + "anm_cw.txt" );
    bnm_cw_ = loadCoefficientSet( path + "bnm_cw.txt" );
}

VMF3MappingModel::VnmWnmMatrix VMF3MappingModel::computeVnmWnmMatrix(
    int nMax, double latitude, double longitude ) const
{
    double theta = mathematical_constants::PI / 2.0 - latitude;
    double x = std::sin( theta ) * std::cos( longitude );
    double y = std::sin( theta ) * std::sin( longitude );
    double z = std::cos( theta );
    std::vector< std::vector< double > > V( nMax + 2, std::vector< double >( nMax + 2, 0.0 ) );
    std::vector< std::vector< double > > W( nMax + 2, std::vector< double >( nMax + 2, 0.0 ) );
    // Base cases
    V[0][0] = 1.0;
    V[1][0] = z * V[0][0];  // V[1][0] = z
    // Fix: manually set V[2][0] = z^2
    V[2][0] = z * V[1][0];
    // Now apply recurrence for n >= 2, m = 0
    for ( int n = 2; n <= nMax; ++n )
    {
        V[n + 1][0] = ( ( 2 * n - 1 ) * z * V[n][0] - ( n - 1 ) * V[n - 1][0] ) / n;
    }
    // Loop over orders m > 0
    for ( int m = 1; m <= nMax; ++m )
    {
        V[m][m] = (2 * m - 1) * ( x * V[m - 1][m - 1] - y * W[m - 1][m - 1] );
        W[m][m] = (2 * m - 1) * ( x * W[m - 1][m - 1] + y * V[m - 1][m - 1] );
        if ( m < nMax )
        {
            V[m + 1][m] = (2 * m + 1) * z * V[m][m];
            W[m + 1][m] = (2 * m + 1) * z * W[m][m];
        }
        for ( int n = m + 2; n <= nMax; ++n )
        {
            V[n + 1][m] = ( (2 * n - 1) * z * V[n][m] - (n + m - 1) * V[n - 1][m] ) / ( n - m );
            W[n + 1][m] = ( (2 * n - 1) * z * W[n][m] - (n + m - 1) * W[n - 1][m] ) / ( n - m );
        }
    }
    return { V, W };
}

double VMF3MappingModel::evaluateSphericalExpansion(
    const Eigen::VectorXd& anm_column,
    const Eigen::VectorXd& bnm_column,
    const std::vector<std::vector<double>>& V,
    const std::vector<std::vector<double>>& W ) const
{
    double result = 0.0;
    int i = 0;
    int nMax = 12;
    for ( int n = 0; n <= nMax; ++n )
    {
        for ( int m = 0; m <= n; ++m, ++i )
        {
            result += anm_column( i ) * V[n][m] + bnm_column( i ) * W[n][m];
        }
    }
    return result;
}

double VMF3MappingModel::evaluateSeasonalCoefficient(
    const Eigen::VectorXd& A0,
    const Eigen::VectorXd& A1,
    const Eigen::VectorXd& B1,
    const Eigen::VectorXd& A2,
    const Eigen::VectorXd& B2,
    const std::vector<std::vector<double>>& V,
    const std::vector<std::vector<double>>& W,
    const double dayOfYear ) const
{
    const double omega1 = 2.0 * mathematical_constants::PI * dayOfYear / 365.25;
    const double omega2 = 2.0 * omega1;
    double sumA0 = evaluateSphericalExpansion( A0, B2, V, W );
    double sumA1 = evaluateSphericalExpansion( A1, B1, V, W );
    double sumA2 = evaluateSphericalExpansion( A2, B2, V, W );
    return sumA0 +
           sumA1 * std::cos( omega1 ) +
           sumA1 * std::sin( omega1 ) +
           sumA2 * std::cos( omega2 ) +
           sumA2 * std::sin( omega2 );
}

double SaastamoinenTroposphericCorrection::computeDryZenithRangeCorrection( const double stationTime )
{
    Eigen::Vector3d stationGeodeticPosition = groundStationGeodeticPositionFunction_( stationTime );
    double altitude = stationGeodeticPosition( 0 );
    double geodeticLatitude = stationGeodeticPosition( 1 );

    // Estefan and Sovers (1994), eq. 13
    double gravitationalAccelerationFactor = 1.0 - 0.00266 * std::cos( 2.0 * geodeticLatitude ) - 2.8e-7 * altitude;

    // Estefan and Sovers (1994), eq. 12
    // 1e-2 factor is conversion of pressure from Pa to mBar
    return 0.0022768 * pressureFunction_( stationTime ) * 1e-2 / gravitationalAccelerationFactor;
}

double SaastamoinenTroposphericCorrection::computeWetZenithRangeCorrection( const double stationTime )
{
    // Estefan and Sovers (1994), eq. 18
    // 1e-2 factor is conversion of pressure from Pa to mBar
    return 0.002277 * waterVaporPartialPressureFunction_( stationTime ) * 1e-2 * ( 1255.0 / temperatureFunction_( stationTime ) + 0.05 );
}

TabulatedIonosphericCorrection::TabulatedIonosphericCorrection(
        std::shared_ptr< TabulatedMediaReferenceCorrectionManager > referenceCorrectionCalculator,
        ObservableType baseObservableType,
        bool isUplinkCorrection,
        double referenceFrequency ):
    LightTimeCorrection( tabulated_ionospheric ), referenceCorrectionCalculator_( referenceCorrectionCalculator ),
    referenceFrequency_( referenceFrequency ),
    isUplinkCorrection_( isUplinkCorrection )
{
    if( isRadiometricObservableType( baseObservableType ) )
    {
        // Moyer (2000), section 10.2.2, 4th paragraph
        if( isGroupVelocityBasedObservableType( baseObservableType ) )
        {
            sign_ = 1;
        }
        else if( isPhaseVelocityBasedObservableType( baseObservableType ) )
        {
            sign_ = -1;
        }
        else
        {
            throw std::runtime_error(
                    "Error when creating tabulated ionospheric correction: radiometric correction not "
                    "recognized." );
        }
    }
    else
    {
        throw std::runtime_error(
                "Error when creating tabulated ionospheric correction: correction is only valid for "
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
    getTransmissionReceptionTimesAndStates( linkEndsStates,
                                            linkEndsTimes,
                                            currentMultiLegTransmitterIndex,
                                            legTransmitterState,
                                            legReceiverState,
                                            legTransmissionTime,
                                            legReceptionTime );

    // Compute light-time correction
    double stationTime = TUDAT_NAN;
    double lightTimeCorrection = 0.0;
    double currentFrequency = TUDAT_NAN;
    if( isUplinkCorrection_ )
    {
        stationTime = legTransmissionTime;
        currentFrequency = ancillarySettings->getIntermediateDoubleData( transmitter_frequency_intermediate, true );
    }
    else
    {
        stationTime = legReceptionTime;
        currentFrequency = ancillarySettings->getIntermediateDoubleData( received_frequency_intermediate, true );
    }
    double stationTimeUtc = sofa_interface::convertTTtoUTC( stationTime );

    if( !std::isnan( currentFrequency ) )
    {
        lightTimeCorrection = ( sign_ * referenceCorrectionCalculator_->computeMediaCorrection( stationTimeUtc ) *
                                std::pow( referenceFrequency_ / currentFrequency, 2.0 ) ) /
                physical_constants::getSpeedOfLight< double >( );
    }

    // Moyer (2000), eq. 10-11
    return lightTimeCorrection;
}

double JakowskiVtecCalculator::calculateVtec( const double time, const Eigen::Vector3d subIonosphericPointGeodeticPosition )
{
    const double subIonosphericLatitude = subIonosphericPointGeodeticPosition( 1 );
    const double subIonosphericLongitude = subIonosphericPointGeodeticPosition( 2 );
    const double sunDeclination = sunDeclinationFunction_( time );

    // Dependency on solar zenith: F1

    // Moyer (2000), eqs. 10-46, 10.47
    const double localTimeHours = std::fmod(
            getUtcTime( time ) / 3600.0 + 12.0 + unit_conversions::convertRadiansToDegrees( subIonosphericLongitude ) / 15.0, 24.0 );

    // Jakowski et al. (2011), eqs. 5, 6, 7
    const double diurnalVariation = 2.0 * mathematical_constants::PI * ( localTimeHours - 14.0 ) / 24.0;
    const double semiDiurnalVariation = 2.0 * mathematical_constants::PI * localTimeHours / 12.0;
    const double terDiurnalVariation = 2.0 * mathematical_constants::PI * localTimeHours / 8.0;

    // Jakowski et al. (2011), eqs. 9, 10, 11
    const double cosChiX = std::cos( subIonosphericLatitude - sunDeclination );
    const double cosChiXx = cosChiX - 2.0 / mathematical_constants::PI * subIonosphericLatitude * std::sin( sunDeclination );
    // From GODOT:
    // "Note: in the paper by Jakowsky (2011) the equation 11 for cxx is wrong (has a square root around cx+0.4)"
    const double cosChiXxx = cosChiX + 0.4;

    // Jakowski et al. (2011), eqs. 8
    const double f1 = cosChiXxx +
            cosChiXx *
                    ( jakowskiCoefficients_.at( 0 ) * std::cos( diurnalVariation ) +
                      jakowskiCoefficients_.at( 1 ) * std::cos( semiDiurnalVariation ) +
                      jakowskiCoefficients_.at( 2 ) * std::sin( semiDiurnalVariation ) +
                      jakowskiCoefficients_.at( 3 ) * std::cos( terDiurnalVariation ) +
                      jakowskiCoefficients_.at( 4 ) * std::sin( terDiurnalVariation ) );

    // Seasonal variation: F2
    const double dayOfYear = sofa_interface::convertSecondsSinceEpochToSecondsOfYear( time ) / physical_constants::JULIAN_DAY;

    // Jakowski et al. (2011), eqs. 13, 14
    const double annualVariation = 2.0 * mathematical_constants::PI * ( dayOfYear - 18.0 ) / 365.25;
    const double semiAnnualVariation = 4.0 * mathematical_constants::PI * ( dayOfYear - 6.0 ) / 365.25;

    // Jakowski et al. (2011), eq. 12
    const double f2 = ( 1.0 + jakowskiCoefficients_.at( 5 ) * std::cos( annualVariation ) +
                        jakowskiCoefficients_.at( 6 ) * std::cos( semiAnnualVariation ) );

    // Geomagnetic field dependency: F3 and F4

    // Geomagnetic latitude, Klobuchar (1975), section 5.3
    const double geomagneticLatitude = std::asin( std::sin( subIonosphericLatitude ) * std::sin( geomagneticPoleLatitude_ ) +
                                                  std::cos( subIonosphericLatitude ) * std::cos( geomagneticPoleLatitude_ ) *
                                                          std::cos( subIonosphericLongitude - geomagneticPoleLongitude_ ) );

    // Jakowski et al. (2011), eq. 15
    const double f3 = 1.0 + jakowskiCoefficients_.at( 7 ) * std::cos( geomagneticLatitude );

    // Jakowski et al. (2011), eqs. 17, 18
    const double phiC1 = unit_conversions::convertDegreesToRadians( 16.0 );
    const double phiC2 = -unit_conversions::convertDegreesToRadians( 10.0 );
    const double sigmaC1 = unit_conversions::convertDegreesToRadians( 12.0 );
    const double sigmaC2 = unit_conversions::convertDegreesToRadians( 13.0 );

    const double ec1 = -std::pow( ( geomagneticLatitude - phiC1 ) / sigmaC1, 2.0 ) / 2.0;
    const double ec2 = -std::pow( ( geomagneticLatitude - phiC2 ) / sigmaC2, 2.0 ) / 2.0;
    // Jakowski et al. (2011), eqs. 16
    const double f4 = ( 1.0 + jakowskiCoefficients_.at( 8 ) * std::exp( ec1 ) + jakowskiCoefficients_.at( 9 ) * std::exp( ec2 ) );

    // Solar activity dependency: F5
    // Jakowski et al. (2011), eq. 19
    const double f5 = jakowskiCoefficients_.at( 10 ) + jakowskiCoefficients_.at( 11 ) * observedSolarRadioFlux107Function_( time );

    // Jakowski et al. (2011), eq. 4
    // 1e16 is conversion factor from TECU to m^-2
    return f1 * f2 * f3 * f4 * f5 * 1e16;
}

MappedVtecIonosphericCorrection::MappedVtecIonosphericCorrection(
        std::shared_ptr< VtecCalculator > vtecCalculator,
        std::function< double( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > elevationFunction,
        std::function< double( Eigen::Vector3d inertialVectorAwayFromStation, double time ) > azimuthFunction,
        std::function< Eigen::Vector3d( double time ) > groundStationGeodeticPositionFunction,
        ObservableType baseObservableType,
        bool isUplinkCorrection,
        double bodyWithAtmosphereMeanEquatorialRadius,
        LightTimeCorrectionType correctionType,
        double firstOrderDelayCoefficient ):
    LightTimeCorrection( correctionType ), vtecCalculator_( vtecCalculator ),
    elevationFunction_( elevationFunction ),
    azimuthFunction_( azimuthFunction ), groundStationGeodeticPositionFunction_( groundStationGeodeticPositionFunction ),
    bodyWithAtmosphereMeanEquatorialRadius_( bodyWithAtmosphereMeanEquatorialRadius ),
    firstOrderDelayCoefficient_( firstOrderDelayCoefficient ), isUplinkCorrection_( isUplinkCorrection )
{
    if( isRadiometricObservableType( baseObservableType ) )
    {
        if( isGroupVelocityBasedObservableType( baseObservableType ) )
        {
            sign_ = 1;
        }
        else if( isPhaseVelocityBasedObservableType( baseObservableType ) )
        {
            sign_ = -1;
        }
        else
        {
            throw std::runtime_error(
                    "Error when creating mapped VTEC ionospheric correction: radiometric correction not "
                    "recognized." );
        }
    }
    else
    {
        throw std::runtime_error(
                "Error when creating mapped VTEC ionospheric correction: correction is only valid for "
                "radiometric types." );
    }
}

double MappedVtecIonosphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
{
    // Retrieve state and time of receiver and transmitter
    Eigen::Vector6d legTransmitterState, legReceiverState;
    double legTransmissionTime, legReceptionTime;
    getTransmissionReceptionTimesAndStates( linkEndsStates,
                                            linkEndsTimes,
                                            currentMultiLegTransmitterIndex,
                                            legTransmitterState,
                                            legReceiverState,
                                            legTransmissionTime,
                                            legReceptionTime );

    double groundStationTime;
    Eigen::Vector6d groundStationState, spacecraftState;
    double currentFrequency = TUDAT_NAN;

    if( isUplinkCorrection_ )
    {
        groundStationTime = legTransmissionTime;
        groundStationState = legTransmitterState;
        spacecraftState = legReceiverState;
        currentFrequency = ancillarySettings->getIntermediateDoubleData( transmitter_frequency_intermediate, true );
    }
    else
    {
        groundStationTime = legReceptionTime;
        groundStationState = legReceiverState;
        spacecraftState = legTransmitterState;
        currentFrequency = ancillarySettings->getIntermediateDoubleData( received_frequency_intermediate, true );
    }

    //    double firstLegTransmissionTime = linkEndsTimes.front( );

    // Retrieve frequency bands
    std::vector< FrequencyBands > frequencyBands = std::vector< FrequencyBands >( { x_band } );

    double elevation = elevationFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ), groundStationTime );
    double azimuth = azimuthFunction_( spacecraftState.segment( 0, 3 ) - groundStationState.segment( 0, 3 ), groundStationTime );
    Eigen::Vector3d groundStationGeodeticPosition = groundStationGeodeticPositionFunction_( groundStationTime );
    double geodeticLatitude = groundStationGeodeticPosition( 1 );
    double geodeticLongitude = groundStationGeodeticPosition( 2 );

    // Moyer (2000), eq. 10-43
    const double zenithAngle = std::asin( bodyWithAtmosphereMeanEquatorialRadius_ /
                                          ( bodyWithAtmosphereMeanEquatorialRadius_ + vtecCalculator_->getReferenceIonosphereHeight( ) ) *
                                          std::cos( elevation ) );

    Eigen::Vector3d subIonosphericPointGeodeticPosition;
    // Altitude
    subIonosphericPointGeodeticPosition( 0 ) = TUDAT_NAN;
    // Latitude: Moyer (2000), eq. 10-44
    subIonosphericPointGeodeticPosition( 1 ) =
            std::asin( std::sin( geodeticLatitude ) * std::sin( elevation + zenithAngle ) +
                       std::cos( geodeticLatitude ) * std::cos( elevation + zenithAngle ) * std::cos( azimuth ) );
    // Longitude: Moyer (2000), eq. 10-45
    subIonosphericPointGeodeticPosition( 2 ) = geodeticLongitude;
    if( std::abs( std::cos( subIonosphericPointGeodeticPosition( 1 ) ) ) > 1e-12 )
    {
        subIonosphericPointGeodeticPosition( 2 ) += std::asin( std::cos( elevation + zenithAngle ) * std::sin( azimuth ) /
                                                               std::cos( subIonosphericPointGeodeticPosition( 1 ) ) );
    }

    // Jakowski et al. (2011), eqs. 1 and 2; IERS conventions 2010, section 9.4
    // Mapping of VTEC to STEC: 1 / cos(zenithAngle)
    return ( sign_ * firstOrderDelayCoefficient_ *
             vtecCalculator_->calculateVtec( groundStationTime, subIonosphericPointGeodeticPosition ) / std::pow( currentFrequency, 2.0 ) /
             std::cos( zenithAngle ) ) /
            physical_constants::getSpeedOfLight< double >( );
}

}  // namespace observation_models

}  // namespace tudat