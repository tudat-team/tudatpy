/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/neQuick2IonosphericCorrection.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/quadrature/gaussianQuadrature.h"

#include <cmath>

namespace tudat
{
namespace observation_models
{

NeQuick2IonosphericCorrection::NeQuick2IonosphericCorrection(
    std::shared_ptr< environment::NeQuick2Model > neQuick2Model,
    std::shared_ptr< environment::IonexConstrainedNeQuick2Model > rescaledModel,
    const ObservableType observableType,
    std::function< Eigen::Vector6d( double ) > earthStateFunction,
    std::function< Eigen::Matrix3d( double ) > earthRotationToBodyFixedFunction,
    double earthEquatorialRadius,
    double firstOrderDelayCoefficient,
    int quadratureOrder )
    : LightTimeCorrection( nequick2_ionospheric ),
      neQuick2Model_( neQuick2Model ),
      rescaledModel_( rescaledModel ),
      earthStateFunction_( earthStateFunction ),
      earthRotationToBodyFixedFunction_( earthRotationToBodyFixedFunction ),
      earthEquatorialRadius_( earthEquatorialRadius ),
      firstOrderDelayCoefficient_( firstOrderDelayCoefficient ),
      quadratureOrder_( quadratureOrder )
{
    if( isRadiometricObservableType( observableType ) )
    {
        if( isGroupVelocityBasedObservableType( observableType ) )
        {
            sign_ = 1;
        }
        else if( isPhaseVelocityBasedObservableType( observableType ) )
        {
            sign_ = -1;
        }
        else
        {
            throw std::runtime_error(
                "Error when creating NeQuick-2 ionospheric correction: radiometric type not recognized." );
        }
    }
    else
    {
        throw std::runtime_error(
            "Error when creating NeQuick-2 ionospheric correction: only valid for radiometric observable types." );
    }
}

double NeQuick2IonosphericCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
    const std::vector< Eigen::Vector6d >& linkEndsStates,
    const std::vector< double >& linkEndsTimes,
    const unsigned int currentMultiLegTransmitterIndex,
    const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings )
{
    // Extract transmitter/receiver states and times
    Eigen::Vector6d legTransmitterState, legReceiverState;
    double legTransmissionTime, legReceptionTime;
    getTransmissionReceptionTimesAndStates( linkEndsStates,
                                            linkEndsTimes,
                                            currentMultiLegTransmitterIndex,
                                            legTransmitterState,
                                            legReceiverState,
                                            legTransmissionTime,
                                            legReceptionTime );

    // Get frequency from ancillary settings (same pattern as solar corona correction)
    double currentFrequency = TUDAT_NAN;
    if( currentMultiLegTransmitterIndex == 0 )
    {
        currentFrequency = ancillarySettings->getIntermediateDoubleData( transmitter_frequency_intermediate, true );
    }
    else if( currentMultiLegTransmitterIndex == 2 )
    {
        currentFrequency = ancillarySettings->getIntermediateDoubleData( received_frequency_intermediate, true );
    }
    else
    {
        throw std::runtime_error(
            "Error when computing NeQuick-2 ionospheric correction: unsupported link index " +
            std::to_string( currentMultiLegTransmitterIndex ) );
    }

    // Use the average time for the correction
    double correctionTime = 0.5 * ( legTransmissionTime + legReceptionTime );

    // Transform positions to Earth body-fixed frame
    Eigen::Vector6d earthState = earthStateFunction_( correctionTime );
    Eigen::Matrix3d rotationToBodyFixed = earthRotationToBodyFixedFunction_( correctionTime );

    Eigen::Vector3d txInertial = legTransmitterState.head< 3 >( ) - earthState.head< 3 >( );
    Eigen::Vector3d rxInertial = legReceiverState.head< 3 >( ) - earthState.head< 3 >( );

    Eigen::Vector3d txBodyFixed = rotationToBodyFixed * txInertial;
    Eigen::Vector3d rxBodyFixed = rotationToBodyFixed * rxInertial;

    // Check if both endpoints are below the ionosphere (~80 km)
    double txAltitude = ( txBodyFixed.norm( ) - earthEquatorialRadius_ ) / 1.0e3;
    double rxAltitude = ( rxBodyFixed.norm( ) - earthEquatorialRadius_ ) / 1.0e3;
    if( txAltitude < 80.0 && rxAltitude < 80.0 )
    {
        return 0.0;
    }

    // Compute slant TEC
    double stec;
    if( rescaledModel_ != nullptr )
    {
        // IONEX-constrained NeQuick-2
        stec = rescaledModel_->computeRescaledSlantTec(
            txBodyFixed, rxBodyFixed, correctionTime, earthEquatorialRadius_, quadratureOrder_ );
    }
    else
    {
        // Free-running NeQuick-2 (no IONEX constraint)
        int month;
        double ut;
        environment::NeQuick2Model::timeToMonthAndUT( correctionTime, month, ut );
        double flx = neQuick2Model_->getSolarFluxFunction( )( correctionTime );
        flx = std::max( 63.0, std::min( flx, 193.0 ) );

        Eigen::Vector3d rayDir = rxBodyFixed - txBodyFixed;
        double pathLength = rayDir.norm( );

        std::function< double( double ) > electronDensityAlongRay =
            [ this, &txBodyFixed, &rayDir, month, flx, ut ]( double t ) {
                Eigen::Vector3d pos = txBodyFixed + t * rayDir;
                double radius = pos.norm( );
                double heightKm = ( radius - earthEquatorialRadius_ ) / 1.0e3;
                if( heightKm < 0.0 ) return 0.0;

                double latRad = std::asin( pos.z( ) / radius );
                double lonRad = std::atan2( pos.y( ), pos.x( ) );
                double latDeg = latRad * 180.0 / mathematical_constants::PI;
                double lonDeg = lonRad * 180.0 / mathematical_constants::PI;

                return neQuick2Model_->computeElectronDensity(
                    heightKm, latDeg, lonDeg, month, flx, ut );
            };

        numerical_quadrature::GaussianQuadrature< double, double > quadrature(
            electronDensityAlongRay, 0.0, 1.0, quadratureOrder_ );

        stec = quadrature.getQuadrature( ) * pathLength;
    }

    // Convert STEC to light-time correction: delay = sign * 40.3 * STEC / f^2 / c
    return sign_ * firstOrderDelayCoefficient_ * stec /
           std::pow( currentFrequency, 2.0 ) /
           physical_constants::getSpeedOfLight< double >( );
}

}  // namespace observation_models
}  // namespace tudat
