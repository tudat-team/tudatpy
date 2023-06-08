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
 *          A.K. Verma, A. Fienga, J. Laskar, K. Issautier, H. Manche, and M. Gastineau (2013), Electron density distribution
 *              and solar plasma correction of radio signals using MGS, MEX, and VEX spacecraft navigation data and its
 *              application to planetary ephemerides, Astronomy and Astrophysics, 550, A124
 *          D. Aksim and D. Pavlov (2022), Improving the Solar Wind Density Model Used in Processing of Spacecraft Ranging
 *              Observations, Monthly Notices of the Royal Astronomical Society, Volume 514, Issue 3
 */

#ifndef TUDAT_SOLARCORONACORRECTION_H
#define TUDAT_SOLARCORONACORRECTION_H

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

// Abstract class
class SolarCoronaCorrection: public LightTimeCorrection
{
public:

    SolarCoronaCorrection(
            const LightTimeCorrectionType lightTimeCorrectionType,
            const ObservableType observableType,
            const std::function< Eigen::Vector6d ( double time ) > sunStateFunction,
            const std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction ):
        LightTimeCorrection( lightTimeCorrectionType ),
        sunStateFunction_( sunStateFunction ),
        transmittedFrequencyFunction_( transmittedFrequencyFunction )
    {
        // Sign according to Moyer (2000), section 10.4.2
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
                throw std::runtime_error( "Error when creating solar corona correction: radiometric correction not "
                                          "recognized." );
            }
        }
        else
        {
            throw std::runtime_error( "Error when creating solar corona correction: correction is only valid for "
                                      "radiometric types." );
        }
    }

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

protected:

    double computeMinimumDistanceOfLineOfSight(
            Eigen::Vector3d transmitterPositionWrtSun,
            Eigen::Vector3d receiverPositionWrtSun );

    double getCurrentFrequency(
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings,
            const double currentTime );

    // Returns NAN for models that don't require defining explicitly the electron density model
    virtual double computeElectronDensity( const Eigen::Vector3d& positionWrtSun, const double time )
    {
        return TUDAT_NAN;
    }

    double computeElectronDensityIntegralNumerically(
            const Eigen::Vector3d& transmitterPositionWrtSun,
            const Eigen::Vector3d& receiverPositionWrtSun,
            const double time );

    // Sign of the correction (+1 or -1)
    int sign_;

    // Function returning the state of the Sun
    const std::function< Eigen::Vector6d ( double time ) > sunStateFunction_;

    // Frequency at the link as a function of the frequency bands per link, and of the current time
    std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction_;

private:

};

// Neglects influence of Sun's latitude
class InversePowerSeriesSolarCoronaCorrection: public SolarCoronaCorrection
{
public:

    InversePowerSeriesSolarCoronaCorrection(
            const ObservableType observableType,
            const std::function< Eigen::Vector6d ( double time ) > sunStateFunction,
            const std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction,
            const std::vector< double >& coefficients =
                    { 1.31 * 5.97e-6 * std::pow( physical_constants::ASTRONOMICAL_UNIT, 2.0 ) / std::pow( 696e6, 2 ) },
            const std::vector< double >& positiveExponents = { 2.0 },
            const double criticalPlasmaDensityDelayCoefficient = 40.3,
            const double sunRadius = 696e6 ):
        SolarCoronaCorrection( inverse_power_series_solar_corona, observableType, sunStateFunction, transmittedFrequencyFunction ),
        coefficients_( coefficients ),
        positiveExponents_( positiveExponents ),
        criticalPlasmaDensityDelayCoefficient_( criticalPlasmaDensityDelayCoefficient ),
        sunRadius_( sunRadius )
    {
        if ( coefficients.size( ) != positiveExponents.size( ) )
        {
            throw std::runtime_error( "Error when creating inverse power series solar corona correction: number of coefficients ("
                + std::to_string( coefficients.size( ) ) + ") and number of exponents (" + std::to_string( positiveExponents.size( ) ) +
                ") are incompatible." );
        }

        exponentsAreIntegers_ = true;
        for ( double exponent : positiveExponents )
        {
            if ( exponent <= 0 )
            {
                throw std::runtime_error( "Error when creating inverse power series solar corona correction: negative exponent was"
                                          "provided (" + std::to_string( exponent ) + "). All provided exponents should be positive." );
            }

            if ( std::fmod( exponent, 1.0 ) != 0 )
            {
                exponentsAreIntegers_ = false;
            }
        }
    }

    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

private:

    double computeElectronDensity( const Eigen::Vector3d& positionWrtSun, const double time ) override;

    double computeSingleTermIntegralAnalytically(
            const Eigen::Vector3d& receiverPositionWrtSun,
            const double sunReceiverTransmitterAngle,
            const double receiverSunTransmitterAngle,
            const unsigned int positiveExponent );

    double computeCosinePowerIntegral(
            const double lowerBound,
            const double upperBound,
            const unsigned int positiveExponent );

    const std::vector< double > coefficients_;

    const std::vector< double > positiveExponents_;

    bool exponentsAreIntegers_;

    const double criticalPlasmaDensityDelayCoefficient_;

    const double sunRadius_;

    std::map< unsigned int, double > cosinePowersIntegrals_;

};


} // namespace observation_models

} // namespace tudat

#endif //TUDAT_SOLARCORONACORRECTION_H
