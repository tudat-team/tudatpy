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
 *          T. Morley and F. Budnik (2007) , EFFECTS ON SPACECRAFT RADIOMETRIC DATA AT SUPERIOR SOLAR CONJUNCTION, 20th
 *              International Symposium on Space Flight Dynamics
 */

#ifndef TUDAT_SOLARCORONACORRECTION_H
#define TUDAT_SOLARCORONACORRECTION_H

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

// Abstract class. Doesn't implement any correction model.
class SolarCoronaCorrection: public LightTimeCorrection
{
public:

    /*!
     * Constructor.
     * @param lightTimeCorrectionType Type of light time correction.
     * @param observableType Observable type associated with the correction.
     * @param sunStateFunction State of the Sun as a function of time.
     * @param transmittedFrequencyFunction Function calculating the frequency at the current link given a vector with
     *     the frequency bands in each link of the model and the transmission time.
     */
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

    /*!
     * Computes the minimum distance along the line of sight to the Sun. Used in some solar correction models.
     * @param transmitterPositionWrtSun Position of the transmitter with respect to the Sun.
     * @param receiverPositionWrtSun Position of the receiver with respect to the Sun.
     * @return Distance
     */
    double computeMinimumDistanceOfLineOfSight(
            Eigen::Vector3d transmitterPositionWrtSun,
            Eigen::Vector3d receiverPositionWrtSun );

    /*!
     * Gets the frequency at the current leg using the ancillary settings.
     * @param ancillarySettings Ancillary settings.
     * @param transmissionTime Time at which the signal was transmitted.
     * @return Frequency at current leg.
     */
    double getCurrentFrequency(
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings,
            const double transmissionTime );


    /*!
     * Computes the electron density at a certain position and time.
     *
     * @param positionWrtSun Position where to compute the electron density.
     * @param time Time at which to compute the electron density.
     */
    virtual double computeElectronDensity( const Eigen::Vector3d& positionWrtSun, const double time )
    {
        return TUDAT_NAN;
    }

    /*!
     * Computes the integral of the electron density along the line of sight, using Gaussian quadrature. Requires the
     * base class to have defined an implementation for computeElectronDensity.
     * @param transmitterPositionWrtSun Position of the transmitter with respect to the Sun.
     * @param receiverPositionWrtSun Position of the receiver with respect to the Sun.
     * @param time Time at which to compute the integral (i.e. assuming that the light time between the receiver and
     *      transmitter is zero). Giving the uncertainty of solar corona models and the associated time scales, it is
     *      unlikely that it ever becomes necessary to make a distinction between the reception and transmission time.
     * @return Integral.
     */
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

// Inverse power series solar correction model, based on Verma et al. (2013).
class InversePowerSeriesSolarCoronaCorrection: public SolarCoronaCorrection
{
public:

    /*!
     * Constructor. Model describes corrections for electron density of the form \sum c r^{-k}, mostly following
     * Verma et al. (2013). r corresponds to the dimensionless distance: distance / solar_radius. Model neglects the
     * effect of the latitude wrt Sun in the electron density.
     *
     * See e.g. Verma et al. (2013) and Morley and Budnik (2007) for examples of coefficients.
     * The default values are taken from Aksim and Pavlov (2022).
     *
     * @param observableType Observable type associated with the correction.
     * @param sunStateFunction State of the Sun as a function of time.
     * @param transmittedFrequencyFunction Function calculating the frequency at the current link given a vector with
     *     the frequency bands in each link of the model and the transmission time.
     * @param coefficients c coefficients of the series.
     * @param positiveExponents k exponents of the series.
     * @param criticalPlasmaDensityDelayCoefficient Coefficient that multiplies the integral of the electron density to
     *      get the value of the correction. Corresponds to
     *      1 / ( 2 * electron_mass * vacuum_permittivity * (2 pi)^2 ) / electron_charge^2 )
     * @param sunRadius Radius of the Sun.
     */
    InversePowerSeriesSolarCoronaCorrection(
            const ObservableType observableType,
            const std::function< Eigen::Vector6d ( double time ) > sunStateFunction,
            const std::function< double ( std::vector< FrequencyBands > frequencyBands, double time ) > transmittedFrequencyFunction,
            const std::vector< double >& coefficients =
                    { 1.31 * 5.97e6 * std::pow( physical_constants::ASTRONOMICAL_UNIT, 2.0 ) / std::pow( 696e6, 2 ) },
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
            // Check if all exponents are positive
            if ( exponent <= 0 )
            {
                throw std::runtime_error( "Error when creating inverse power series solar corona correction: negative exponent was"
                                          "provided (" + std::to_string( exponent ) + "). All provided exponents should be positive." );
            }

            // Check if all exponents are integers
            if ( std::fmod( exponent, 1.0 ) != 0 )
            {
                exponentsAreIntegers_ = false;
            }
        }
    }

    /*!
    * Function to compute the light-time correction, assuming an inverse power series model for the electron density
    * distribution, according to Verma et al. (2013).
    * @param linkEndsStates List of states at each link end during observation.
    * @param linkEndsTimes List of times at each link end during observation.
    * @param currentMultiLegTransmitterIndex Index in the linkEndsStates and linkEndsTimes of the transmitter in the current link.
    * @param ancillarySettings Observation ancillary simulation settings.
    * @return
    */
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

private:

    /*!
     * Computes the electron density, according to Verma et al. (2013).
     *
     * @param positionWrtSun Position where to compute the electron density.
     * @param time Time at which to compute the electron density.
     */
    double computeElectronDensity( const Eigen::Vector3d& positionWrtSun, const double time ) override;

    /*!
     * Computes the integral of a single r^{-k} electron density term over the line of sight. Does so analytically, which
     * is only valid if k is an integer >= 0.
     *
     * @param receiverPositionWrtSun Position of the receiver with respect to the Sun.
     * @param sunReceiverTransmitterAngle Angle between Receiver-Sun and Receiver-Transmitter vectors.
     * @param receiverSunTransmitterAngle Angle between Sun-Receiver and Sun-Transmitter vectors.
     * @param positiveExponent Value of the k exponent. Should be an integer.
     * @return Integral value.
     */
    double computeSingleTermIntegralAnalytically(
            const Eigen::Vector3d& receiverPositionWrtSun,
            const double sunReceiverTransmitterAngle,
            const double receiverSunTransmitterAngle,
            const unsigned int positiveExponent );

    /*!
     * Computes the integral of (cos(x))^k analytically, for an integer k >= 0.
     * @param lowerBound Lower bound of the integral.
     * @param upperBound Upper bound of the integral.
     * @param positiveExponent Exponent k.
     * @return Integral value.
     */
    double computeCosinePowerIntegral(
            const double lowerBound,
            const double upperBound,
            const unsigned int positiveExponent );

    // Vector containing the c coefficients of the electron density model (\sum c r^{-k}).
    const std::vector< double > coefficients_;

    // Vector containing the k exponents of the electron density model (\sum c r^{-k}).
    const std::vector< double > positiveExponents_;

    // Boolean indicating whether all exponents are integers. If so, correction can be calculated fully analytically
    bool exponentsAreIntegers_;

    // Solar corona correction coefficient.
    const double criticalPlasmaDensityDelayCoefficient_;

    // Radius of the Sun.
    const double sunRadius_;

    // Cache containing the integral of the already computed cosine-powers integrals.
    std::map< unsigned int, double > cosinePowersIntegralsCache_;

};


} // namespace observation_models

} // namespace tudat

#endif //TUDAT_SOLARCORONACORRECTION_H
