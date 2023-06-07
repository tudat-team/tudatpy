/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATELIGHTTIMECORRECTION_H
#define TUDAT_CREATELIGHTTIMECORRECTION_H

#include <Eigen/Core>

#include <memory>
#include <functional>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/astro/observation_models/corrections/tabulatedMediaCorrection.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/solarActivityData.h"

namespace tudat
{

namespace observation_models
{

//! Typedef for function calculating light-time correction in light-time calculation loop.
typedef std::function< double(
        const Eigen::Vector6d&, const Eigen::Vector6d&,
        const double, const double ) > LightTimeCorrectionFunction;

//! Base class for light-time correction settings.
/*!
 *  Base class for light-time correction settings. This class is not used for calculations of corrections,
 *  but is used for the purpose of defining the light time correction properties.
 *  The createLightTimeCorrections function produces the classes that calculate the actual corrections, based on settings
 *  in instances of derived LightTimeCorrectionSettings classes.
 */
class LightTimeCorrectionSettings
{
public:

    //! Constructor, takes light-time correction type.
    /*!
     *  \param correctionType Type of light-time correction that is to be created
     */
    LightTimeCorrectionSettings( const LightTimeCorrectionType correctionType ):
        correctionType_( correctionType ){ }

    //! Default destructor.
    virtual ~LightTimeCorrectionSettings( ){ }

    //! Function returning the type of light-time correction that is to be created
    /*!
     *  Function returning the type of light-time correction that is to be created
     *  \return Type of light-time correction that is to be created
     */
    LightTimeCorrectionType getCorrectionType( )
    {
        return correctionType_;
    }

protected:

    //! Type of light-time correction that is to be created
    LightTimeCorrectionType correctionType_;
};

//! Typedef for a list of list time correction settings per link end
typedef std::map< LinkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > > > LightTimeCorrectionSettingsMap;

//! Class to defining settings for first-order relativistic light time correction (Shapiro time delay)  due to a
//! set of point masses
class FirstOrderRelativisticLightTimeCorrectionSettings: public LightTimeCorrectionSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param perturbingBodies List of bodies for which the point masses are used to compute the light-time correction.
     */
    FirstOrderRelativisticLightTimeCorrectionSettings( const std::vector< std::string >& perturbingBodies ):
        LightTimeCorrectionSettings( first_order_relativistic ), perturbingBodies_( perturbingBodies ){ }

    //! Destructor
    ~FirstOrderRelativisticLightTimeCorrectionSettings( ){ }

    //! Function returning the list of bodies for which the point masses are used to compute the light-time correction.
    /*!
     *  Function returning the list of bodies for which the point masses are used to compute the light-time correction.
     *  \return List of bodies for which the point masses are used to compute the light-time correction.
     */
    std::vector< std::string > getPerturbingBodies( ){ return perturbingBodies_; }

private:

    //! List of bodies for which the point masses are used to compute the light-time correction.
    std::vector< std::string > perturbingBodies_;

};

// Class defining  settings for tabulated tropospheric corrections
class TabulatedTroposphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    TabulatedTroposphericCorrectionSettings(
            const AtmosphericCorrectionPerStationAndSpacecraftType& troposphericDryCorrectionAdjustment,
            const AtmosphericCorrectionPerStationAndSpacecraftType& troposphericWetCorrectionAdjustment,
            const std::string& bodyWithAtmosphere = "Earth",
            const TroposphericMappingModel troposphericMappingModel = niell,
            const AtmosphericCorrectionPerStationAndSpacecraftType& troposphericDryCorrection =
                    input_output::createDefaultTroposphericDryCorrection( ),
            const AtmosphericCorrectionPerStationAndSpacecraftType& troposphericWetCorrection =
                    input_output::createDefaultTroposphericWetCorrection( ) ):
        LightTimeCorrectionSettings( tabulated_tropospheric ),
        troposphericDryCorrectionAdjustment_( troposphericDryCorrectionAdjustment ),
        troposphericWetCorrectionAdjustment_( troposphericWetCorrectionAdjustment ),
        troposphericDryCorrection_( troposphericDryCorrection ),
        troposphericWetCorrection_( troposphericWetCorrection ),
        bodyWithAtmosphere_( bodyWithAtmosphere ),
        troposphericMappingModelType_( troposphericMappingModel )
    { }

    AtmosphericCorrectionPerStationAndSpacecraftType getTroposphericDryCorrectionAdjustment( )
    {
        return troposphericDryCorrectionAdjustment_;
    }

    AtmosphericCorrectionPerStationAndSpacecraftType getTroposphericWetCorrectionAdjustment( )
    {
        return troposphericWetCorrectionAdjustment_;
    }

    AtmosphericCorrectionPerStationAndSpacecraftType getTroposphericDryCorrection( )
    {
        return troposphericDryCorrection_;
    }

    AtmosphericCorrectionPerStationAndSpacecraftType getTroposphericWetCorrection( )
    {
        return troposphericWetCorrection_;
    }

    TroposphericMappingModel getTroposphericMappingModelType( )
    {
        return troposphericMappingModelType_;
    }

    std::string getBodyWithAtmosphere( )
    {
        return bodyWithAtmosphere_;
    }

private:

    AtmosphericCorrectionPerStationAndSpacecraftType troposphericDryCorrectionAdjustment_;

    AtmosphericCorrectionPerStationAndSpacecraftType troposphericWetCorrectionAdjustment_;

    AtmosphericCorrectionPerStationAndSpacecraftType troposphericDryCorrection_;

    AtmosphericCorrectionPerStationAndSpacecraftType troposphericWetCorrection_;

    std::string bodyWithAtmosphere_;

    TroposphericMappingModel troposphericMappingModelType_;

};

class SaastamoinenTroposphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    SaastamoinenTroposphericCorrectionSettings(
            const std::string& bodyWithAtmosphere = "Earth",
            const TroposphericMappingModel troposphericMappingModel = niell,
            const WaterVaporPartialPressureModel waterVaporPartialPressureModel = tabulated ):
        LightTimeCorrectionSettings( saastamoinen_tropospheric ),
        bodyWithAtmosphere_( bodyWithAtmosphere ),
        troposphericMappingModelType_( troposphericMappingModel ),
        waterVaporPartialPressureModelType_( waterVaporPartialPressureModel )
    { }

    std::string getBodyWithAtmosphere( )
    {
        return bodyWithAtmosphere_;
    }

    TroposphericMappingModel getTroposphericMappingModelType( )
    {
        return troposphericMappingModelType_;
    }

    WaterVaporPartialPressureModel getWaterVaporPartialPressureModelType( )
    {
        return waterVaporPartialPressureModelType_;
    }

private:

    std::string bodyWithAtmosphere_;

    TroposphericMappingModel troposphericMappingModelType_;

    WaterVaporPartialPressureModel waterVaporPartialPressureModelType_;

};

// Class defining settings for tabulated ionospheric corrections
class TabulatedIonosphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    TabulatedIonosphericCorrectionSettings(
            const AtmosphericCorrectionPerStationAndSpacecraftType& referenceRangeCorrection,
            const double referenceFrequency = 2295e6,
            const std::string& bodyWithAtmosphere = "Earth" ):
        LightTimeCorrectionSettings( tabulated_ionospheric ),
        referenceRangeCorrection_( referenceRangeCorrection ),
        referenceFrequency_( referenceFrequency ),
        bodyWithAtmosphere_( bodyWithAtmosphere )
    { }

    AtmosphericCorrectionPerStationAndSpacecraftType getReferenceRangeCorrection( )
    {
        return referenceRangeCorrection_;
    }

    double getReferenceFrequency( )
    {
        return referenceFrequency_;
    }

    std::string getBodyWithAtmosphere( )
    {
        return bodyWithAtmosphere_;
    }

private:

    AtmosphericCorrectionPerStationAndSpacecraftType referenceRangeCorrection_;

    double referenceFrequency_;

    std::string bodyWithAtmosphere_;

};

// Class defining settings for Jakowski ionospheric corrections
class JakowskiIonosphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    JakowskiIonosphericCorrectionSettings(
            const double ionosphereHeight = 400.0e3,
            const double firstOrderDelayCoefficient = 40.3,
            const input_output::solar_activity::SolarActivityDataMap& solarActivityData =
                    input_output::solar_activity::readSolarActivityData(
                            paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt" ),
            const double geomagneticPoleLatitude = unit_conversions::convertDegreesToRadians( 80.9 ),
            const double geomagneticPoleLongitude = unit_conversions::convertDegreesToRadians( -72.6 ),
            const bool useUtcTimeForLocalTimeComputation = false,
            const std::string& bodyWithAtmosphere = "Earth"):
        LightTimeCorrectionSettings( jakowski_vtec_ionospheric ),
        ionosphereHeight_( ionosphereHeight ),
        firstOrderDelayCoefficient_( firstOrderDelayCoefficient ),
        solarActivityData_( solarActivityData ),
        geomagneticPoleLatitude_( geomagneticPoleLatitude ),
        geomagneticPoleLongitude_( geomagneticPoleLongitude ),
        useUtcTimeForLocalTime_( useUtcTimeForLocalTimeComputation ),
        bodyWithAtmosphere_( bodyWithAtmosphere )
    { }

    double getIonosphereHeight( )
    {
        return ionosphereHeight_;
    }

    double getFirstOrderDelayCoefficient( )
    {
        return firstOrderDelayCoefficient_;
    }

    input_output::solar_activity::SolarActivityDataMap getSolarActivityData( )
    {
        return solarActivityData_;
    }

    double getGeomagneticPoleLatitude( )
    {
        return geomagneticPoleLatitude_;
    }

    double getGeomagneticPoleLongitude( )
    {
        return geomagneticPoleLongitude_;
    }

    bool getUseUtcTimeForLocalTime( )
    {
        return useUtcTimeForLocalTime_;
    }

    std::string getBodyWithAtmosphere( )
    {
        return bodyWithAtmosphere_;
    }

private:

    const double ionosphereHeight_;

    const double firstOrderDelayCoefficient_;

    input_output::solar_activity::SolarActivityDataMap solarActivityData_;

    const double geomagneticPoleLatitude_;

    const double geomagneticPoleLongitude_;

    const bool useUtcTimeForLocalTime_;

    const std::string bodyWithAtmosphere_;

};

// Class defining settings for tabulated ionospheric corrections
class InversePowerSeriesSolarCoronaCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    InversePowerSeriesSolarCoronaCorrectionSettings(
            const std::vector< double >& coefficients = { 1.31 * 5.97e-6 },
            const std::vector< double >& positiveExponents = { 2.0 },
            const double criticalPlasmaDensityDelayCoefficient = 40.3,
            const std::string& sunBodyName = "Sun" ):
        LightTimeCorrectionSettings( inverse_power_series_solar_corona ),
        coefficients_( coefficients ),
        positiveExponents_( positiveExponents ),
        criticalPlasmaDensityDelayCoefficient_( criticalPlasmaDensityDelayCoefficient ),
        sunBodyName_( sunBodyName )
    { }

    std::vector< double > getCoefficients( )
    {
        return coefficients_;
    }

    std::vector< double > getPositiveExponents( )
    {
        return positiveExponents_;
    }

    double getCriticalPlasmaDensityDelayCoefficient( )
    {
        return criticalPlasmaDensityDelayCoefficient_;
    }

    std::string getSunBodyName( )
    {
        return sunBodyName_;
    }


private:

    const std::vector< double > coefficients_;

    const std::vector< double > positiveExponents_;

    const double criticalPlasmaDensityDelayCoefficient_;

    const std::string& sunBodyName_;

};

inline std::shared_ptr< LightTimeCorrectionSettings > firstOrderRelativisticLightTimeCorrectionSettings(
        const std::vector< std::string >& perturbingBodies )
{
    return std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies );
}

//! Function to create object that computes a single (type of) correction to the light-time
/*!
 * Function to create object that computes a single (type of) correction to the light-time. Returns nullptr if the
 * selected corrections aren't applicable to the specified pair of transmitter and receiver.
 * \param correctionSettings User-defined settings for the light-time correction that is to be created
 * \param bodies List of body objects that constitutes the environment
 * \param transmitter Id of transmitting body/reference point (first/second)
 * \param receiver Id of receiving body/reference point (first/second)
 * \return Object for computing required light-time correction
 */
std::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const std::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType,
        const ObservableType observableType = undefined_observation_model);

std::shared_ptr< TroposhericElevationMapping > createTroposphericElevationMapping(
        const TroposphericMappingModel troposphericMappingModelType,
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEndId& transmitter,
        const LinkEndId& receiver,
        const bool isUplinkCorrection );

std::function< double ( std::vector< FrequencyBands >, double ) > createLinkFrequencyFunction(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType );

} // namespace observation_models

} // namespace tudat


#endif // TUDAT_CREATELIGHTTIMECORRECTION_H
