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

// Class defining tabulated tropospheric corrections
class TabulatedTroposphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    TabulatedTroposphericCorrectionSettings(
            const AtmosphericCorrectionPerStationAndSpacecraftType& troposphericDryCorrection,
            const AtmosphericCorrectionPerStationAndSpacecraftType& troposphericWetCorrection,
            const std::string& bodyWithAtmosphere = "Earth",
            const TroposphericMappingModel troposphericMappingModel = niell ):
        LightTimeCorrectionSettings( tabulated_tropospheric ),
        troposphericDryCorrection_( troposphericDryCorrection ),
        troposphericWetCorrection_( troposphericWetCorrection ),
        bodyWithAtmosphere_( bodyWithAtmosphere ),
        troposphericMappingModelType_( troposphericMappingModel )
    { }

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

    AtmosphericCorrectionPerStationAndSpacecraftType troposphericDryCorrection_;

    AtmosphericCorrectionPerStationAndSpacecraftType troposphericWetCorrection_;

    std::string bodyWithAtmosphere_;

    TroposphericMappingModel troposphericMappingModelType_;

};

// Class defining tabulated ionospheric corrections
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
        const LinkEndId& transmitter,
        const LinkEndId& receiver,
        const ObservableType observableType = undefined_observation_model);

std::shared_ptr< TroposhericElevationMapping > createTroposphericElevationMapping(
        const TroposphericMappingModel troposphericMappingModelType,
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEndId& transmitter,
        const LinkEndId& receiver,
        const bool isUplinkCorrection );

} // namespace observation_models

} // namespace tudat


#endif // TUDAT_CREATELIGHTTIMECORRECTION_H
