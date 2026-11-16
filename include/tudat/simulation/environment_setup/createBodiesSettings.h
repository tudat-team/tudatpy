/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODIESSETTINGS_H
#define TUDAT_CREATEBODIESSETTINGS_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/simulation/environment_setup/createBodyShapeModel.h"
#include "tudat/simulation/environment_setup/createClimateModel.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createCameras.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"

namespace tudat
{

namespace simulation_setup
{

class SpaceTimeMetricSettings;

//! Struct holding settings for a body to be created.
/*!
 *  Struct holding settings for a body to be created. From the settings, a Body object is
 *  created by the createBodies function. Default values can be generated from the function in
 *  defaultBodies.h.
 */
struct BodySettings {
    //! Constant mass.
    double constantMass = TUDAT_NAN;

    //! Settings for the atmosphere model that the body is to contain.
    std::shared_ptr< AtmosphereSettings > atmosphereSettings;

    //! Settings for the ephemeris model that the body is to contain.
    std::shared_ptr< EphemerisSettings > ephemerisSettings;

    //! Settings for the gravity field model that the body is to contain.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    //! Settings for the rotation model that the body is to contain.
    std::shared_ptr< RotationModelSettings > rotationModelSettings;

    //! Settings for the shape model that the body is to contain.
    std::shared_ptr< BodyShapeSettings > shapeModelSettings;

    //! Settings for the radiation source model that the body is to contain.
    std::shared_ptr< RadiationSourceModelSettings > radiationSourceModelSettings;

    //! Settings for the radiation pressure target model that the body is to contain.
    std::shared_ptr< RadiationPressureTargetModelSettings > radiationPressureTargetModelSettings;

    //! Settings for the aerodynamic coefficients that the body is to contain.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    std::shared_ptr< RigidBodyPropertiesSettings > rigidBodyPropertiesSettings;

    std::shared_ptr< FullPanelledBodySettings > bodyExteriorPanelSettings_;

    //! Settings for variations of the gravity field of the body.
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

    std::vector< std::shared_ptr< BodyDeformationSettings > > bodyDeformationSettings;

    std::vector< std::shared_ptr< GroundStationSettings > > groundStationSettings;

    std::shared_ptr< ClimateModelSettings > climateModelSettings;

    std::vector< std::shared_ptr< CameraSettings > > cameraSettings;
};

class SpaceTimePropertiesSettings
{
public:
    SpaceTimePropertiesSettings( const std::shared_ptr< SpaceTimeMetricSettings >& metricSettings = nullptr,
                                 const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSet =
                                         std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 ),
                                 const double equivalencePrincipleLpiViolationParameter = 0.0 ):
        metricSettings_( metricSettings ), ppnParameterSet_( ppnParameterSet ),
        equivalencePrincipleLpiViolationParameter_( equivalencePrincipleLpiViolationParameter )
    {
        if( ppnParameterSet_ == nullptr )
        {
            ppnParameterSet_ = std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 );
        }
    }

    std::shared_ptr< SpaceTimeMetricSettings > getMetricSettings( ) const
    {
        return metricSettings_;
    }

    void setMetricSettings( const std::shared_ptr< SpaceTimeMetricSettings >& metricSettings )
    {
        metricSettings_ = metricSettings;
    }

    std::shared_ptr< relativity::PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    void setPpnParameterSet( const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSet )
    {
        if( ppnParameterSet == nullptr )
        {
            throw std::runtime_error( "Error when setting PPN parameter set, input is nullptr." );
        }
        ppnParameterSet_ = ppnParameterSet;
    }

    double getEquivalencePrincipleLpiViolationParameter( ) const
    {
        return equivalencePrincipleLpiViolationParameter_;
    }

    void setEquivalencePrincipleLpiViolationParameter( const double equivalencePrincipleLpiViolationParameter )
    {
        equivalencePrincipleLpiViolationParameter_ = equivalencePrincipleLpiViolationParameter;
    }

private:
    std::shared_ptr< SpaceTimeMetricSettings > metricSettings_;

    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

    double equivalencePrincipleLpiViolationParameter_;
};

//! Factory function for PPN parameter settings.
/*!
 *  \param parameterGamma First-order PPN parameter gamma.
 *  \param parameterBeta First-order PPN parameter beta.
 *  \param parameterEpsilon Second-order post-Newtonian parameter epsilon.
 *  \return Shared pointer to PPN parameter settings.
 */
inline std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet( const double parameterGamma = 1.0,
                                                                       const double parameterBeta = 1.0,
                                                                       const double parameterEpsilon = 0.0 )
{
    return std::make_shared< relativity::PPNParameterSet >( parameterGamma, parameterBeta, 0.0, parameterEpsilon );
}

//! Factory function for system-level space-time properties settings.
/*!
 *  \param metricSettings Optional settings defining which metric model to build.
 *  \param ppnParameterSet Optional PPN parameter settings. If null, GR defaults are used.
 *  \param equivalencePrincipleLpiViolationParameter Equivalence-principle local-position-invariance violation parameter.
 *  \return Shared pointer to space-time properties settings.
 */
inline std::shared_ptr< SpaceTimePropertiesSettings > spaceTimePropertiesSettings(
        const std::shared_ptr< SpaceTimeMetricSettings >& metricSettings = nullptr,
        const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSet = nullptr,
        const double equivalencePrincipleLpiViolationParameter = 0.0 )
{
    return std::make_shared< SpaceTimePropertiesSettings >( metricSettings, ppnParameterSet, equivalencePrincipleLpiViolationParameter );
}

class BodyListSettings
{
public:
    BodyListSettings( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" ):
        bodySettings_( std::map< std::string, std::shared_ptr< BodySettings > >( ) ), frameOrigin_( frameOrigin ),
        frameOrientation_( frameOrientation ), spaceTimePropertiesSettings_( std::make_shared< SpaceTimePropertiesSettings >( ) )
    {}

    BodyListSettings( const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings,
                      const std::string frameOrigin = "SSB",
                      const std::string frameOrientation = "ECLIPJ2000" ):
        bodySettings_( bodySettings ), frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ),
        spaceTimePropertiesSettings_( std::make_shared< SpaceTimePropertiesSettings >( ) )
    {}

    std::shared_ptr< BodySettings > at( const std::string& bodyName ) const
    {
        if( bodySettings_.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving body settings for " + bodyName + ", no such body exist in this object." );
        }

        return bodySettings_.at( bodyName );
    }

    std::shared_ptr< BodySettings > get( const std::string& bodyName ) const
    {
        return at( bodyName );
    }

    int count( const std::string& bodyName ) const
    {
        return bodySettings_.count( bodyName );
    }

    void addSettings( std::shared_ptr< BodySettings > settingsToAdd, const std::string bodyName )
    {
        bodySettings_[ bodyName ] = settingsToAdd;
    }

    void addSettings( const std::string bodyName )
    {
        bodySettings_[ bodyName ] = std::make_shared< BodySettings >( );
    }

    void clear( )
    {
        bodySettings_.clear( );
    }

    void resetFrames( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" )
    {
        frameOrigin_ = frameOrigin;
        frameOrientation_ = frameOrientation;
    }

    std::string getFrameOrigin( ) const
    {
        return frameOrigin_;
    }

    std::string getFrameOrientation( ) const
    {
        return frameOrientation_;
    }

    std::map< std::string, std::shared_ptr< BodySettings > > getMap( ) const
    {
        return bodySettings_;
    }

    std::shared_ptr< SpaceTimePropertiesSettings > getSpaceTimeSettings( ) const
    {
        return spaceTimePropertiesSettings_;
    }

    void setSpaceTimeSettings( const std::shared_ptr< SpaceTimePropertiesSettings >& spaceTimePropertiesSettings )
    {
        if( spaceTimePropertiesSettings == nullptr )
        {
            throw std::runtime_error( "Error when setting space-time settings, input is nullptr." );
        }
        spaceTimePropertiesSettings_ = spaceTimePropertiesSettings;
    }

private:
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings_;

    std::string frameOrigin_;

    std::string frameOrientation_;

    std::shared_ptr< SpaceTimePropertiesSettings > spaceTimePropertiesSettings_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEBODIESSETTINGS_H
