/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODIES_H
#define TUDAT_CREATEBODIES_H

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/createBodyShapeModel.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createRadiationPressureInterface.h"
#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"

namespace tudat
{

namespace simulation_setup
{

//! Struct holding settings for a body to be created.
/*!
 *  Struct holding settings for a body to be created. From the settings, a Body object is
 *  created by the createBodies function. Default values can be generated from the function in
 *  defaultBodies.h.
 */
struct BodySettings
{
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

    //! Settings for the radiations pressure interfaces that the body is to contain (source body as key).
    // RP-OLD
    std::map< std::string,
    std::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;

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


};

void addAerodynamicCoefficientInterface(
        const SystemOfBodies& bodies, const std::string bodyName,
        const std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings );

// RP-OLD
void addRadiationPressureInterface(
        const SystemOfBodies& bodies, const std::string bodyName,
        const std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings );

void addRadiationPressureTargetModel(
    const SystemOfBodies& bodies, const std::string bodyName,
    const std::shared_ptr< RadiationPressureTargetModelSettings > radiationPressureSettings );


void addRotationModel(
        const SystemOfBodies& bodies, const std::string bodyName,
        const std::shared_ptr< RotationModelSettings > rotationModelSettings );

void addGravityFieldModel(
    const SystemOfBodies& bodies, const std::string bodyName,
    const std::shared_ptr< GravityFieldSettings > gravityFieldSettings,
    const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings =
        std::vector< std::shared_ptr< GravityFieldVariationSettings > >( ));

void addRigidBodyProperties(
    const SystemOfBodies& bodies, const std::string bodyName,
    const std::shared_ptr< RigidBodyPropertiesSettings > rigidBodyProperties );



class BodyListSettings
{
public:

    BodyListSettings( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" ):
        bodySettings_( std::map< std::string, std::shared_ptr< BodySettings > >( ) ),
        frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ){ }

    BodyListSettings( const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings,
                      const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" ):
        bodySettings_( bodySettings ), frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ){ }

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

    std::string getFrameOrigin( ) const { return frameOrigin_; }

    std::string getFrameOrientation( ) const { return frameOrientation_; }

    std::map< std::string, std::shared_ptr< BodySettings > > getMap( ) const { return bodySettings_; }


private:

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings_;

    std::string frameOrigin_;

    std::string frameOrientation_;
};

void setSimpleRotationSettingsFromSpice(
        const BodyListSettings& bodySettings, const std::string& bodyName, const double  spiceEvaluation );


//! Function that determines the order in which bodies are to be created
/*!
 * Function that determines the order in which bodies are to be created, to ensure that any dependency between body models are
 * correctly handled. Currently, no dependencies exist that force any particular creation order.
 * \param bodySettings Map of body settings (with map key the body name)
 * \return List of pairs: name and body settings of that body
 */
std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings );

//! Function to create a map of bodies objects.
/*!
 *  Function to create a map of body objects based on model-specific settings for the bodies,
 *  containing settings for each relevant environment model.
 *  \param bodySettings List of settings for the bodies that are to be created, defined as a map of
 *  pointers to an object of class BodySettings
 *  \return List of bodies created according to settings in bodySettings.
 */
template< typename StateScalarType = double , typename TimeType = double >
SystemOfBodies createSystemOfBodies(
        const BodyListSettings& bodySettings )
{
    std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > orderedBodySettings
            = determineBodyCreationOrder( bodySettings.getMap( ) );

    // Declare map of bodies that is to be returned.
    SystemOfBodies bodyList = SystemOfBodies(
                bodySettings.getFrameOrigin( ), bodySettings.getFrameOrientation( ) );

    // Create empty body objects.
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        bodyList.createEmptyBody< StateScalarType, TimeType >( orderedBodySettings.at( i ).first, false );
    }

    // Define constant mass for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        const double constantMass = orderedBodySettings.at( i ).second->constantMass;
        if ( constantMass == constantMass )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setConstantBodyMass( constantMass );
        }
    }

    // Create ephemeris objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->ephemerisSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setEphemeris(
                        createBodyEphemeris< StateScalarType, TimeType >(
                            orderedBodySettings.at( i ).second->ephemerisSettings,
                                             orderedBodySettings.at( i ).first ) );
        }
    }

    // Create atmosphere model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->atmosphereSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setAtmosphereModel(
                        createAtmosphereModel( orderedBodySettings.at( i ).second->atmosphereSettings,
                                               orderedBodySettings.at( i ).first ) );
        }
    }

    // Create body shape model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->shapeModelSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setShapeModel(
                        createBodyShapeModel( orderedBodySettings.at( i ).second->shapeModelSettings,
                                              orderedBodySettings.at( i ).first ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->rotationModelSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRotationalEphemeris(
                        createRotationModel( orderedBodySettings.at( i ).second->rotationModelSettings,
                                             orderedBodySettings.at( i ).first, bodyList ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->bodyExteriorPanelSettings_ != nullptr )
        {
            addBodyExteriorPanelledShape( orderedBodySettings.at( i ).second->bodyExteriorPanelSettings_ ,
                                          orderedBodySettings.at( i ).first, bodyList );
        }
    }
    std::vector< std::shared_ptr< BodyPanelSettings > > bodyExteriorPanelSettings_;


    // Create gravity field model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->gravityFieldSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setGravityFieldModel(
                        createGravityFieldModel( orderedBodySettings.at( i ).second->gravityFieldSettings,
                                                 orderedBodySettings.at( i ).first, bodyList,
                                                 orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }


    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->gravityFieldVariationSettings.size( ) > 0 )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setGravityFieldVariationSet(
                        createGravityFieldModelVariationsSet(
                            orderedBodySettings.at( i ).first, bodyList,
                            orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    // Create gravity field model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->rigidBodyPropertiesSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setMassProperties(
                createRigidBodyProperties( orderedBodySettings.at( i ).second->rigidBodyPropertiesSettings,
                                          orderedBodySettings.at( i ).first, bodyList ) );
        }
    }


    // Create aerodynamic coefficient interface objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface(
                            orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings,
                            orderedBodySettings.at( i ).first, bodyList ) );
        }
    }


    // Create radiation pressure coefficient objects for each body (if required).
    // RP-OLD
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        std::map< std::string, std::shared_ptr< RadiationPressureInterfaceSettings > >
                radiationPressureSettings
                = orderedBodySettings.at( i ).second->radiationPressureSettings;
        for( std::map< std::string, std::shared_ptr< RadiationPressureInterfaceSettings > >::iterator
             radiationPressureSettingsIterator = radiationPressureSettings.begin( );
             radiationPressureSettingsIterator != radiationPressureSettings.end( );
             radiationPressureSettingsIterator++ )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRadiationPressureInterface(
                        radiationPressureSettingsIterator->first,
                        createRadiationPressureInterface(
                            radiationPressureSettingsIterator->second,
                            orderedBodySettings.at( i ).first, bodyList ) );
        }

    }

    // Create radiation source model objects for each body (if required).
    // Point sources are created first since paneled sources may depend on them as original sources
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        auto radiationSourceModelSettings = orderedBodySettings.at( i ).second->radiationSourceModelSettings;
        if( radiationSourceModelSettings != nullptr &&
            radiationSourceModelSettings->getRadiationSourceModelType() == RadiationSourceModelType::isotropic_point_source )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRadiationSourceModel(
                        createRadiationSourceModel(
                            orderedBodySettings.at( i ).second->radiationSourceModelSettings,
                            orderedBodySettings.at( i ).first, bodyList ));
        }
    }
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        auto radiationSourceModelSettings = orderedBodySettings.at( i ).second->radiationSourceModelSettings;
        if( radiationSourceModelSettings != nullptr &&
            radiationSourceModelSettings->getRadiationSourceModelType() == RadiationSourceModelType::extended_source )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRadiationSourceModel(
                        createRadiationSourceModel(
                            orderedBodySettings.at( i ).second->radiationSourceModelSettings,
                            orderedBodySettings.at( i ).first, bodyList ));
        }
    }

    // Create radiation pressure target model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->radiationPressureTargetModelSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRadiationPressureTargetModels(
                        createRadiationPressureTargetModel(
                            orderedBodySettings.at( i ).second->radiationPressureTargetModelSettings,
                            orderedBodySettings.at( i ).first, bodyList ) );
        }
    }

    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->bodyDeformationSettings.size( ) > 0 )
        {
            for( unsigned int j = 0; j < orderedBodySettings.at( i ).second->bodyDeformationSettings.size( ); j++ )
            {
                bodyList.at( orderedBodySettings.at( i ).first )->addBodyDeformationModel(
                            createBodyDeformationModel(
                                orderedBodySettings.at( i ).second->bodyDeformationSettings.at( j ),
                                orderedBodySettings.at( i ).first, bodyList ) );
            }
        }
    }


    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        for( unsigned int j = 0; j < orderedBodySettings.at( i ).second->groundStationSettings.size( ); j++ )
        {
            createGroundStation( bodyList.at( orderedBodySettings.at( i ).first ),
                     orderedBodySettings.at( i ).second->groundStationSettings.at( j ) );
        }
    }

    bodyList.processBodyFrameDefinitions< StateScalarType, TimeType >( );

    return bodyList;

}

//! Function to create a simplified system of bodies
/*!
 * Bodies created: Sun, all planets of solar system, Pluto
 * All bodies with Gtop ephemerides and point mass gravity
 * Earth with spherical shape model and simple rotation model
 * @param secondsSinceJ2000 Initial time of the simulation, expressed in seconds since J2000
 * @return List of bodies created
 */
simulation_setup::SystemOfBodies createSimplifiedSystemOfBodies(const double secondsSinceJ2000 = 0 );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEBODIES_H
