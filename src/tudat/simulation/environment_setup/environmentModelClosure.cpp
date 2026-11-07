/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/environmentModelClosure.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "tudat/astro/aerodynamics/panelledAerodynamicCoefficientInterface.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/customRotationalEphemeris.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/astro/ephemerides/iauRotationModel.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"
#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/gravitation/tabulatedGravityFieldVariations.h"
#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"
#include "tudat/astro/system_models/engineModel.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/basics/timeType.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

template< typename ModelType, typename BaseType >
bool isModelOfType( const std::shared_ptr< BaseType > model )
{
    return std::dynamic_pointer_cast< ModelType >( model ) != nullptr;
}

void throwUnknownModelTypeError( const std::string& modelType )
{
    throw std::runtime_error( "Error when determining whether " + modelType +
                              " requires closure, type of model is not known (add it in environmentModelClosure.cpp)" );
}

}  // namespace

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel )
{
    using namespace ephemerides;

    if( rotationModel == nullptr )
    {
        return false;
    }
    else if( isModelOfType< SynchronousRotationalEphemeris >( rotationModel ) ||
             isModelOfType< AerodynamicAngleRotationalEphemeris >( rotationModel ) )
    {
        return true;
    }
    else if( isModelOfType< DirectionBasedRotationalEphemeris >( rotationModel ) )
    {
        std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculator =
                std::dynamic_pointer_cast< DirectionBasedRotationalEphemeris >( rotationModel )->getInertialBodyAxisDirectionCalculator( );
        return isModelOfType< StateBasedBodyFixedDirectionCalculator >( directionCalculator ) ||
                directionCalculator->getDirectionFrame( ).first != inertial_satellite_based_frame;
    }
    else if( isModelOfType< SimpleRotationalEphemeris >( rotationModel ) || isModelOfType< SpiceRotationalEphemeris >( rotationModel ) ||
             isModelOfType< ConstantRotationalEphemeris >( rotationModel ) || isModelOfType< CustomRotationalEphemeris >( rotationModel ) ||
             isModelOfType< IauRotationModel >( rotationModel ) || isModelOfType< PlanetaryRotationModel >( rotationModel ) ||
             isModelOfType< GcrsToItrsRotationModel >( rotationModel ) ||
             isModelOfType< TabulatedRotationalEphemeris< double, double > >( rotationModel ) ||
             isModelOfType< TabulatedRotationalEphemeris< long double, double > >( rotationModel ) ||
             isModelOfType< TabulatedRotationalEphemeris< double, Time > >( rotationModel ) ||
             isModelOfType< TabulatedRotationalEphemeris< long double, Time > >( rotationModel ) )
    {
        return false;
    }
    throwUnknownModelTypeError( "rotation model" );
    return false;
}

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< gravitation::GravityFieldVariationsSet > variationSet )
{
    using namespace gravitation;

    if( variationSet == nullptr )
    {
        return false;
    }

    bool requiresClosure = false;
    for( std::shared_ptr< GravityFieldVariations > variation : variationSet->getVariationObjects( ) )
    {
        if( isModelOfType< SolidBodyTideGravityFieldVariations >( variation ) )
        {
            requiresClosure = true;
        }
        else if( !isModelOfType< PeriodicGravityFieldVariations >( variation ) &&
                 !isModelOfType< PolynomialGravityFieldVariations >( variation ) &&
                 !isModelOfType< TabulatedGravityFieldVariations >( variation ) )
        {
            throwUnknownModelTypeError( "gravity field variation" );
        }
    }
    return requiresClosure;
}

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions )
{
    return flightConditions != nullptr;
}

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface )
{
    using namespace aerodynamics;

    if( coefficientInterface == nullptr )
    {
        return false;
    }
    else if( isModelOfType< PanelledAerodynamicCoefficientInterface >( coefficientInterface ) )
    {
        return true;
    }
    else if( isModelOfType< ScaledAerodynamicCoefficientInterface >( coefficientInterface ) )
    {
        return coefficientInterface->getAddForceContributionToMoments( ) ||
                doesEnvironmentModelRequireClosure(
                        std::dynamic_pointer_cast< ScaledAerodynamicCoefficientInterface >( coefficientInterface )
                                ->getBaseCoefficientInterface( ) );
    }
    else if( isModelOfType< CustomAerodynamicCoefficientInterface >( coefficientInterface ) ||
             isModelOfType< HypersonicLocalInclinationAnalysis >( coefficientInterface ) )
    {
        return coefficientInterface->getAddForceContributionToMoments( );
    }
    throwUnknownModelTypeError( "aerodynamic coefficient interface" );
    return false;
}

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< system_models::VehicleSystems > vehicleSystems )
{
    if( vehicleSystems == nullptr )
    {
        return false;
    }
    for( std::shared_ptr< system_models::VehicleExteriorPanel > panel : vehicleSystems->getAllPanels( ) )
    {
        if( panel->getTrackedBody( ) != "" )
        {
            return true;
        }
    }
    for( auto engineModel : vehicleSystems->getEngineModels( ) )
    {
        if( isModelOfType< propulsion::ParameterizedThrustMagnitudeWrapper >( engineModel.second->getThrustMagnitudeWrapper( ) ) )
        {
            return true;
        }
    }
    return false;
}

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< electromagnetism::RadiationSourceModel > radiationSourceModel )
{
    using namespace electromagnetism;

    if( radiationSourceModel == nullptr || isModelOfType< IsotropicPointRadiationSourceModel >( radiationSourceModel ) )
    {
        return false;
    }
    else if( isModelOfType< PaneledRadiationSourceModel >( radiationSourceModel ) )
    {
        return true;
    }
    throwUnknownModelTypeError( "radiation source model" );
    return false;
}

bool doesEnvironmentModelRequireClosure( const std::shared_ptr< electromagnetism::RadiationPressureTargetModel > targetModel )
{
    using namespace electromagnetism;

    if( targetModel == nullptr || isModelOfType< CannonballRadiationPressureTargetModel >( targetModel ) )
    {
        return false;
    }
    else if( isModelOfType< PaneledRadiationPressureTargetModel >( targetModel ) )
    {
        return true;
    }
    throwUnknownModelTypeError( "radiation pressure target model" );
    return false;
}

void registerEnvironmentModelClosures( const std::shared_ptr< Body > body )
{
    if( doesEnvironmentModelRequireClosure( body->getRotationalEphemeris( ) ) )
    {
        body->setEnvironmentModelRequiresClosure( rotation_model_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getGravityFieldVariationSet( ) ) )
    {
        body->setEnvironmentModelRequiresClosure( gravity_field_variations_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getFlightConditions( ) ) )
    {
        body->setEnvironmentModelRequiresClosure( flight_conditions_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getAerodynamicCoefficientInterface( ) ) )
    {
        body->setEnvironmentModelRequiresClosure( aerodynamic_coefficients_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getVehicleSystems( ) ) )
    {
        body->setEnvironmentModelRequiresClosure( vehicle_systems_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getRadiationSourceModel( ) ) )
    {
        body->setEnvironmentModelRequiresClosure( radiation_source_closure );
    }
    for( auto targetModel : body->getRadiationPressureTargetModels( ) )
    {
        if( doesEnvironmentModelRequireClosure( targetModel ) )
        {
            body->setEnvironmentModelRequiresClosure( radiation_pressure_target_closure );
        }
    }
}

bool applyEnvironmentModelClosure( const EnvironmentModelClosureType closureType,
                                   const std::string& bodyName,
                                   const SystemOfBodies& bodies )
{
    std::shared_ptr< Body > body = bodies.at( bodyName );
    switch( closureType )
    {
        case rotation_model_closure:
            return applyEnvironmentModelClosure( body->getRotationalEphemeris( ), bodyName, bodies );
        case gravity_field_variations_closure:
            return applyEnvironmentModelClosure( body->getGravityFieldVariationSet( ), bodyName, bodies );
        case flight_conditions_closure:
            return applyEnvironmentModelClosure( body->getFlightConditions( ), bodyName, bodies );
        case aerodynamic_coefficients_closure:
            return applyEnvironmentModelClosure( body->getAerodynamicCoefficientInterface( ), bodyName, bodies );
        case vehicle_systems_closure:
            return applyEnvironmentModelClosure( body->getVehicleSystems( ), bodyName, bodies );
        case radiation_source_closure:
            return applyEnvironmentModelClosure( body->getRadiationSourceModel( ), bodyName, bodies );
        case radiation_pressure_target_closure: {
            bool isClosureApplied = true;
            for( auto targetModel : body->getRadiationPressureTargetModels( ) )
            {
                isClosureApplied = applyEnvironmentModelClosure( targetModel, bodyName, bodies ) && isClosureApplied;
            }
            return isClosureApplied;
        }
        default:
            throw std::runtime_error( "Error when applying closure of environment model of body " + bodyName + ", closure type " +
                                      std::to_string( closureType ) + " not recognized" );
    }
}

std::map< std::string, std::vector< EnvironmentModelClosureType > > getEnvironmentModelsWithoutClosure( const SystemOfBodies& bodies )
{
    std::map< std::string, std::vector< EnvironmentModelClosureType > > modelsWithoutClosure;
    for( auto bodyIterator : bodies.getMap( ) )
    {
        for( auto closureState : bodyIterator.second->getEnvironmentModelClosureStates( ) )
        {
            if( !closureState.second )
            {
                modelsWithoutClosure[ bodyIterator.first ].push_back( closureState.first );
            }
        }
    }
    return modelsWithoutClosure;
}

std::map< std::string, std::vector< EnvironmentModelClosureType > > closeSystemOfBodies( const SystemOfBodies& bodies )
{
    bool isAnyClosureApplied = true;
    while( isAnyClosureApplied )
    {
        isAnyClosureApplied = false;
        for( auto modelsWithoutClosure : getEnvironmentModelsWithoutClosure( bodies ) )
        {
            for( EnvironmentModelClosureType closureType : modelsWithoutClosure.second )
            {
                if( applyEnvironmentModelClosure( closureType, modelsWithoutClosure.first, bodies ) )
                {
                    bodies.at( modelsWithoutClosure.first )->setEnvironmentModelClosureApplied( closureType );
                    isAnyClosureApplied = true;
                }
            }
        }
    }
    return getEnvironmentModelsWithoutClosure( bodies );
}

}  // namespace simulation_setup

}  // namespace tudat
