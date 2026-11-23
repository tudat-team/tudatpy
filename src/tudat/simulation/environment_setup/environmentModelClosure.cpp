/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

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

void throwUnknownModelTypeError( const std::string& modelType )
{
    throw std::runtime_error( "Error when determining whether " + modelType +
                              " requires closure, type of model is not known (add it in environmentModelClosure.cpp)" );
}

//! Function to determine whether a rotation model requires closure (does not check whether it has been applied)
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

//! Function to determine whether a set of gravity field variations requires closure (does not check whether it has been applied)
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

//! Function to determine whether flight conditions require closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions )
{
    return flightConditions != nullptr;
}

//! Function to determine whether an aerodynamic coefficient interface requires closure (does not check whether it has been applied)
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

//! Function to determine whether vehicle systems require closure (does not check whether it has been applied)
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
    return false;
}

//! Function to determine whether a radiation source model requires closure (does not check whether it has been applied)
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

//! Function to determine whether a radiation pressure target model requires closure (does not check whether it has been applied)
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

//! Function to apply the closure of a single type of environment model of a body. Returns false if the closure cannot yet be
//! applied (environment models it depends on are not yet closed, or do not yet exist).
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

//! Function to determine which environment models (types) of a body require closure
std::vector< EnvironmentModelClosureType > getEnvironmentModelClosureTypes( const std::shared_ptr< Body > body )
{
    std::vector< EnvironmentModelClosureType > closureTypes;
    if( doesEnvironmentModelRequireClosure( body->getRotationalEphemeris( ) ) )
    {
        closureTypes.push_back( rotation_model_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getGravityFieldVariationSet( ) ) )
    {
        closureTypes.push_back( gravity_field_variations_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getFlightConditions( ) ) )
    {
        closureTypes.push_back( flight_conditions_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getAerodynamicCoefficientInterface( ) ) )
    {
        closureTypes.push_back( aerodynamic_coefficients_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getVehicleSystems( ) ) )
    {
        closureTypes.push_back( vehicle_systems_closure );
    }
    if( doesEnvironmentModelRequireClosure( body->getRadiationSourceModel( ) ) )
    {
        closureTypes.push_back( radiation_source_closure );
    }
    for( auto targetModel : body->getRadiationPressureTargetModels( ) )
    {
        if( doesEnvironmentModelRequireClosure( targetModel ) )
        {
            closureTypes.push_back( radiation_pressure_target_closure );
        }
    }
    // Closure types of models of which a body may have several are added once
    closureTypes.erase( std::unique( closureTypes.begin( ), closureTypes.end( ) ), closureTypes.end( ) );
    return closureTypes;
}

}  // namespace

void closeEnvironmentModels( const std::string& bodyName,
                             const SystemOfBodies& bodies,
                             const std::vector< EnvironmentModelClosureType >& closureTypes )
{
    std::vector< EnvironmentModelClosureType > closuresToApply;
    for( EnvironmentModelClosureType closureType : getEnvironmentModelClosureTypes( bodies.at( bodyName ) ) )
    {
        if( closureTypes.size( ) == 0 || std::find( closureTypes.begin( ), closureTypes.end( ), closureType ) != closureTypes.end( ) )
        {
            closuresToApply.push_back( closureType );
        }
    }

    // Closures are applied repeatedly, since some closures require other closures (of the same body) to be applied first
    bool isAnyClosureApplied = true;
    while( isAnyClosureApplied )
    {
        isAnyClosureApplied = false;
        for( auto closureIterator = closuresToApply.begin( ); closureIterator != closuresToApply.end( ); )
        {
            if( applyEnvironmentModelClosure( *closureIterator, bodyName, bodies ) )
            {
                closureIterator = closuresToApply.erase( closureIterator );
                isAnyClosureApplied = true;
            }
            else
            {
                closureIterator++;
            }
        }
    }
}

void closeSystemOfBodies( const SystemOfBodies& bodies )
{
    for( auto bodyIterator : bodies.getMap( ) )
    {
        closeEnvironmentModels( bodyIterator.first, bodies );
    }
}

}  // namespace simulation_setup

}  // namespace tudat
