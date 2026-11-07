/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <typeinfo>

#include "tudat/simulation/environment_setup/copySystemOfBodies.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/basics/timeType.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

void printSharedModelWarning( const bool printWarnings,
                              const std::string& bodyName,
                              const std::string& modelType,
                              const std::string& reason )
{
    if( printWarnings )
    {
        std::cerr << "Warning when copying system of bodies, the " << modelType << " of body " << bodyName
                  << " is shared between the original and the copy (" << reason
                  << "). It must not be used in multi-threaded propagation, unless it is thread-safe." << std::endl;
    }
}

template< typename StateScalarType, typename TimeType >
std::shared_ptr< ephemerides::Ephemeris > copyTabulatedEphemeris( const std::shared_ptr< ephemerides::Ephemeris > ephemeris )
{
    std::shared_ptr< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > > tabulatedEphemeris =
            std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >( ephemeris );
    if( tabulatedEphemeris == nullptr )
    {
        return nullptr;
    }
    typename ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType >::StateInterpolatorPointer interpolatorCopy;
    if( tabulatedEphemeris->getInterpolator( ) != nullptr )
    {
        interpolatorCopy = tabulatedEphemeris->getInterpolator( )->clone( );
    }
    return std::make_shared< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
            interpolatorCopy, tabulatedEphemeris->getReferenceFrameOrigin( ), tabulatedEphemeris->getReferenceFrameOrientation( ) );
}

std::shared_ptr< ephemerides::Ephemeris > copyEphemeris( const std::shared_ptr< ephemerides::Ephemeris > ephemeris,
                                                         const std::string& bodyName,
                                                         const bool printWarnings )
{
    using namespace ephemerides;

    // Models without data that changes during propagation are shared
    if( ephemeris == nullptr || std::dynamic_pointer_cast< ConstantEphemeris >( ephemeris ) != nullptr ||
        std::dynamic_pointer_cast< SpiceEphemeris >( ephemeris ) != nullptr )
    {
        return ephemeris;
    }

    std::shared_ptr< MultiArcEphemeris > multiArcEphemeris = std::dynamic_pointer_cast< MultiArcEphemeris >( ephemeris );
    if( multiArcEphemeris != nullptr )
    {
        std::vector< std::shared_ptr< Ephemeris > > arcEphemeridesCopy;
        for( std::shared_ptr< Ephemeris > arcEphemeris : multiArcEphemeris->getSingleArcEphemerides( ) )
        {
            arcEphemeridesCopy.push_back( copyEphemeris( arcEphemeris, bodyName, printWarnings ) );
        }
        std::shared_ptr< MultiArcEphemeris > multiArcEphemerisCopy = std::make_shared< MultiArcEphemeris >(
                std::map< double, std::shared_ptr< Ephemeris > >( ),
                multiArcEphemeris->getReferenceFrameOrigin( ),
                multiArcEphemeris->getReferenceFrameOrientation( ),
                copyEphemeris( multiArcEphemeris->getDefaultEphemeris( ), bodyName, printWarnings ) );
        multiArcEphemerisCopy->resetSingleArcEphemerides( arcEphemeridesCopy, multiArcEphemeris->getArcStartTimes( ) );
        return multiArcEphemerisCopy;
    }

    std::shared_ptr< Ephemeris > ephemerisCopy;
    try
    {
        ephemerisCopy = copyTabulatedEphemeris< double, double >( ephemeris );
        if( ephemerisCopy == nullptr )
        {
            ephemerisCopy = copyTabulatedEphemeris< long double, double >( ephemeris );
        }
        if( ephemerisCopy == nullptr )
        {
            ephemerisCopy = copyTabulatedEphemeris< double, Time >( ephemeris );
        }
        if( ephemerisCopy == nullptr )
        {
            ephemerisCopy = copyTabulatedEphemeris< long double, Time >( ephemeris );
        }
    }
    catch( const std::runtime_error& )
    {
        printSharedModelWarning( printWarnings, bodyName, "ephemeris", "copying its interpolator type is not implemented" );
        return ephemeris;
    }

    if( ephemerisCopy == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "ephemeris", "copying this ephemeris type is not implemented" );
        return ephemeris;
    }
    return ephemerisCopy;
}

template< typename StateScalarType, typename TimeType >
std::shared_ptr< ephemerides::RotationalEphemeris > copyTabulatedRotationalEphemeris(
        const std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris )
{
    std::shared_ptr< ephemerides::TabulatedRotationalEphemeris< StateScalarType, TimeType > > tabulatedEphemeris =
            std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< StateScalarType, TimeType > >( rotationalEphemeris );
    if( tabulatedEphemeris == nullptr )
    {
        return nullptr;
    }
    else if( tabulatedEphemeris->getInterpolator( ) == nullptr )
    {
        throw std::runtime_error( "Error when copying tabulated rotation model, no interpolator found" );
    }
    return std::make_shared< ephemerides::TabulatedRotationalEphemeris< StateScalarType, TimeType > >(
            tabulatedEphemeris->getInterpolator( )->clone( ),
            tabulatedEphemeris->getBaseFrameOrientation( ),
            tabulatedEphemeris->getTargetFrameOrientation( ) );
}

std::shared_ptr< ephemerides::RotationalEphemeris > copyRotationalEphemeris(
        const std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        const std::string& bodyName,
        const bool printWarnings )
{
    using namespace ephemerides;

    // Models without data that changes during propagation are shared
    if( rotationalEphemeris == nullptr || std::dynamic_pointer_cast< SimpleRotationalEphemeris >( rotationalEphemeris ) != nullptr ||
        std::dynamic_pointer_cast< SpiceRotationalEphemeris >( rotationalEphemeris ) != nullptr )
    {
        return rotationalEphemeris;
    }

    std::shared_ptr< RotationalEphemeris > rotationalEphemerisCopy;
    try
    {
        rotationalEphemerisCopy = copyTabulatedRotationalEphemeris< double, double >( rotationalEphemeris );
        if( rotationalEphemerisCopy == nullptr )
        {
            rotationalEphemerisCopy = copyTabulatedRotationalEphemeris< long double, double >( rotationalEphemeris );
        }
        if( rotationalEphemerisCopy == nullptr )
        {
            rotationalEphemerisCopy = copyTabulatedRotationalEphemeris< double, Time >( rotationalEphemeris );
        }
        if( rotationalEphemerisCopy == nullptr )
        {
            rotationalEphemerisCopy = copyTabulatedRotationalEphemeris< long double, Time >( rotationalEphemeris );
        }
    }
    catch( const std::runtime_error& )
    {
        printSharedModelWarning( printWarnings, bodyName, "rotation model", "copying its interpolator type is not implemented" );
        return rotationalEphemeris;
    }

    if( rotationalEphemerisCopy == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "rotation model", "copying this rotation model type is not implemented" );
        return rotationalEphemeris;
    }
    return rotationalEphemerisCopy;
}

std::shared_ptr< gravitation::GravityFieldModel > copyGravityFieldModel(
        const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel,
        const std::string& bodyName,
        const bool printWarnings )
{
    // Point-mass field has no data that changes during propagation, and is shared
    if( gravityFieldModel == nullptr || typeid( *gravityFieldModel ) == typeid( gravitation::GravityFieldModel ) )
    {
        return gravityFieldModel;
    }
    else if( typeid( *gravityFieldModel ) == typeid( gravitation::SphericalHarmonicsGravityField ) )
    {
        std::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicsField =
                std::static_pointer_cast< gravitation::SphericalHarmonicsGravityField >( gravityFieldModel );
        return std::make_shared< gravitation::SphericalHarmonicsGravityField >( sphericalHarmonicsField->getGravitationalParameter( ),
                                                                                sphericalHarmonicsField->getReferenceRadius( ),
                                                                                sphericalHarmonicsField->getCosineCoefficients( ),
                                                                                sphericalHarmonicsField->getSineCoefficients( ),
                                                                                sphericalHarmonicsField->getFixedReferenceFrame( ),
                                                                                sphericalHarmonicsField->getScaledMeanMomentOfInertia( ) );
    }
    printSharedModelWarning( printWarnings, bodyName, "gravity field model", "copying this gravity field type is not implemented" );
    return gravityFieldModel;
}

std::shared_ptr< aerodynamics::AtmosphereModel > copyAtmosphereModel(
        const std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
        const std::string& bodyName,
        const bool printWarnings )
{
    if( atmosphereModel == nullptr )
    {
        return atmosphereModel;
    }
    else if( typeid( *atmosphereModel ) == typeid( aerodynamics::ExponentialAtmosphere ) )
    {
        std::shared_ptr< aerodynamics::ExponentialAtmosphere > exponentialAtmosphere =
                std::static_pointer_cast< aerodynamics::ExponentialAtmosphere >( atmosphereModel );
        if( exponentialAtmosphere->getBaseDensityFunction( ) != nullptr || exponentialAtmosphere->getScaleHeightFunction( ) != nullptr )
        {
            printSharedModelWarning(
                    printWarnings, bodyName, "time-dependent base density/scale height function of the atmosphere", "custom function" );
        }
        return std::make_shared< aerodynamics::ExponentialAtmosphere >( *exponentialAtmosphere );
    }
    printSharedModelWarning( printWarnings, bodyName, "atmosphere model", "copying this atmosphere type is not implemented" );
    return atmosphereModel;
}

std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > copyAerodynamicCoefficientInterface(
        const std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface,
        const std::string& bodyName,
        const bool printWarnings )
{
    // Only constant coefficients (custom interface without independent variables) can be copied
    if( typeid( *coefficientInterface ) == typeid( aerodynamics::CustomAerodynamicCoefficientInterface ) &&
        coefficientInterface->getNumberOfIndependentVariables( ) == 0 &&
        coefficientInterface->getMomentContributionInterface( ) == nullptr && coefficientInterface->getControlSurfaceNames( ).size( ) == 0 )
    {
        return std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                *std::static_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >( coefficientInterface ) );
    }
    printSharedModelWarning( printWarnings,
                             bodyName,
                             "aerodynamic coefficient interface",
                             "copying is only implemented for constant coefficients without moment contribution or control surfaces" );
    return coefficientInterface;
}

std::shared_ptr< RigidBodyProperties > copyMassProperties( const std::shared_ptr< RigidBodyProperties > massProperties,
                                                           const std::string& bodyName,
                                                           const bool printWarnings )
{
    if( typeid( *massProperties ) == typeid( TimeDependentRigidBodyProperties ) )
    {
        std::shared_ptr< TimeDependentRigidBodyProperties > timeDependentProperties =
                std::static_pointer_cast< TimeDependentRigidBodyProperties >( massProperties );
        if( timeDependentProperties->getMassFunction( ) != nullptr )
        {
            printSharedModelWarning( printWarnings, bodyName, "mass function", "custom function" );
        }
        return std::make_shared< TimeDependentRigidBodyProperties >( *timeDependentProperties );
    }
    printSharedModelWarning( printWarnings, bodyName, "mass properties", "copying this mass property type is not implemented" );
    return massProperties;
}

// Copy (or share) all environment models of a body that do not require links to other models in the copied system of bodies
void copyUnlinkedBodyModels( const std::shared_ptr< Body > body,
                             const std::shared_ptr< Body > bodyCopy,
                             const std::string& bodyName,
                             const bool printWarnings )
{
    bodyCopy->setEphemeris( copyEphemeris( body->getEphemeris( ), bodyName, printWarnings ) );
    bodyCopy->setRotationalEphemeris( copyRotationalEphemeris( body->getRotationalEphemeris( ), bodyName, printWarnings ) );
    bodyCopy->setAtmosphereModel( copyAtmosphereModel( body->getAtmosphereModel( ), bodyName, printWarnings ) );

    // Shape models without data that changes during propagation are shared
    std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel = body->getShapeModel( );
    if( shapeModel != nullptr && std::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( shapeModel ) == nullptr &&
        std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ) == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "shape model", "copying this shape model type is not implemented" );
    }
    bodyCopy->setShapeModel( shapeModel );

    // Mass properties computed from the gravity field are recreated by Body::setGravityFieldModel
    if( body->getMassProperties( ) != nullptr &&
        std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( body->getMassProperties( ) ) == nullptr )
    {
        bodyCopy->setMassProperties( copyMassProperties( body->getMassProperties( ), bodyName, printWarnings ) );
    }
    if( body->getGravityFieldModel( ) != nullptr )
    {
        bodyCopy->setGravityFieldModel( copyGravityFieldModel( body->getGravityFieldModel( ), bodyName, printWarnings ) );
    }

    if( body->getGravityFieldVariationSet( ) != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "gravity field variations", "copying is not implemented" );
        bodyCopy->setGravityFieldVariationSet( body->getGravityFieldVariationSet( ) );
    }

    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = body->getVehicleSystems( );
    if( vehicleSystems->getEngineModels( ).size( ) == 0 && vehicleSystems->getVehicleExteriorPanels( ).size( ) == 0 &&
        vehicleSystems->getReferencePoints( ).size( ) == 0 )
    {
        bodyCopy->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( *vehicleSystems ) );
    }
    else
    {
        printSharedModelWarning(
                printWarnings, bodyName, "vehicle systems", "copying is not implemented for engines, exterior panels or reference points" );
        bodyCopy->setVehicleSystems( vehicleSystems );
    }

    if( body->getAerodynamicCoefficientInterface( ) != nullptr )
    {
        bodyCopy->setAerodynamicCoefficientInterface(
                copyAerodynamicCoefficientInterface( body->getAerodynamicCoefficientInterface( ), bodyName, printWarnings ) );
    }

    if( body->getRadiationSourceModel( ) != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "radiation source model", "copying is not implemented" );
        bodyCopy->setRadiationSourceModel( body->getRadiationSourceModel( ) );
    }
    if( body->getRadiationPressureTargetModels( ).size( ) > 0 )
    {
        printSharedModelWarning( printWarnings, bodyName, "radiation pressure target models", "copying is not implemented" );
        bodyCopy->setRadiationPressureTargetModels( body->getRadiationPressureTargetModels( ) );
    }

    // Models used for observations (evaluated outside of the propagation) are shared
    for( auto groundStation : body->getGroundStationMap( ) )
    {
        bodyCopy->addGroundStation( groundStation.first, groundStation.second );
    }
    for( auto deformationModel : body->getBodyDeformationModels( ) )
    {
        bodyCopy->addBodyDeformationModel( deformationModel );
    }
    bodyCopy->setIonosphereModel( body->getIonosphereModel( ) );
    bodyCopy->setClimateModel( body->getClimateModel( ) );
    bodyCopy->setTimeScaleConverter( body->getTimeScaleConverter( ) );
}

// Recreate the flight conditions of a body, linked to the models in the copied system of bodies
void copyFlightConditions( const std::shared_ptr< Body > body,
                           const SystemOfBodies& bodiesCopy,
                           const std::string& bodyName,
                           const bool printWarnings )
{
    std::shared_ptr< aerodynamics::FlightConditions > flightConditions = body->getFlightConditions( );
    if( flightConditions == nullptr )
    {
        return;
    }

    std::string centralBodyName = flightConditions->getCentralBody( );
    std::shared_ptr< aerodynamics::FlightConditions > flightConditionsCopy;
    if( typeid( *flightConditions ) == typeid( aerodynamics::AtmosphericFlightConditions ) )
    {
        flightConditionsCopy =
                createAtmosphericFlightConditions( bodiesCopy.at( bodyName ), bodiesCopy.at( centralBodyName ), bodyName, centralBodyName );
    }
    else
    {
        flightConditionsCopy =
                createFlightConditions( bodiesCopy.at( bodyName ), bodiesCopy.at( centralBodyName ), bodyName, centralBodyName );
    }
    bodiesCopy.at( bodyName )->setFlightConditions( flightConditionsCopy );

    // Angle interfaces other than the one obtained from the body's rotation model (set by Body::setFlightConditions) are shared
    std::shared_ptr< reference_frames::BodyFixedAerodynamicAngleInterface > angleInterface =
            flightConditions->getAerodynamicAngleCalculator( )->getBodyFixedAngleInterface( );
    if( angleInterface != nullptr &&
        std::dynamic_pointer_cast< reference_frames::FromGenericEphemerisAerodynamicAngleInterface >( angleInterface ) == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "aerodynamic angle interface", "custom or linked aerodynamic angles" );
        flightConditionsCopy->getAerodynamicAngleCalculator( )->setBodyFixedAngleInterface( angleInterface );
    }
}

}  // namespace

SystemOfBodies copySystemOfBodiesWithoutFrameProcessing( const SystemOfBodies& bodies, const bool printWarnings )
{
    SystemOfBodies bodiesCopy( bodies.getFrameOrigin( ), bodies.getFrameOrientation( ) );
    bodiesCopy.setSpaceTimeProperties( bodies.getSpaceTimeProperties( ) );

    for( auto bodyIterator : bodies.getMap( ) )
    {
        bodiesCopy.createEmptyBody( bodyIterator.first, false );
    }

    for( auto bodyIterator : bodies.getMap( ) )
    {
        copyUnlinkedBodyModels( bodyIterator.second, bodiesCopy.at( bodyIterator.first ), bodyIterator.first, printWarnings );
    }

    for( auto bodyIterator : bodies.getMap( ) )
    {
        copyFlightConditions( bodyIterator.second, bodiesCopy, bodyIterator.first, printWarnings );
    }

    return bodiesCopy;
}

}  // namespace simulation_setup

}  // namespace tudat
