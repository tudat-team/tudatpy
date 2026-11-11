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
#include <typeinfo>

#include "tudat/simulation/environment_setup/copySystemOfBodies.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/astro/aerodynamics/comaWindModel.h"
#include "tudat/astro/aerodynamics/customConstantTemperatureAtmosphere.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "tudat/astro/aerodynamics/marsDtmAtmosphereModel.h"
#include "tudat/astro/aerodynamics/panelledAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/astro/aerodynamics/tabulatedAtmosphere.h"
#include "tudat/astro/aerodynamics/windModel.h"
#include "tudat/astro/electromagnetism/luminosityModel.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositionsCircularCoplanar.h"
#include "tudat/astro/ephemerides/compositeEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/customRotationalEphemeris.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/astro/ephemerides/iauRotationModel.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/gravitation/ringGravityField.h"
#include "tudat/astro/gravitation/tabulatedGravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"
#include "tudat/astro/system_models/engineModel.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/basics/timeType.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/environmentModelClosure.h"

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

// Copy a model with its copy constructor if it is of one of the given types (not of a derived type); returns nullptr otherwise
template< typename... ModelTypes, typename BaseType >
std::shared_ptr< BaseType > copyModelOfType( const std::shared_ptr< BaseType > model )
{
    std::shared_ptr< BaseType > modelCopy;
    if( model != nullptr )
    {
        ( ( typeid( *model ) == typeid( ModelTypes ) &&
            ( modelCopy = std::make_shared< ModelTypes >( *std::static_pointer_cast< ModelTypes >( model ) ) ) ) ||
          ... );
    }
    return modelCopy;
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

    // SPICE is not thread-safe, and cannot be used during multi-threaded propagation
    if( isModelOfType< SpiceEphemeris >( ephemeris ) || isModelOfType< TleEphemeris >( ephemeris ) )
    {
        throw std::runtime_error( "Error when copying ephemeris of body " + bodyName +
                                  ", ephemeris calls SPICE, which is not thread-safe. Use, for instance, an interpolated SPICE "
                                  "ephemeris for multi-threaded propagation." );
    }

    // Models without data that changes during propagation are shared
    if( ephemeris == nullptr || isModelOfType< ConstantEphemeris >( ephemeris ) || isModelOfType< ApproximateGtopEphemeris >( ephemeris ) )
    {
        return ephemeris;
    }
    else if( isModelOfType< CustomEphemeris<> >( ephemeris ) )
    {
        printSharedModelWarning( printWarnings, bodyName, "ephemeris", "custom function" );
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

    std::shared_ptr< ScaledEphemeris > scaledEphemeris = std::dynamic_pointer_cast< ScaledEphemeris >( ephemeris );
    if( scaledEphemeris != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "state scaling function of the ephemeris", "custom function" );
        return std::make_shared< ScaledEphemeris >( copyEphemeris( scaledEphemeris->getBaseEphemeris( ), bodyName, printWarnings ),
                                                    scaledEphemeris->getStateScalingFunction( ),
                                                    scaledEphemeris->getIsScalingAbsolute( ) );
    }

    // Approximate JPL ephemerides (which modify members when evaluated) are recreated
    std::shared_ptr< ApproximateJplEphemeris > approximateJplEphemeris = std::dynamic_pointer_cast< ApproximateJplEphemeris >( ephemeris );
    if( approximateJplEphemeris != nullptr )
    {
        return std::make_shared< ApproximateJplEphemeris >( approximateJplEphemeris->getBodyName( ),
                                                            approximateJplEphemeris->getSunGravitationalParameter( ) );
    }
    std::shared_ptr< ApproximateJplCircularCoplanarEphemeris > circularCoplanarEphemeris =
            std::dynamic_pointer_cast< ApproximateJplCircularCoplanarEphemeris >( ephemeris );
    if( circularCoplanarEphemeris != nullptr )
    {
        return std::make_shared< ApproximateJplCircularCoplanarEphemeris >( circularCoplanarEphemeris->getBodyName( ),
                                                                            circularCoplanarEphemeris->getSunGravitationalParameter( ) );
    }

    // Composite ephemerides have no data that changes when evaluated; the functions from which it is composed are shared
    std::shared_ptr< Ephemeris > ephemerisCopy = copyModelOfType< CompositeEphemeris< double, double >,
                                                                  CompositeEphemeris< double, long double >,
                                                                  CompositeEphemeris< Time, double >,
                                                                  CompositeEphemeris< Time, long double > >( ephemeris );
    if( ephemerisCopy != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "functions composing the ephemeris", "composite ephemeris" );
        return ephemerisCopy;
    }

    ephemerisCopy = copyModelOfType< KeplerEphemeris >( ephemeris );
    try
    {
        if( ephemerisCopy == nullptr )
        {
            ephemerisCopy = copyTabulatedEphemeris< double, double >( ephemeris );
        }
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

    // SPICE is not thread-safe, and cannot be used during multi-threaded propagation
    if( isModelOfType< SpiceRotationalEphemeris >( rotationalEphemeris ) )
    {
        throw std::runtime_error( "Error when copying rotation model of body " + bodyName +
                                  ", rotation model calls SPICE, which is not thread-safe. Use, for instance, a tabulated rotation "
                                  "model for multi-threaded propagation." );
    }

    // Models without data that changes during propagation are shared
    if( rotationalEphemeris == nullptr || isModelOfType< SimpleRotationalEphemeris >( rotationalEphemeris ) ||
        isModelOfType< ConstantRotationalEphemeris >( rotationalEphemeris ) )
    {
        return rotationalEphemeris;
    }
    else if( isModelOfType< CustomRotationalEphemeris >( rotationalEphemeris ) )
    {
        printSharedModelWarning( printWarnings, bodyName, "rotation model", "custom function" );
        return rotationalEphemeris;
    }

    // Rotation models with links to other models are copied without closure (which is applied to the copied system of bodies)
    std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotation =
            std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( rotationalEphemeris );
    if( synchronousRotation != nullptr )
    {
        std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotationCopy =
                std::make_shared< SynchronousRotationalEphemeris >( nullptr,
                                                                    synchronousRotation->getCentralBodyName( ),
                                                                    synchronousRotation->getBaseFrameOrientation( ),
                                                                    synchronousRotation->getTargetFrameOrientation( ) );
        std::shared_ptr< LongitudeLibrationCalculator > librationCalculator = synchronousRotation->getLongitudeLibrationCalculator( );
        if( librationCalculator != nullptr )
        {
            std::shared_ptr< LongitudeLibrationCalculator > librationCalculatorCopy =
                    copyModelOfType< DirectLongitudeLibrationCalculator >( librationCalculator );
            if( librationCalculatorCopy == nullptr )
            {
                printSharedModelWarning(
                        printWarnings, bodyName, "libration calculator", "copying this libration calculator type is not implemented" );
                librationCalculatorCopy = librationCalculator;
            }
            synchronousRotationCopy->setLibrationCalculation( librationCalculatorCopy );
        }
        return synchronousRotationCopy;
    }

    std::shared_ptr< AerodynamicAngleRotationalEphemeris > aerodynamicAngleRotation =
            std::dynamic_pointer_cast< AerodynamicAngleRotationalEphemeris >( rotationalEphemeris );
    if( aerodynamicAngleRotation != nullptr )
    {
        // Angle function of trimmed model is recreated by its closure
        std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction;
        if( !aerodynamicAngleRotation->getIsPitchTrimmed( ) && aerodynamicAngleRotation->getAerodynamicAngleFunction( ) != nullptr )
        {
            printSharedModelWarning( printWarnings, bodyName, "aerodynamic angle function of the rotation model", "custom function" );
            aerodynamicAngleFunction = aerodynamicAngleRotation->getAerodynamicAngleFunction( );
        }
        std::shared_ptr< AerodynamicAngleRotationalEphemeris > aerodynamicAngleRotationCopy =
                std::make_shared< AerodynamicAngleRotationalEphemeris >( aerodynamicAngleRotation->getCentralBodyName( ),
                                                                         aerodynamicAngleRotation->getBaseFrameOrientation( ),
                                                                         aerodynamicAngleRotation->getTargetFrameOrientation( ),
                                                                         aerodynamicAngleFunction );
        aerodynamicAngleRotationCopy->setIsPitchTrimmed( aerodynamicAngleRotation->getIsPitchTrimmed( ),
                                                         aerodynamicAngleRotation->getTrimSideslipAndBankAngleFunction( ) );
        return aerodynamicAngleRotationCopy;
    }

    std::shared_ptr< DirectionBasedRotationalEphemeris > directionBasedRotation =
            std::dynamic_pointer_cast< DirectionBasedRotationalEphemeris >( rotationalEphemeris );
    if( directionBasedRotation != nullptr )
    {
        std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculator =
                directionBasedRotation->getInertialBodyAxisDirectionCalculator( );
        std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculatorCopy;
        std::shared_ptr< StateBasedBodyFixedDirectionCalculator > stateBasedDirectionCalculator =
                std::dynamic_pointer_cast< StateBasedBodyFixedDirectionCalculator >( directionCalculator );
        if( stateBasedDirectionCalculator != nullptr )
        {
            directionCalculatorCopy = std::make_shared< StateBasedBodyFixedDirectionCalculator >(
                    stateBasedDirectionCalculator->getCentralBody( ),
                    stateBasedDirectionCalculator->getIsColinearWithVelocity( ),
                    stateBasedDirectionCalculator->getDirectionIsOppositeToVector( ),
                    nullptr );
        }
        else
        {
            directionCalculatorCopy = copyModelOfType< CustomBodyFixedDirectionCalculator >( directionCalculator );
            if( directionCalculatorCopy == nullptr )
            {
                printSharedModelWarning( printWarnings,
                                         bodyName,
                                         "direction calculator of the rotation model",
                                         "copying this calculator type is not implemented" );
                directionCalculatorCopy = directionCalculator;
            }
            else
            {
                printSharedModelWarning( printWarnings, bodyName, "direction function of the rotation model", "custom function" );
                directionCalculatorCopy->applyEnvironmentModelClosure( nullptr );
            }
        }
        directionCalculatorCopy->setDirectionFrame( directionCalculator->getDirectionFrame( ) );

        if( directionBasedRotation->getFreeRotationAngleFunction( ) != nullptr )
        {
            printSharedModelWarning( printWarnings, bodyName, "free rotation angle function of the rotation model", "custom function" );
        }
        return std::make_shared< DirectionBasedRotationalEphemeris >( directionCalculatorCopy,
                                                                      directionBasedRotation->getAssociatedBodyFixedDirection( ),
                                                                      directionBasedRotation->getBaseFrameOrientation( ),
                                                                      directionBasedRotation->getTargetFrameOrientation( ),
                                                                      directionBasedRotation->getFreeRotationAngleFunction( ) );
    }

    std::shared_ptr< RotationalEphemeris > rotationalEphemerisCopy =
            copyModelOfType< IauRotationModel, PlanetaryRotationModel, GcrsToItrsRotationModel >( rotationalEphemeris );

    try
    {
        if( rotationalEphemerisCopy == nullptr )
        {
            rotationalEphemerisCopy = copyTabulatedRotationalEphemeris< double, double >( rotationalEphemeris );
        }
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
    using namespace gravitation;

    // Point-mass field has no data that changes during propagation, and is shared
    if( gravityFieldModel == nullptr )
    {
        return gravityFieldModel;
    }
    const GravityFieldModel& gravityField = *gravityFieldModel;
    if( typeid( gravityField ) == typeid( GravityFieldModel ) )
    {
        return gravityFieldModel;
    }

    // Field variations of time-dependent field are reset by copyGravityFieldVariations
    std::shared_ptr< GravityFieldModel > gravityFieldModelCopy =
            copyModelOfType< SphericalHarmonicsGravityField, TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel );
    if( gravityFieldModelCopy != nullptr )
    {
        return gravityFieldModelCopy;
    }

    // Polyhedron and ring fields are recreated, since their cache would be shared by their copy constructor
    std::shared_ptr< PolyhedronGravityField > polyhedronField = std::dynamic_pointer_cast< PolyhedronGravityField >( gravityFieldModel );
    if( polyhedronField != nullptr )
    {
        return std::make_shared< PolyhedronGravityField >( polyhedronField->getGravitationalParameter( ),
                                                           polyhedronField->getVerticesCoordinates( ),
                                                           polyhedronField->getVerticesDefiningEachFacet( ),
                                                           polyhedronField->getFixedReferenceFrame( ) );
    }
    std::shared_ptr< RingGravityField > ringField = std::dynamic_pointer_cast< RingGravityField >( gravityFieldModel );
    if( ringField != nullptr )
    {
        return std::make_shared< RingGravityField >( ringField->getGravitationalParameter( ),
                                                     ringField->getRingRadius( ),
                                                     ringField->getEllipticIntegralSFromDAndB( ),
                                                     ringField->getFixedReferenceFrame( ) );
    }

    printSharedModelWarning( printWarnings, bodyName, "gravity field model", "copying this gravity field type is not implemented" );
    return gravityFieldModel;
}

// Recreate the gravity field variations of a body, linked to the models in the copied system of bodies
void copyGravityFieldVariations( const std::shared_ptr< Body > body,
                                 const SystemOfBodies& bodiesCopy,
                                 const std::string& bodyName,
                                 const bool printWarnings )
{
    using namespace gravitation;

    std::shared_ptr< GravityFieldVariationsSet > variationSet = body->getGravityFieldVariationSet( );
    std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > gravityFieldCopy =
            std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodiesCopy.at( bodyName )->getGravityFieldModel( ) );
    if( variationSet == nullptr )
    {
        return;
    }
    else if( gravityFieldCopy == nullptr || gravityFieldCopy == body->getGravityFieldModel( ) )
    {
        printSharedModelWarning( printWarnings, bodyName, "gravity field variations", "gravity field model is shared" );
        bodiesCopy.at( bodyName )->setGravityFieldVariationSet( variationSet );
        return;
    }

    std::vector< std::shared_ptr< GravityFieldVariations > > variationObjects = variationSet->getVariationObjects( );
    std::vector< std::shared_ptr< GravityFieldVariations > > variationObjectsCopy;
    for( unsigned int i = 0; i < variationObjects.size( ); i++ )
    {
        std::shared_ptr< GravityFieldVariations > variation = variationObjects.at( i );
        std::shared_ptr< GravityFieldVariations > variationCopy;

        // Tidal variations are recreated without closure (their copy constructor would keep functions bound to the original
        // object), and are linked to the bodies in the copy by closure of the copied system of bodies
        if( isModelOfType< BasicSolidBodyTideGravityFieldVariations >( variation ) )
        {
            std::shared_ptr< BasicSolidBodyTideGravityFieldVariations > tidalVariation =
                    std::static_pointer_cast< BasicSolidBodyTideGravityFieldVariations >( variation );
            std::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > tidalVariationSettings =
                    std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( tidalVariation->getDeformingBodies( ),
                                                                                     tidalVariation->getLoveNumbers( ) );
            tidalVariationSettings->setMeanTidalForcingTerms( tidalVariation->getMeanForcingCosineTerms( ),
                                                              tidalVariation->getMeanForcingSineTerms( ) );
            variationCopy = createGravityFieldVariationsModel( tidalVariationSettings, bodyName, bodiesCopy, false );
        }
        else if( isModelOfType< ModeCoupledSolidBodyTideGravityFieldVariations >( variation ) )
        {
            std::shared_ptr< ModeCoupledSolidBodyTideGravityFieldVariations > tidalVariation =
                    std::static_pointer_cast< ModeCoupledSolidBodyTideGravityFieldVariations >( variation );
            variationCopy =
                    createGravityFieldVariationsModel( std::make_shared< ModeCoupledSolidBodyGravityFieldVariationSettings >(
                                                               tidalVariation->getDeformingBodies( ), tidalVariation->getLoveNumbers( ) ),
                                                       bodyName,
                                                       bodiesCopy,
                                                       false );
        }
        else
        {
            variationCopy =
                    copyModelOfType< PeriodicGravityFieldVariations, PolynomialGravityFieldVariations, TabulatedGravityFieldVariations >(
                            variation );
        }

        if( variationCopy == nullptr )
        {
            printSharedModelWarning(
                    printWarnings, bodyName, "gravity field variation", "copying this gravity field variation type is not implemented" );
            variationCopy = variation;
        }
        variationObjectsCopy.push_back( variationCopy );
    }

    std::shared_ptr< GravityFieldVariationsSet > variationSetCopy = std::make_shared< GravityFieldVariationsSet >( *variationSet );
    variationSetCopy->resetVariationObjects( variationObjectsCopy );
    bodiesCopy.at( bodyName )->setGravityFieldVariationSet( variationSetCopy );

    // Correction functions are updated when processing the frame definitions of the copy
    gravityFieldCopy->setFieldVariationSettings( variationSetCopy, false );
}

std::shared_ptr< aerodynamics::AtmosphereModel > copyAtmosphereModel(
        const std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
        const std::string& bodyName,
        const bool printWarnings )
{
    using namespace aerodynamics;

    if( atmosphereModel == nullptr )
    {
        return atmosphereModel;
    }
    else if( isModelOfType< CustomConstantTemperatureAtmosphere >( atmosphereModel ) ||
             isModelOfType< CustomNumberDensityAtmosphere >( atmosphereModel ) )
    {
        printSharedModelWarning( printWarnings, bodyName, "atmosphere model", "custom function" );
        return atmosphereModel;
    }

    // Scaled atmosphere is recreated from a copy of its base atmosphere; its wind model is copied below
    std::shared_ptr< AtmosphereModel > atmosphereModelCopy;
    std::shared_ptr< ScaledAtmosphereModel > scaledAtmosphere = std::dynamic_pointer_cast< ScaledAtmosphereModel >( atmosphereModel );
    if( scaledAtmosphere != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "density scaling function of the atmosphere", "custom function" );
        atmosphereModelCopy = std::make_shared< ScaledAtmosphereModel >(
                copyAtmosphereModel( scaledAtmosphere->getBaseAtmosphere( ), bodyName, printWarnings ),
                scaledAtmosphere->getDensityScalingFunction( ),
                scaledAtmosphere->getIsScalingAbsolute( ) );
    }

    std::shared_ptr< ExponentialAtmosphere > exponentialAtmosphere = std::dynamic_pointer_cast< ExponentialAtmosphere >( atmosphereModel );
    if( exponentialAtmosphere != nullptr &&
        ( exponentialAtmosphere->getBaseDensityFunction( ) != nullptr || exponentialAtmosphere->getScaleHeightFunction( ) != nullptr ) )
    {
        printSharedModelWarning(
                printWarnings, bodyName, "time-dependent base density/scale height function of the atmosphere", "custom function" );
    }

    // Wind model (constant, empty or custom) has no data that changes during propagation, and is shared by the copy constructors.
    // Links of the coma model to the bodies are set by the closure of the copy
    if( atmosphereModelCopy == nullptr )
    {
        atmosphereModelCopy =
                copyModelOfType< ExponentialAtmosphere, TabulatedAtmosphere, NRLMSISE00Atmosphere, MarsDtmAtmosphereModel, ComaModel >(
                        atmosphereModel );
    }

    if( atmosphereModelCopy == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "atmosphere model", "copying this atmosphere type is not implemented" );
        return atmosphereModel;
    }

    std::shared_ptr< WindModel > windModel = atmosphereModel->getWindModel( );
    if( scaledAtmosphere != nullptr )
    {
        atmosphereModelCopy->setWindModel( windModel );
    }
    if( isModelOfType< ComaWindModel >( windModel ) )
    {
        // Links to the (copied) coma model and the bodies are set by the closure of the copy
        atmosphereModelCopy->setWindModel( copyModelOfType< ComaWindModel >( windModel ) );
    }
    else if( isModelOfType< CustomWindModel >( windModel ) )
    {
        printSharedModelWarning( printWarnings, bodyName, "wind model", "custom function" );
    }
    else if( windModel != nullptr && !isModelOfType< ConstantWindModel >( windModel ) && !isModelOfType< EmptyWindModel >( windModel ) )
    {
        printSharedModelWarning( printWarnings, bodyName, "wind model", "copying this wind model type is not implemented" );
    }
    return atmosphereModelCopy;
}

std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > copyAerodynamicCoefficientInterface(
        const std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface,
        const std::shared_ptr< system_models::VehicleSystems > vehicleSystemsCopy,
        const std::string& bodyName,
        const bool printWarnings )
{
    using namespace aerodynamics;

    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterfaceCopy;

    // Panelled interface is recreated, using the panels of the copied vehicle systems
    std::shared_ptr< PanelledAerodynamicCoefficientInterface > panelledInterface =
            std::dynamic_pointer_cast< PanelledAerodynamicCoefficientInterface >( coefficientInterface );
    std::shared_ptr< ScaledAerodynamicCoefficientInterface > scaledInterface =
            std::dynamic_pointer_cast< ScaledAerodynamicCoefficientInterface >( coefficientInterface );
    if( panelledInterface != nullptr )
    {
        std::shared_ptr< GasSurfaceInteractionModel > gasSurfaceInteractionModel = panelledInterface->getGasSurfaceInteractionModel( );
        std::shared_ptr< GasSurfaceInteractionModel > gasSurfaceInteractionModelCopy =
                createGasSurfaceInteractionModel( gasSurfaceInteractionModel->getGasSurfaceInteractionModelType( ),
                                                  vehicleSystemsCopy->getAllPanels( ),
                                                  gasSurfaceInteractionModel->getReferenceArea( ),
                                                  gasSurfaceInteractionModel->getMaximumNumberOfPixels( ),
                                                  gasSurfaceInteractionModel->getOnlyDrag( ) );
        gasSurfaceInteractionModelCopy->setConstantAerodynamicCoefficients(
                gasSurfaceInteractionModel->getConstantAerodynamicCoefficients( ) );
        coefficientInterfaceCopy =
                std::make_shared< PanelledAerodynamicCoefficientInterface >( gasSurfaceInteractionModelCopy,
                                                                             coefficientInterface->getIndependentVariableNames( ),
                                                                             coefficientInterface->getReferenceArea( ),
                                                                             coefficientInterface->getForceCoefficientsFrame( ) );
    }
    else if( scaledInterface != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "scaling functions of the aerodynamic coefficients", "custom function" );
        coefficientInterfaceCopy = std::make_shared< ScaledAerodynamicCoefficientInterface >(
                copyAerodynamicCoefficientInterface(
                        scaledInterface->getBaseCoefficientInterface( ), vehicleSystemsCopy, bodyName, printWarnings ),
                scaledInterface->getForceCoefficientScalingFunction( ),
                scaledInterface->getMomentCoefficientScalingFunction( ),
                scaledInterface->getIsScalingRelative( ) );
    }
    else
    {
        // Coefficient functions (constant or interpolated) and hypersonic coefficient tables are shared
        coefficientInterfaceCopy =
                copyModelOfType< CustomAerodynamicCoefficientInterface, HypersonicLocalInclinationAnalysis >( coefficientInterface );
        if( coefficientInterfaceCopy == nullptr )
        {
            printSharedModelWarning( printWarnings,
                                     bodyName,
                                     "aerodynamic coefficient interface",
                                     "copying this coefficient interface type is not implemented" );
            return coefficientInterface;
        }
    }

    // Moment contribution of the force is set by the closure of the copy
    coefficientInterfaceCopy->setAddForceContributionToMoments( coefficientInterface->getAddForceContributionToMoments( ) );

    std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > > controlSurfaceInterfaces =
            coefficientInterface->getControlSurfaceIncrementInterfaces( );
    if( controlSurfaceInterfaces.size( ) > 0 )
    {
        for( auto& controlSurfaceInterface : controlSurfaceInterfaces )
        {
            std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceInterfaceCopy =
                    copyModelOfType< CustomControlSurfaceIncrementAerodynamicInterface >( controlSurfaceInterface.second );
            if( controlSurfaceInterfaceCopy == nullptr )
            {
                printSharedModelWarning( printWarnings,
                                         bodyName,
                                         "control surface increment interface " + controlSurfaceInterface.first,
                                         "copying this control surface interface type is not implemented" );
            }
            else
            {
                controlSurfaceInterface.second = controlSurfaceInterfaceCopy;
            }
        }
        coefficientInterfaceCopy->setControlSurfaceIncrements( controlSurfaceInterfaces, false );
    }
    return coefficientInterfaceCopy;
}

std::shared_ptr< RigidBodyProperties > copyMassProperties( const std::shared_ptr< RigidBodyProperties > massProperties,
                                                           const std::string& bodyName,
                                                           const bool printWarnings )
{
    std::shared_ptr< TimeDependentRigidBodyProperties > timeDependentProperties =
            std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties );
    if( timeDependentProperties != nullptr && timeDependentProperties->getMassFunction( ) != nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "mass function", "custom function" );
    }

    std::shared_ptr< RigidBodyProperties > massPropertiesCopy =
            copyModelOfType< TimeDependentRigidBodyProperties, MassDependentRigidBodyProperties >( massProperties );
    if( massPropertiesCopy == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "mass properties", "copying this mass property type is not implemented" );
        return massProperties;
    }
    return massPropertiesCopy;
}

std::shared_ptr< electromagnetism::RadiationSourceModel > copyRadiationSourceModel(
        const std::shared_ptr< electromagnetism::RadiationSourceModel > radiationSourceModel,
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceBodyShapeModelCopy,
        const std::string& bodyName,
        const bool printWarnings )
{
    using namespace electromagnetism;

    // Paneled source is copied without closure (link to original sources), with copies of the radiosity models of its panels
    std::shared_ptr< DynamicallyPaneledRadiationSourceModel > paneledSourceModel =
            std::dynamic_pointer_cast< DynamicallyPaneledRadiationSourceModel >( radiationSourceModel );
    if( paneledSourceModel != nullptr && paneledSourceModel->getPanels( ).size( ) > 0 )
    {
        std::vector< std::unique_ptr< SourcePanelRadiosityModel > > radiosityModels;
        for( auto& radiosityModel : paneledSourceModel->getPanels( ).at( 0 ).getRadiosityModels( ) )
        {
            radiosityModels.push_back( radiosityModel->clone( ) );
        }
        std::shared_ptr< DynamicallyPaneledRadiationSourceModel > paneledSourceModelCopy =
                std::make_shared< DynamicallyPaneledRadiationSourceModel >( sourceBodyShapeModelCopy,
                                                                            nullptr,
                                                                            radiosityModels,
                                                                            paneledSourceModel->getNumberOfPanelsPerRing( ),
                                                                            paneledSourceModel->getSourceName( ) );
        paneledSourceModelCopy->setOriginalSources( paneledSourceModel->getOriginalSourceBodyNames( ),
                                                    paneledSourceModel->getOriginalSourceToSourceOccultingBodies( ) );
        return paneledSourceModelCopy;
    }

    // Statically paneled source is copied without closure, with copies of its panels
    std::shared_ptr< StaticallyPaneledRadiationSourceModel > staticPaneledSourceModel =
            std::dynamic_pointer_cast< StaticallyPaneledRadiationSourceModel >( radiationSourceModel );
    if( staticPaneledSourceModel != nullptr )
    {
        std::vector< RadiationSourcePanel > panelsCopy;
        for( const RadiationSourcePanel& panel : staticPaneledSourceModel->getPanels( ) )
        {
            std::vector< std::unique_ptr< SourcePanelRadiosityModel > > radiosityModels;
            for( auto& radiosityModel : panel.getRadiosityModels( ) )
            {
                radiosityModels.push_back( radiosityModel->clone( ) );
            }
            panelsCopy.emplace_back(
                    panel.getArea( ), panel.getRelativeCenter( ), panel.getSurfaceNormal( ), std::move( radiosityModels ) );
        }
        std::shared_ptr< StaticallyPaneledRadiationSourceModel > staticPaneledSourceModelCopy =
                std::make_shared< StaticallyPaneledRadiationSourceModel >(
                        nullptr, std::move( panelsCopy ), staticPaneledSourceModel->getSourceName( ) );

        // Original sources, retrieved from the radiosity model updater if the source was not created from settings
        std::vector< std::string > originalSourceNames = staticPaneledSourceModel->getOriginalSourceBodyNames( );
        std::map< std::string, std::vector< std::string > > originalSourceToSourceOccultingBodies =
                staticPaneledSourceModel->getOriginalSourceToSourceOccultingBodies( );
        if( originalSourceNames.size( ) == 0 && staticPaneledSourceModel->getSourcePanelRadiosityModelUpdater( ) != nullptr )
        {
            originalSourceNames = staticPaneledSourceModel->getSourcePanelRadiosityModelUpdater( )->getOriginalSourceBodyNames( );
            for( auto occultationModel :
                 staticPaneledSourceModel->getSourcePanelRadiosityModelUpdater( )->getOriginalSourceToSourceOccultationModels( ) )
            {
                originalSourceToSourceOccultingBodies[ occultationModel.first ] = occultationModel.second->getOccultingBodyNames( );
            }
        }
        staticPaneledSourceModelCopy->setOriginalSources( originalSourceNames, originalSourceToSourceOccultingBodies );
        return staticPaneledSourceModelCopy;
    }

    std::shared_ptr< IsotropicPointRadiationSourceModel > pointSourceModel =
            std::dynamic_pointer_cast< IsotropicPointRadiationSourceModel >( radiationSourceModel );
    if( pointSourceModel == nullptr )
    {
        printSharedModelWarning(
                printWarnings, bodyName, "radiation source model", "copying this radiation source model type is not implemented" );
        return radiationSourceModel;
    }

    std::shared_ptr< LuminosityModel > luminosityModel = pointSourceModel->getLuminosityModel( );
    std::shared_ptr< LuminosityModel > luminosityModelCopy = copyModelOfType< ConstantLuminosityModel >( luminosityModel );
    if( luminosityModelCopy == nullptr )
    {
        luminosityModelCopy = copyModelOfType< VariableLuminosityModel >( luminosityModel );
        if( luminosityModelCopy != nullptr )
        {
            printSharedModelWarning( printWarnings, bodyName, "luminosity function", "custom function" );
        }
    }
    if( luminosityModelCopy == nullptr )
    {
        printSharedModelWarning( printWarnings, bodyName, "luminosity model", "copying this luminosity model type is not implemented" );
        luminosityModelCopy = luminosityModel;
    }
    return std::make_shared< IsotropicPointRadiationSourceModel >( luminosityModelCopy, pointSourceModel->getSourceName( ) );
}

std::shared_ptr< system_models::VehicleSystems > copyVehicleSystems( const std::shared_ptr< system_models::VehicleSystems > vehicleSystems,
                                                                     const std::string& bodyName,
                                                                     const bool printWarnings )
{
    using namespace system_models;

    std::shared_ptr< VehicleSystems > vehicleSystemsCopy = std::make_shared< VehicleSystems >( *vehicleSystems );

    // Engines are recreated with copied thrust magnitude models
    for( auto engineModel : vehicleSystems->getEngineModels( ) )
    {
        std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitude = engineModel.second->getThrustMagnitudeWrapper( );
        std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeCopy =
                copyModelOfType< propulsion::ConstantThrustMagnitudeWrapper >( thrustMagnitude );
        if( thrustMagnitudeCopy == nullptr )
        {
            thrustMagnitudeCopy = copyModelOfType< propulsion::CustomThrustMagnitudeWrapper >( thrustMagnitude );
            if( thrustMagnitudeCopy != nullptr )
            {
                printSharedModelWarning(
                        printWarnings, bodyName, "thrust magnitude function of engine " + engineModel.first, "custom function" );
            }
        }
        if( thrustMagnitudeCopy == nullptr )
        {
            // Input variables from the environment are set by the closure of the copy
            thrustMagnitudeCopy = copyModelOfType< propulsion::ParameterizedThrustMagnitudeWrapper >( thrustMagnitude );
            if( thrustMagnitudeCopy != nullptr )
            {
                printSharedModelWarning( printWarnings,
                                         bodyName,
                                         "thrust magnitude and guidance input functions of engine " + engineModel.first,
                                         "custom function" );
            }
        }
        if( thrustMagnitudeCopy == nullptr )
        {
            printSharedModelWarning( printWarnings,
                                     bodyName,
                                     "thrust magnitude model of engine " + engineModel.first,
                                     "copying this thrust magnitude type is not implemented" );
            thrustMagnitudeCopy = thrustMagnitude;
        }
        vehicleSystemsCopy->setEngineModel(
                std::make_shared< EngineModel >(
                        thrustMagnitudeCopy, engineModel.first, engineModel.second->getBodyFixedThrustDirectionFunction( ) ),
                false );
    }

    for( auto referencePoint : vehicleSystems->getReferencePoints( ) )
    {
        vehicleSystemsCopy->setReferencePointPosition(
                referencePoint.first, copyEphemeris( referencePoint.second, bodyName + "_" + referencePoint.first, printWarnings ) );
    }

    // Exterior panels (which store their current orientation) and rotation models of vehicle parts are copied
    std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > exteriorPanels =
            vehicleSystems->getVehicleExteriorPanels( );
    if( exteriorPanels.size( ) > 0 )
    {
        for( auto& partPanels : exteriorPanels )
        {
            for( auto& exteriorPanel : partPanels.second )
            {
                exteriorPanel = std::make_shared< VehicleExteriorPanel >( *exteriorPanel );

                // Surface normal of panel tracking a body is set by closure of vehicle systems
                if( exteriorPanel->getTrackedBody( ) != "" )
                {
                    exteriorPanel->applyEnvironmentModelClosure( nullptr );
                }
            }
        }
        vehicleSystemsCopy->setVehicleExteriorPanels( exteriorPanels );

        std::map< std::string, std::shared_ptr< ephemerides::RotationalEphemeris > > partOrientations =
                vehicleSystems->getVehiclePartOrientation( );
        for( auto& partOrientation : partOrientations )
        {
            partOrientation.second =
                    copyRotationalEphemeris( partOrientation.second, bodyName + "_" + partOrientation.first, printWarnings );
        }
        vehicleSystemsCopy->setVehiclePartOrientation( partOrientations );
        vehicleSystemsCopy->resetTime( );
    }
    return vehicleSystemsCopy;
}

std::shared_ptr< ground_stations::StationMotionModel > copyStationMotionModel(
        const std::shared_ptr< ground_stations::StationMotionModel > stationMotionModel )
{
    using namespace ground_stations;

    std::shared_ptr< CombinedStationMotionModel > combinedStationMotion =
            std::dynamic_pointer_cast< CombinedStationMotionModel >( stationMotionModel );
    if( combinedStationMotion != nullptr )
    {
        std::vector< std::shared_ptr< StationMotionModel > > modelListCopy;
        for( auto model : combinedStationMotion->getModelList( ) )
        {
            modelListCopy.push_back( copyStationMotionModel( model ) );
        }
        return std::make_shared< CombinedStationMotionModel >( modelListCopy );
    }

    // Models linked to bodies are copied (links are set by closure), others have no data that changes during propagation
    std::shared_ptr< StationMotionModel > stationMotionModelCopy =
            copyModelOfType< BodyDeformationStationMotionModel, BodyCentricToBarycentricRelativisticStationMotion >( stationMotionModel );
    return stationMotionModelCopy == nullptr ? stationMotionModel : stationMotionModelCopy;
}

std::shared_ptr< ground_stations::GroundStation > copyGroundStation( const std::shared_ptr< ground_stations::GroundStation > groundStation )
{
    using namespace ground_stations;

    std::shared_ptr< GroundStationState > stationStateCopy =
            std::make_shared< GroundStationState >( *groundStation->getNominalStationState( ) );
    if( stationStateCopy->getStationMotionModel( ) != nullptr )
    {
        stationStateCopy->setStationMotionModel( copyStationMotionModel( stationStateCopy->getStationMotionModel( ) ) );
    }

    std::shared_ptr< GroundStation > groundStationCopy = std::make_shared< GroundStation >(
            stationStateCopy,
            std::make_shared< PointingAnglesCalculator >( nullptr, nullptr ),
            groundStation->getStationId( ),
            groundStation->hasFrequencyCalculator( ) ? groundStation->getTransmittingFrequencyCalculator( ) : nullptr );
    groundStationCopy->setTimingSystem( groundStation->getTimingSystem( ) );
    groundStationCopy->setVehicleSystems( groundStation->getVehicleSystems( ) );
    groundStationCopy->setMeteoData( groundStation->getMeteoData( ) );
    groundStationCopy->setTroposphereData( groundStation->getTroposphereData( ) );
    return groundStationCopy;
}

// Copy (or share) all environment models of a body, without closure (links to other models), which is applied afterwards
void copyBodyModels( const std::shared_ptr< Body > body,
                     const SystemOfBodies& bodiesCopy,
                     const std::string& bodyName,
                     const bool printWarnings )
{
    std::shared_ptr< Body > bodyCopy = bodiesCopy.at( bodyName );
    bodyCopy->setEphemeris( copyEphemeris( body->getEphemeris( ), bodyName, printWarnings ) );
    bodyCopy->setRotationalEphemeris( copyRotationalEphemeris( body->getRotationalEphemeris( ), bodyName, printWarnings ) );
    bodyCopy->setAtmosphereModel( copyAtmosphereModel( body->getAtmosphereModel( ), bodyName, printWarnings ) );

    // Shape models have no data that changes during propagation, and are shared
    bodyCopy->setShapeModel( body->getShapeModel( ) );

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

    if( body->getVehicleSystems( ) != nullptr )
    {
        bodyCopy->setVehicleSystems( copyVehicleSystems( body->getVehicleSystems( ), bodyName, printWarnings ) );
    }

    if( body->getAerodynamicCoefficientInterface( ) != nullptr )
    {
        bodyCopy->setAerodynamicCoefficientInterface( copyAerodynamicCoefficientInterface(
                body->getAerodynamicCoefficientInterface( ), bodyCopy->getVehicleSystems( ), bodyName, printWarnings ) );
    }

    if( body->getRadiationSourceModel( ) != nullptr )
    {
        bodyCopy->setRadiationSourceModel(
                copyRadiationSourceModel( body->getRadiationSourceModel( ), bodyCopy->getShapeModel( ), bodyName, printWarnings ) );
    }

    // Radiation pressure targets (paneled targets without closure: link to panels of vehicle systems)
    std::vector< std::shared_ptr< electromagnetism::RadiationPressureTargetModel > > targetModelsCopy;
    for( auto targetModel : body->getRadiationPressureTargetModels( ) )
    {
        std::shared_ptr< electromagnetism::RadiationPressureTargetModel > targetModelCopy =
                copyModelOfType< electromagnetism::CannonballRadiationPressureTargetModel >( targetModel );
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > paneledTargetModel =
                std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >( targetModel );
        if( paneledTargetModel != nullptr )
        {
            targetModelCopy = std::make_shared< electromagnetism::PaneledRadiationPressureTargetModel >(
                    paneledTargetModel->getSourceToTargetOccultingBodies( ), paneledTargetModel->getMaximumNumberOfPixelsPerSource( ) );
        }
        if( targetModelCopy == nullptr )
        {
            printSharedModelWarning( printWarnings,
                                     bodyName,
                                     "radiation pressure target model",
                                     "copying this radiation pressure target model type is not implemented" );
            targetModelCopy = targetModel;
        }
        targetModelsCopy.push_back( targetModelCopy );
    }
    bodyCopy->setRadiationPressureTargetModels( targetModelsCopy );

    // Flight conditions (links to environment models are reset by closure)
    std::shared_ptr< aerodynamics::FlightConditions > flightConditions = body->getFlightConditions( );
    if( flightConditions != nullptr )
    {
        std::shared_ptr< aerodynamics::FlightConditions > flightConditionsCopy =
                copyModelOfType< aerodynamics::AtmosphericFlightConditions, aerodynamics::FlightConditions >( flightConditions );
        if( flightConditionsCopy == nullptr )
        {
            printSharedModelWarning(
                    printWarnings, bodyName, "flight conditions", "copying this flight conditions type is not implemented" );
            flightConditionsCopy = flightConditions;
        }
        else
        {
            flightConditionsCopy->applyEnvironmentModelClosure( nullptr, nullptr );
        }
        bodyCopy->setFlightConditions( flightConditionsCopy );
    }
    // Ground stations and deformation models (links to the bodies are set by the closure of the copy)
    for( auto groundStation : body->getGroundStationMap( ) )
    {
        bodyCopy->addGroundStation( groundStation.first, copyGroundStation( groundStation.second ) );
    }
    for( auto deformationModel : body->getBodyDeformationModels( ) )
    {
        std::shared_ptr< basic_astrodynamics::BodyDeformationModel > deformationModelCopy =
                copyModelOfType< basic_astrodynamics::Iers2010EarthDeformation, basic_astrodynamics::BasicTidalBodyDeformation >(
                        deformationModel );
        // Pole and ocean tide models have no data that changes during propagation
        bodyCopy->addBodyDeformationModel( deformationModelCopy == nullptr ? deformationModel : deformationModelCopy );
    }
    bodyCopy->setIonosphereModel( body->getIonosphereModel( ) );
    bodyCopy->setClimateModel( body->getClimateModel( ) );
    bodyCopy->setTimeScaleConverter( body->getTimeScaleConverter( ) );
}

}  // namespace

SystemOfBodies copySystemOfBodiesWithoutFrameProcessing( const SystemOfBodies& bodies, const bool printWarnings )
{
    SystemOfBodies bodiesCopy( bodies.getFrameOrigin( ), bodies.getFrameOrientation( ) );

    // Space-time properties of the copy, with its own PPN parameters (metric is recreated after closure of the copied bodies)
    std::shared_ptr< SpaceTimeProperties > spaceTimeProperties = bodies.getSpaceTimeProperties( );
    bodiesCopy.setSpaceTimeProperties( std::make_shared< SpaceTimeProperties >(
            std::make_shared< relativity::PPNParameterSet >( *spaceTimeProperties->getPpnParameterSet( ) ),
            spaceTimeProperties->getEquivalencePrincipleLpiViolationParameter( ) ) );
    bodiesCopy.getSpaceTimeProperties( )->setBaseMetricSettings( spaceTimeProperties->getBaseMetricSettings( ) );

    for( auto bodyIterator : bodies.getMap( ) )
    {
        bodiesCopy.createEmptyBody( bodyIterator.first, false );
    }

    // Copy models (without closure)
    for( auto bodyIterator : bodies.getMap( ) )
    {
        copyBodyModels( bodyIterator.second, bodiesCopy, bodyIterator.first, printWarnings );
    }

    // Gravity field variations (which need the copied gravity fields of all bodies)
    for( auto bodyIterator : bodies.getMap( ) )
    {
        copyGravityFieldVariations( bodyIterator.second, bodiesCopy, bodyIterator.first, printWarnings );
    }

    // Apply closure of copied environment models (links to models in the copied system of bodies)
    for( auto bodyIterator : bodiesCopy.getMap( ) )
    {
        registerEnvironmentModelClosures( bodyIterator.second );
    }
    std::map< std::string, std::vector< EnvironmentModelClosureType > > modelsWithoutClosure = closeSystemOfBodies( bodiesCopy, false );
    std::map< std::string, std::vector< EnvironmentModelClosureType > > originalModelsWithoutClosure =
            getEnvironmentModelsWithoutClosure( bodies );
    for( auto bodyModelsWithoutClosure : modelsWithoutClosure )
    {
        for( EnvironmentModelClosureType closureType : bodyModelsWithoutClosure.second )
        {
            bool isOriginalWithoutClosure = originalModelsWithoutClosure.count( bodyModelsWithoutClosure.first ) > 0 &&
                    std::find( originalModelsWithoutClosure.at( bodyModelsWithoutClosure.first ).begin( ),
                               originalModelsWithoutClosure.at( bodyModelsWithoutClosure.first ).end( ),
                               closureType ) != originalModelsWithoutClosure.at( bodyModelsWithoutClosure.first ).end( );
            if( !isOriginalWithoutClosure )
            {
                throw std::runtime_error( "Error when copying system of bodies, closure of environment model (type " +
                                          std::to_string( closureType ) + ") of body " + bodyModelsWithoutClosure.first +
                                          " could not be applied to the copy" );
            }
        }
    }

    // Custom aerodynamic angle interfaces are shared with the copy (other interfaces are set by closure, or by the propagation
    // setup for the angles from the body)
    for( auto bodyIterator : bodies.getMap( ) )
    {
        std::shared_ptr< aerodynamics::FlightConditions > flightConditions = bodyIterator.second->getFlightConditions( );
        std::shared_ptr< aerodynamics::FlightConditions > flightConditionsCopy =
                bodiesCopy.at( bodyIterator.first )->getFlightConditions( );
        if( flightConditions == nullptr || flightConditionsCopy == flightConditions ||
            flightConditions->getAerodynamicAngleCalculator( ) == nullptr )
        {
            continue;
        }
        std::shared_ptr< reference_frames::BodyFixedAerodynamicAngleInterface > angleInterface =
                flightConditions->getAerodynamicAngleCalculator( )->getBodyFixedAngleInterface( );
        if( angleInterface != nullptr && angleInterface->getAngleSource( ) == reference_frames::custom_body_fixed_angles )
        {
            printSharedModelWarning( printWarnings, bodyIterator.first, "aerodynamic angle interface", "custom function" );
            if( flightConditionsCopy->getAerodynamicAngleCalculator( )->getBodyFixedAngleInterface( ) == nullptr )
            {
                flightConditionsCopy->getAerodynamicAngleCalculator( )->setBodyFixedAngleInterface( angleInterface );
            }
        }
    }

    // Space-time metric of the copy (linked to the copied bodies)
    if( spaceTimeProperties->getBaseMetricSettings( ) != nullptr )
    {
        bodiesCopy.getSpaceTimeProperties( )->setBaseMetric(
                createSpaceTimeMetric( spaceTimeProperties->getBaseMetricSettings( ), bodiesCopy ) );
    }
    else if( spaceTimeProperties->getBaseMetric( ) != nullptr )
    {
        printSharedModelWarning( printWarnings, "(global)", "space-time metric", "metric was not created from settings" );
        bodiesCopy.getSpaceTimeProperties( )->setBaseMetric( spaceTimeProperties->getBaseMetric( ) );
    }

    return bodiesCopy;
}

}  // namespace simulation_setup

}  // namespace tudat
