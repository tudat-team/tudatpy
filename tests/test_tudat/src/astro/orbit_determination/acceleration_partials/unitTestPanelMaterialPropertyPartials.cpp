/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Validation of the estimation of the four panel (macro-model) gas-surface interaction material
 *    properties: energy / normal / tangential accommodation coefficients and the normal velocity at
 *    wall ratio. Three layers are tested:
 *      1. The analytical partial of the body-frame force coefficient vector w.r.t. each property
 *         (GasSurfaceInteractionModel::computeAerodynamicCoefficientsPartial) against a central
 *         finite difference of the forward model, for the Storch, Sentman and Cook models across a
 *         sweep of incoming directions and coefficient values.
 *      2. The full analytical aerodynamic-acceleration partial (AerodynamicAccelerationPartial)
 *         against a numerical finite difference, for a panelled vehicle in atmospheric flight.
 *      3. The behaviour of the PanelMaterialPropertyParameter estimatable-parameter class itself.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/aerodynamics/gasSurfaceInteractionModel.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/panelMaterialProperty.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::estimatable_parameters;
using namespace tudat::system_models;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_panel_material_property_partials )

//! Compare two 3-vectors with a tolerance that is scaled by the overall vector magnitude, so that
//! components that are (near-)zero in both vectors do not trigger spurious fractional-tolerance failures.
void checkVectorClose( const Eigen::Vector3d& analytical,
                       const Eigen::Vector3d& numerical,
                       const double tolerance,
                       const double absoluteTolerance = 1.0E-12 )
{
    const double scale = std::max( numerical.cwiseAbs( ).maxCoeff( ), analytical.cwiseAbs( ).maxCoeff( ) );
    for( int k = 0; k < 3; k++ )
    {
        BOOST_CHECK_SMALL( std::fabs( analytical( k ) - numerical( k ) ), tolerance * scale + absoluteTolerance );
    }
}

//! Create a single body-fixed exterior panel with geometry defined and all four material properties set.
std::shared_ptr< VehicleExteriorPanel > makeMaterialPanel( const Eigen::Vector3d& surfaceNormal,
                                                           const double panelArea,
                                                           const double panelTemperature,
                                                           const double energyAccommodation,
                                                           const double normalAccommodation,
                                                           const double tangentialAccommodation,
                                                           const double normalVelocityRatio,
                                                           const std::string& panelTypeId )
{
    std::function< Eigen::Vector3d( ) > normalFunction = [ = ]( ) { return surfaceNormal; };
    std::function< Eigen::Vector3d( ) > positionFunction = [ = ]( ) { return surfaceNormal; };
    // A (non-degenerate) triangle is required to mark the panel geometry as defined; with self-shadowing
    // disabled (maximumNumberOfPixels = 0) the triangle does not enter the force computation.
    Triangle3d triangle3d( Eigen::Vector3d( 0.0, 0.0, 0.0 ), Eigen::Vector3d( 0.0, 1.0, 0.0 ), Eigen::Vector3d( 0.0, 0.0, 1.0 ) );

    std::shared_ptr< VehicleExteriorPanel > panel = std::make_shared< VehicleExteriorPanel >(
            normalFunction, positionFunction, panelArea, panelTemperature, "", nullptr, triangle3d, Eigen::Vector3d::Zero( ), true );
    panel->setEnergyAccomodationCoefficient( energyAccommodation );
    panel->setNormalAccomodationCoefficient( normalAccommodation );
    panel->setTangentialAccomodationCoefficient( tangentialAccommodation );
    panel->setNormalVelocityAtWallRatio( normalVelocityRatio );
    panel->setPanelTypeId( panelTypeId );
    panel->updatePanel( Eigen::Quaterniond::Identity( ) );
    return panel;
}

//! Set the material property of all panels in a group on a gas-surface interaction model and re-evaluate it.
Eigen::Vector3d evaluateCoefficientsForPropertyValue( const std::shared_ptr< GasSurfaceInteractionModel > model,
                                                      const std::vector< std::shared_ptr< VehicleExteriorPanel > >& panels,
                                                      const PanelMaterialPropertyType propertyType,
                                                      const double value,
                                                      const Eigen::Vector3d& incomingDirection )
{
    for( auto panel : panels )
    {
        switch( propertyType )
        {
            case energy_accommodation_property:
                panel->setEnergyAccomodationCoefficient( value );
                break;
            case normal_accommodation_property:
                panel->setNormalAccomodationCoefficient( value );
                break;
            case tangential_accommodation_property:
                panel->setTangentialAccomodationCoefficient( value );
                break;
            case normal_velocity_ratio_property:
                panel->setNormalVelocityAtWallRatio( value );
                break;
        }
    }
    model->setIncomingDirection( incomingDirection );
    return model->computeAerodynamicCoefficients( );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Layer 1: GSI analytical coefficient partial vs finite difference.
//////////////////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( testGasSurfaceInteractionCoefficientPartials )
{
    const double PI = mathematical_constants::PI;
    const std::string groupId = "Bus";
    const double panelArea = 0.5;
    const double referenceArea = 1.0;
    const double panelTemperature = 300.0;
    const double freeStreamTemperature = 500.0;
    const double airSpeed = 3.5E3;
    const double specificGasConstant = 400.0;

    // A back-facing panel (normal along +x) and a forward-facing panel (normal along -z) so that both
    // an illuminated and a shadowed panel are exercised within the same group.
    std::vector< Eigen::Vector3d > panelNormals = { Eigen::Vector3d( 1.0, 0.0, 0.0 ), Eigen::Vector3d( 0.0, 0.0, -1.0 ) };

    // Sweep over incoming directions, including grazing (PI/2) cases.
    std::vector< double > angleOfAttack = { 0.0, PI / 10.0, PI / 4.0, PI / 3.0, PI / 2.0 };
    std::vector< double > angleOfSideslip = { 0.0, PI / 8.0, PI / 3.0 };

    // Nominal material property values, including small and large coefficients to probe the regime robustness.
    std::vector< double > nominalCoefficientValues = { 0.05, 0.5, 0.95 };

    for( const double nominalValue : nominalCoefficientValues )
    {
        std::vector< std::shared_ptr< VehicleExteriorPanel > > panels;
        for( const Eigen::Vector3d& normal : panelNormals )
        {
            panels.push_back( makeMaterialPanel(
                    normal, panelArea, panelTemperature, nominalValue, nominalValue, nominalValue, nominalValue, groupId ) );
        }

        std::map< GasSurfaceInteractionModelType, std::vector< PanelMaterialPropertyType > > modelProperties = {
            { storch, { normal_accommodation_property, tangential_accommodation_property, normal_velocity_ratio_property } },
            { sentman, { energy_accommodation_property } },
            { cook, { energy_accommodation_property } }
        };

        for( const auto& modelEntry : modelProperties )
        {
            std::shared_ptr< GasSurfaceInteractionModel > model =
                    createGasSurfaceInteractionModel( modelEntry.first, panels, referenceArea, 0, false );
            if( modelEntry.first == sentman || modelEntry.first == cook )
            {
                model->setFreeStreamTemperature( freeStreamTemperature );
            }
            if( modelEntry.first == sentman )
            {
                model->setAirSpeed( airSpeed );
                model->setSpecifiGasConstant( specificGasConstant );
            }

            for( const double aoa : angleOfAttack )
            {
                for( const double aos : angleOfSideslip )
                {
                    Eigen::Vector3d incomingDirection(
                            -std::cos( aoa ) * std::cos( aos ), -std::sin( aos ), -std::sin( aoa ) * std::cos( aos ) );

                    for( const PanelMaterialPropertyType propertyType : modelEntry.second )
                    {
                        // Forward evaluation at nominal value to set the current geometry state, then analytical partial.
                        evaluateCoefficientsForPropertyValue( model, panels, propertyType, nominalValue, incomingDirection );
                        Eigen::Vector3d analyticalPartial = model->computeAerodynamicCoefficientsPartial( propertyType, groupId );

                        // Central finite difference of the forward model w.r.t. the property.
                        const double perturbation = 1.0E-6;
                        Eigen::Vector3d upperturbed = evaluateCoefficientsForPropertyValue(
                                model, panels, propertyType, nominalValue + perturbation, incomingDirection );
                        Eigen::Vector3d downperturbed = evaluateCoefficientsForPropertyValue(
                                model, panels, propertyType, nominalValue - perturbation, incomingDirection );
                        Eigen::Vector3d numericalPartial = ( upperturbed - downperturbed ) / ( 2.0 * perturbation );

                        // Restore the nominal value.
                        evaluateCoefficientsForPropertyValue( model, panels, propertyType, nominalValue, incomingDirection );

                        checkVectorClose( analyticalPartial, numericalPartial, 1.0E-6, 1.0E-10 );
                    }

                    // Properties that do not enter a given model must yield an exactly zero partial.
                    std::vector< PanelMaterialPropertyType > allProperties = { energy_accommodation_property,
                                                                               normal_accommodation_property,
                                                                               tangential_accommodation_property,
                                                                               normal_velocity_ratio_property };
                    evaluateCoefficientsForPropertyValue( model, panels, energy_accommodation_property, nominalValue, incomingDirection );
                    for( const PanelMaterialPropertyType propertyType : allProperties )
                    {
                        if( std::find( modelEntry.second.begin( ), modelEntry.second.end( ), propertyType ) == modelEntry.second.end( ) )
                        {
                            Eigen::Vector3d unusedPartial = model->computeAerodynamicCoefficientsPartial( propertyType, groupId );
                            BOOST_CHECK_SMALL( unusedPartial.norm( ), std::numeric_limits< double >::epsilon( ) );
                        }
                    }

                    // A partial w.r.t. a non-existent panel group must be exactly zero.
                    Eigen::Vector3d emptyGroupPartial =
                            model->computeAerodynamicCoefficientsPartial( modelEntry.second.front( ), "DoesNotExist" );
                    BOOST_CHECK_SMALL( emptyGroupPartial.norm( ), std::numeric_limits< double >::epsilon( ) );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testGasSurfaceInteractionPartialGroupSelectionAndModelOptions )
{
    const std::string busGroupId = "Bus";
    const std::string solarArrayGroupId = "SolarArray";
    const double referenceArea = 4.0;
    const double freeStreamTemperature = 500.0;
    const Eigen::Vector3d incomingDirection = Eigen::Vector3d( -1.0, -0.5, 0.2 ).normalized( );

    std::shared_ptr< VehicleExteriorPanel > busPanel =
            makeMaterialPanel( Eigen::Vector3d::UnitX( ), 1.2, 290.0, 0.35, 0.55, 0.65, 0.75, busGroupId );
    std::shared_ptr< VehicleExteriorPanel > solarArrayPanel =
            makeMaterialPanel( Eigen::Vector3d::UnitY( ), 2.1, 340.0, 0.80, 0.45, 0.60, 0.70, solarArrayGroupId );
    std::vector< std::shared_ptr< VehicleExteriorPanel > > allPanels = { busPanel, solarArrayPanel };
    std::vector< std::shared_ptr< VehicleExteriorPanel > > busPanels = { busPanel };
    std::vector< std::shared_ptr< VehicleExteriorPanel > > solarArrayPanels = { solarArrayPanel };

    std::shared_ptr< GasSurfaceInteractionModel > cookModel = createGasSurfaceInteractionModel( cook, allPanels, referenceArea, 0, false );
    cookModel->setFreeStreamTemperature( freeStreamTemperature );
    evaluateCoefficientsForPropertyValue( cookModel, busPanels, energy_accommodation_property, 0.35, incomingDirection );

    const Eigen::Vector3d busAnalyticalPartial =
            cookModel->computeAerodynamicCoefficientsPartial( energy_accommodation_property, busGroupId );
    const Eigen::Vector3d solarArrayAnalyticalPartial =
            cookModel->computeAerodynamicCoefficientsPartial( energy_accommodation_property, solarArrayGroupId );
    BOOST_CHECK( busAnalyticalPartial.norm( ) > 0.0 );
    BOOST_CHECK( solarArrayAnalyticalPartial.norm( ) > 0.0 );
    BOOST_CHECK( ( busAnalyticalPartial - solarArrayAnalyticalPartial ).norm( ) > 1.0E-6 );

    // The analytical partial must use the illumination fraction cached by the most recent forward evaluation.
    cookModel->getIlluminatedPanelFractions( ).at( 0 ) = 0.25;
    const Eigen::Vector3d partiallyIlluminatedBusPartial =
            cookModel->computeAerodynamicCoefficientsPartial( energy_accommodation_property, busGroupId );
    checkVectorClose( partiallyIlluminatedBusPartial, 0.25 * busAnalyticalPartial, 1.0E-14 );
    cookModel->getIlluminatedPanelFractions( ).at( 0 ) = 1.0;

    const double perturbation = 1.0E-6;
    const Eigen::Vector3d busUp = evaluateCoefficientsForPropertyValue(
            cookModel, busPanels, energy_accommodation_property, 0.35 + perturbation, incomingDirection );
    const Eigen::Vector3d busDown = evaluateCoefficientsForPropertyValue(
            cookModel, busPanels, energy_accommodation_property, 0.35 - perturbation, incomingDirection );
    checkVectorClose( busAnalyticalPartial, ( busUp - busDown ) / ( 2.0 * perturbation ), 1.0E-6, 1.0E-10 );

    const Eigen::Vector3d solarArrayUp = evaluateCoefficientsForPropertyValue(
            cookModel, solarArrayPanels, energy_accommodation_property, 0.80 + perturbation, incomingDirection );
    const Eigen::Vector3d solarArrayDown = evaluateCoefficientsForPropertyValue(
            cookModel, solarArrayPanels, energy_accommodation_property, 0.80 - perturbation, incomingDirection );
    checkVectorClose( solarArrayAnalyticalPartial, ( solarArrayUp - solarArrayDown ) / ( 2.0 * perturbation ), 1.0E-6, 1.0E-10 );

    // Exercise the projection applied when a gas-surface interaction model is configured to return drag only.
    std::shared_ptr< GasSurfaceInteractionModel > dragOnlyCookModel =
            createGasSurfaceInteractionModel( cook, allPanels, referenceArea, 0, true );
    dragOnlyCookModel->setFreeStreamTemperature( freeStreamTemperature );
    evaluateCoefficientsForPropertyValue( dragOnlyCookModel, busPanels, energy_accommodation_property, 0.35, incomingDirection );
    const Eigen::Vector3d dragOnlyAnalyticalPartial =
            dragOnlyCookModel->computeAerodynamicCoefficientsPartial( energy_accommodation_property, busGroupId );
    const Eigen::Vector3d dragOnlyUp = evaluateCoefficientsForPropertyValue(
            dragOnlyCookModel, busPanels, energy_accommodation_property, 0.35 + perturbation, incomingDirection );
    const Eigen::Vector3d dragOnlyDown = evaluateCoefficientsForPropertyValue(
            dragOnlyCookModel, busPanels, energy_accommodation_property, 0.35 - perturbation, incomingDirection );
    checkVectorClose( dragOnlyAnalyticalPartial, ( dragOnlyUp - dragOnlyDown ) / ( 2.0 * perturbation ), 1.0E-6, 1.0E-10 );
    BOOST_CHECK_SMALL( dragOnlyAnalyticalPartial.cross( incomingDirection ).norm( ), 1.0E-14 );

    // Models that do not consume panel material properties inherit an exact zero partial.
    std::shared_ptr< GasSurfaceInteractionModel > newtonModel =
            createGasSurfaceInteractionModel( newton, allPanels, referenceArea, 0, false );
    newtonModel->setIncomingDirection( incomingDirection );
    newtonModel->computeAerodynamicCoefficients( );
    BOOST_CHECK_SMALL( newtonModel->computeAerodynamicCoefficientsPartial( energy_accommodation_property, busGroupId ).norm( ),
                       std::numeric_limits< double >::epsilon( ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Layer 2: full aerodynamic acceleration partial vs numerical finite difference.
//////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
//! Reset and re-evaluate the vehicle flight conditions (used as the environment update for the numerical partial).
void updateFlightConditionsWithPerturbedState( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions,
                                               const double timeToUpdate )
{
    flightConditions->resetCurrentTime( );
    flightConditions->updateConditions( timeToUpdate );
}
}  // namespace

BOOST_AUTO_TEST_CASE( testAerodynamicAccelerationMaterialPropertyPartials )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::string groupId = "Bus";

    std::vector< std::pair< EstimatebleParametersEnum, double > > allMaterialParameters = { { energy_accommodation_coefficient, 1.0E-5 },
                                                                                            { normal_accommodation_coefficient, 1.0E-5 },
                                                                                            { tangential_accommodation_coefficient,
                                                                                              1.0E-5 },
                                                                                            { normal_velocity_at_wall_ratio, 1.0E-5 } };
    std::map< GasSurfaceInteractionModelType, std::vector< EstimatebleParametersEnum > > activeModelParameters = {
        { storch, { normal_accommodation_coefficient, tangential_accommodation_coefficient, normal_velocity_at_wall_ratio } },
        { sentman, { energy_accommodation_coefficient } },
        { cook, { energy_accommodation_coefficient } }
    };
    std::vector< Eigen::Vector3d > componentScalingCases = { Eigen::Vector3d::Ones( ), Eigen::Vector3d( 1.35, 0.75, 1.60 ) };

    for( const auto& modelEntry : activeModelParameters )
    {
        for( const Eigen::Vector3d& componentScalings : componentScalingCases )
        {
            // Build the environment freshly for each gas-surface interaction model and scaling case.
            BodyListSettings defaultBodySettings = getDefaultBodySettings( { "Earth" } );
            defaultBodySettings.at( "Earth" )->ephemerisSettings =
                    std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
            SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

            const double vehicleMass = 5.0E3;
            const double referenceArea = 4.0;
            bodies.createEmptyBody( "Vehicle" );
            bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );
            bodies.at( "Vehicle" )
                    ->setRotationalEphemeris(
                            std::make_shared< SimpleRotationalEphemeris >( 0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );

            std::vector< Eigen::Vector3d > panelNormals = { Eigen::Vector3d( 1.0, 0.2, -0.1 ).normalized( ),
                                                            Eigen::Vector3d( -0.3, 1.0, 0.4 ).normalized( ),
                                                            Eigen::Vector3d( 0.5, -0.2, 1.0 ).normalized( ),
                                                            Eigen::Vector3d( -1.0, -0.5, 0.3 ).normalized( ) };
            std::vector< double > panelAreas = { 1.2, 0.7, 1.5, 0.9 };
            std::vector< double > panelTemperatures = { 290.0, 315.0, 340.0, 305.0 };
            std::vector< std::shared_ptr< VehicleExteriorPanel > > panels;
            for( unsigned int i = 0; i < panelNormals.size( ); i++ )
            {
                panels.push_back( makeMaterialPanel(
                        panelNormals.at( i ), panelAreas.at( i ), panelTemperatures.at( i ), 0.85, 0.9, 0.8, 0.95, groupId ) );
            }
            std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > panelMap;
            panelMap[ groupId ] = panels;
            std::shared_ptr< VehicleSystems > vehicleSystems = std::make_shared< VehicleSystems >( vehicleMass );
            vehicleSystems->setVehicleExteriorPanels( panelMap );
            bodies.at( "Vehicle" )->setVehicleSystems( vehicleSystems );

            std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                    std::make_shared< PanelledAerodynamicCoefficientSettings >(
                            modelEntry.first, referenceArea, 0, false, body_fixed_frame_coefficients );
            bodies.at( "Vehicle" )
                    ->setAerodynamicCoefficientInterface(
                            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );

            // Spherical entry state at 120 km altitude.
            Eigen::Vector6d vehicleSphericalEntryState;
            vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
                    spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
            vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
            vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
            vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
            vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) = -0.9 * mathematical_constants::PI / 180.0;
            vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;
            Eigen::Vector6d systemInitialState =
                    tudat::ephemerides::transformStateToTargetFrame( convertSphericalOrbitalToCartesianState( vehicleSphericalEntryState ),
                                                                     0.0,
                                                                     bodies.at( "Earth" )->getRotationalEphemeris( ) );

            bodies.at( "Earth" )->setStateFromEphemeris( 0.0 );
            bodies.at( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( 0.0 );
            bodies.at( "Vehicle" )->setState( systemInitialState );

            std::shared_ptr< AccelerationModel3d > accelerationModel =
                    simulation_setup::createAerodynamicAcceleratioModel( bodies.at( "Vehicle" ), bodies.at( "Earth" ), "Vehicle", "Earth" );
            std::shared_ptr< AerodynamicAcceleration > aerodynamicAcceleration =
                    std::dynamic_pointer_cast< AerodynamicAcceleration >( accelerationModel );
            BOOST_REQUIRE( aerodynamicAcceleration != nullptr );
            aerodynamicAcceleration->setDragComponentScaling( componentScalings( 0 ) );
            aerodynamicAcceleration->setSideComponentScaling( componentScalings( 1 ) );
            aerodynamicAcceleration->setLiftComponentScaling( componentScalings( 2 ) );

            bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( 0.0 );
            accelerationModel->updateMembers( 0.0 );

            std::shared_ptr< AccelerationPartial > aerodynamicAccelerationPartial =
                    createAnalyticalAccelerationPartial( accelerationModel,
                                                         std::make_pair( "Vehicle", bodies.at( "Vehicle" ) ),
                                                         std::make_pair( "Earth", bodies.at( "Earth" ) ),
                                                         bodies );
            aerodynamicAccelerationPartial->update( 0.0 );

            std::function< void( ) > environmentUpdateFunction =
                    std::bind( &updateFlightConditionsWithPerturbedState, bodies.at( "Vehicle" )->getFlightConditions( ), 0.0 );

            for( const auto& parameterEntry : allMaterialParameters )
            {
                BOOST_TEST_CONTEXT( "GSI model " << modelEntry.first << ", parameter " << parameterEntry.first << ", scalings "
                                                 << componentScalings.transpose( ) )
                {
                    std::shared_ptr< EstimatableParameter< double > > parameter = createDoubleParameterToEstimate< double, double >(
                            std::make_shared< EstimatableParameterSettings >( "Vehicle", parameterEntry.first, groupId ), bodies );

                    Eigen::Vector3d analyticalPartial = aerodynamicAccelerationPartial->wrtParameter( parameter );
                    Eigen::Vector3d numericalPartial = calculateAccelerationWrtParameterPartials(
                            parameter, accelerationModel, parameterEntry.second, environmentUpdateFunction );

                    const bool isActiveParameter =
                            std::find( modelEntry.second.begin( ), modelEntry.second.end( ), parameterEntry.first ) !=
                            modelEntry.second.end( );
                    if( isActiveParameter )
                    {
                        BOOST_CHECK( numericalPartial.norm( ) > 1.0E-12 );
                        checkVectorClose( analyticalPartial, numericalPartial, 1.0E-5 );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( analyticalPartial.norm( ), 1.0E-14 );
                        BOOST_CHECK_SMALL( numericalPartial.norm( ), 1.0E-12 );
                    }
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Layer 3: PanelMaterialPropertyParameter behaviour.
//////////////////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( testPanelMaterialPropertyParameterClass )
{
    const std::string groupId = "Bus";

    // Consistent group: exercise the getter/setter dispatch for every supported property.
    {
        std::vector< std::shared_ptr< VehicleExteriorPanel > > panels = {
            makeMaterialPanel( Eigen::Vector3d::UnitX( ), 1.0, 300.0, 0.1, 0.2, 0.3, 0.4, groupId ),
            makeMaterialPanel( Eigen::Vector3d::UnitY( ), 1.0, 300.0, 0.1, 0.2, 0.3, 0.4, groupId )
        };
        struct PropertyTestCase {
            EstimatebleParametersEnum parameterType;
            double initialValue;
            double updatedValue;
            std::function< double( const std::shared_ptr< VehicleExteriorPanel >& ) > getter;
        };
        const std::vector< PropertyTestCase > propertyTestCases = {
            { energy_accommodation_coefficient, 0.1, 0.51, []( const auto& panel ) { return panel->getEnergyAccomodationCoefficient( ); } },
            { normal_accommodation_coefficient, 0.2, 0.52, []( const auto& panel ) { return panel->getNormalAccomodationCoefficient( ); } },
            { tangential_accommodation_coefficient,
              0.3,
              0.53,
              []( const auto& panel ) { return panel->getTangentialAccomodationCoefficient( ); } },
            { normal_velocity_at_wall_ratio, 0.4, 0.54, []( const auto& panel ) { return panel->getNormalVelocityAtWallRatio( ); } }
        };

        for( const PropertyTestCase& propertyTestCase : propertyTestCases )
        {
            PanelMaterialPropertyParameter parameter( panels, "Vehicle", groupId, propertyTestCase.parameterType );
            BOOST_CHECK( isDoubleParameter( propertyTestCase.parameterType ) );
            BOOST_CHECK_EQUAL( parameter.getParameterSize( ), 1 );
            BOOST_CHECK( parameter.getParameterDescription( ).find( groupId ) != std::string::npos );
            BOOST_CHECK_CLOSE( parameter.getParameterValue( ), propertyTestCase.initialValue, 1.0E-12 );

            parameter.setParameterValue( propertyTestCase.updatedValue );
            BOOST_CHECK_CLOSE( parameter.getParameterValue( ), propertyTestCase.updatedValue, 1.0E-12 );
            for( const auto& panel : panels )
            {
                BOOST_CHECK_CLOSE( propertyTestCase.getter( panel ), propertyTestCase.updatedValue, 1.0E-12 );
            }
        }
    }

    // Inconsistent group: normalizeValue resets every panel to the average value.
    {
        std::vector< std::shared_ptr< VehicleExteriorPanel > > panels = {
            makeMaterialPanel( Eigen::Vector3d::UnitX( ), 1.0, 300.0, 0.5, 0.6, 0.5, 0.5, groupId ),
            makeMaterialPanel( Eigen::Vector3d::UnitY( ), 1.0, 300.0, 0.5, 0.8, 0.5, 0.5, groupId )
        };

        PanelMaterialPropertyParameter parameter( panels, "Vehicle", groupId, normal_accommodation_coefficient );
        BOOST_CHECK_CLOSE( parameter.getParameterValue( ), 0.7, 1.0E-12 );
        for( auto panel : panels )
        {
            BOOST_CHECK_CLOSE( panel->getNormalAccomodationCoefficient( ), 0.7, 1.0E-12 );
        }
    }

    // Constructor must reject an unsupported parameter type.
    {
        std::vector< std::shared_ptr< VehicleExteriorPanel > > panels = { makeMaterialPanel(
                Eigen::Vector3d::UnitX( ), 1.0, 300.0, 0.7, 0.7, 0.7, 0.7, groupId ) };
        BOOST_CHECK_THROW( PanelMaterialPropertyParameter( panels, "Vehicle", groupId, specular_reflectivity ), std::runtime_error );
    }

    // Constructor must reject an empty panel list.
    {
        std::vector< std::shared_ptr< VehicleExteriorPanel > > emptyPanels;
        BOOST_CHECK_THROW( PanelMaterialPropertyParameter( emptyPanels, "Vehicle", groupId, energy_accommodation_coefficient ),
                           std::runtime_error );
    }
}

BOOST_AUTO_TEST_CASE( testPanelMaterialPropertyParameterFactoryGroupSelection )
{
    const std::string busGroupId = "Bus";
    const std::string solarArrayGroupId = "SolarArray";
    std::vector< std::shared_ptr< VehicleExteriorPanel > > busPanels = {
        makeMaterialPanel( Eigen::Vector3d::UnitX( ), 1.0, 300.0, 0.2, 0.3, 0.4, 0.5, busGroupId ),
        makeMaterialPanel( Eigen::Vector3d::UnitY( ), 1.0, 300.0, 0.2, 0.3, 0.4, 0.5, busGroupId )
    };
    std::vector< std::shared_ptr< VehicleExteriorPanel > > solarArrayPanels = { makeMaterialPanel(
            Eigen::Vector3d::UnitZ( ), 1.0, 300.0, 0.8, 0.7, 0.6, 0.5, solarArrayGroupId ) };

    SystemOfBodies bodies;
    bodies.createEmptyBody( "Vehicle" );
    std::shared_ptr< VehicleSystems > vehicleSystems = std::make_shared< VehicleSystems >( 100.0 );
    vehicleSystems->setVehicleExteriorPanels( { { busGroupId, busPanels }, { solarArrayGroupId, solarArrayPanels } } );
    bodies.at( "Vehicle" )->setVehicleSystems( vehicleSystems );

    std::shared_ptr< EstimatableParameter< double > > parameter =
            createDoubleParameterToEstimate< double, double >( energyAccommodationCoefficient( "Vehicle", busGroupId ), bodies );
    BOOST_CHECK( parameter->getParameterName( ).first == energy_accommodation_coefficient );
    BOOST_CHECK_EQUAL( parameter->getParameterName( ).second.first, "Vehicle" );
    BOOST_CHECK_EQUAL( parameter->getParameterName( ).second.second, busGroupId );

    parameter->setParameterValue( 0.45 );
    for( const auto& panel : busPanels )
    {
        BOOST_CHECK_CLOSE( panel->getEnergyAccomodationCoefficient( ), 0.45, 1.0E-12 );
    }
    BOOST_CHECK_CLOSE( solarArrayPanels.front( )->getEnergyAccomodationCoefficient( ), 0.8, 1.0E-12 );

    BOOST_CHECK_THROW(
            ( createDoubleParameterToEstimate< double, double >( energyAccommodationCoefficient( "Vehicle", "MissingGroup" ), bodies ) ),
            std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testMaterialOnlyPanelGroupIdIsPreserved )
{
    const std::string groupId = "Bus";
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )
            ->setRotationalEphemeris(
                    std::make_shared< SimpleRotationalEphemeris >( 0.0, 0.0, 0.0, 0.0, 0.0, "ECLIPJ2000", "VehicleFixed" ) );

    std::shared_ptr< BodyPanelSettings > panelSettings = bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitX( ), 1.0 ),
                                                                            nullptr,
                                                                            groupId,
                                                                            materialProperties( -1.0, -1.0, 0.8, 0.7, 0.6, 0.5 ) );
    addBodyExteriorPanelledShape( fullPanelledBodySettings( { panelSettings } ), "Vehicle", bodies );

    std::shared_ptr< EstimatableParameter< double > > parameter =
            createDoubleParameterToEstimate< double, double >( energyAccommodationCoefficient( "Vehicle", groupId ), bodies );
    BOOST_CHECK_EQUAL( parameter->getParameterName( ).second.second, groupId );
    BOOST_CHECK_CLOSE( parameter->getParameterValue( ), 0.8, 1.0E-12 );
}

BOOST_AUTO_TEST_CASE( testPanelMaterialPropertyParameterFactoryWithoutVehicleSystems )
{
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Vehicle" );
    BOOST_CHECK_THROW( ( createDoubleParameterToEstimate< double, double >( energyAccommodationCoefficient( "Vehicle", "Bus" ), bodies ) ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
