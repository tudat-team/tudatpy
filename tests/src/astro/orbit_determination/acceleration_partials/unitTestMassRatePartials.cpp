/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantDragCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/empiricalAccelerationCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantThrust.h"
#include "tudat/astro/orbit_determination/massDerivativePartial.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/relativity/relativisticAccelerationCorrection.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/estimation_setup/createStateDerivativePartials.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electromagnetism;
using namespace tudat::basic_astrodynamics;
using namespace tudat::propulsion;
using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_acceleration_partials )



BOOST_AUTO_TEST_CASE( testMassRatePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    for( unsigned int i = 0; i < 2; i++ )
    {
        // Create empty bodies, earth and sun.
        std::shared_ptr< Body > vehicle = std::make_shared< Body >( );

        SystemOfBodies bodies;
        bodies.addBody( vehicle, "Vehicle" );

        double vehicleMass = 5.0E3;
        bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

        Eigen::Vector3d thrustDirection;
        thrustDirection << -1.4, 2.4, 5.6;

        std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                [=](const double){ return thrustDirection; };
        bodies.at( "Vehicle" )->setRotationalEphemeris(
                    createRotationModel(
                        std::make_shared< BodyFixedDirectionBasedRotationSettings >(
                            thrustDirectionFunction, "ECLIPJ2000", "VehicleFixed" ),
                        "Vehicle", bodies ) );


        double thrustMagnitude1 = 1.0E3;
        double specificImpulse1 = 250.0;
        addEngineModel( "Vehicle", "Engine1",
                        std::make_shared< ConstantThrustMagnitudeSettings >(
                            thrustMagnitude1, specificImpulse1 ), bodies );

        double thrustMagnitude2 = 20.0E3;
        double specificImpulse2 = 50.0;
        addEngineModel( "Vehicle", "Engine2",
                        std::make_shared< ConstantThrustMagnitudeSettings >(
                            thrustMagnitude2, specificImpulse2 ), bodies, Eigen::Vector3d::UnitY( ) );

        // Create acceleration due to sun on earth.
        std::shared_ptr< ThrustAcceleration > thrustAcceleration;
        if( i == 0 )
        {
            thrustAcceleration = std::dynamic_pointer_cast< ThrustAcceleration >(
                        tudat::simulation_setup::createThrustAcceleratioModel(
                            std::make_shared< ThrustAccelerationSettings >( "Engine1" ), bodies, "Vehicle" ) );
        }
        else if( i == 1 )
        {
            thrustAcceleration = std::dynamic_pointer_cast< ThrustAcceleration >(
                        tudat::simulation_setup::createThrustAcceleratioModel(
                            std::make_shared< ThrustAccelerationSettings >(
                                std::vector< std::string >( { "Engine1", "Engine2" } ) ), bodies, "Vehicle" ) );
        }

        std::shared_ptr< FromThrustMassRateModel > massRateModel =  std::make_shared< FromThrustMassRateModel >(
                    thrustAcceleration );

        std::shared_ptr< EstimatableParameter< double > > constantThrustParameter1 = std::make_shared<
                ConstantThrustMagnitudeParameter >(
                    std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >(
                        vehicle->getVehicleSystems( )->getEngineModels( ).at( "Engine1" )->getThrustMagnitudeWrapper( ) ),
                    "Vehicle", "Engine1" );
        std::shared_ptr< EstimatableParameter< double > > constantThrustParameter2 = std::make_shared<
                ConstantThrustMagnitudeParameter >(
                    std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >(
                        vehicle->getVehicleSystems( )->getEngineModels( ).at( "Engine2" )->getThrustMagnitudeWrapper( ) ),
                    "Vehicle", "Engine2" );

        std::shared_ptr< EstimatableParameter< double > > constantSpecificImpulseParameter1 = std::make_shared<
                ConstantSpecificImpulseParameter< ConstantThrustMagnitudeWrapper > >(
                    std::dynamic_pointer_cast< ConstantThrustMagnitudeWrapper >(
                        vehicle->getVehicleSystems( )->getEngineModels( ).at( "Engine1" )->getThrustMagnitudeWrapper( ) ),
                    "Vehicle", "Engine1" );
        std::shared_ptr< EstimatableParameter< double > > constantSpecificImpulseParameter2 = std::make_shared<
                ConstantSpecificImpulseParameter< ConstantThrustMagnitudeWrapper > >(
                    std::dynamic_pointer_cast< ConstantThrustMagnitudeWrapper >(
                        vehicle->getVehicleSystems( )->getEngineModels( ).at( "Engine2" )->getThrustMagnitudeWrapper( ) ),
                    "Vehicle", "Engine2" );


        // Create central gravity partial.
        std::shared_ptr< FromThrustMassRatePartial > massRatePartial =
                std::dynamic_pointer_cast< FromThrustMassRatePartial >(
                    createAnalyticalMassRatePartial( massRateModel, std::make_pair( "Vehicle", vehicle ),
                                                     bodies ) );

        BOOST_CHECK_EQUAL( massRatePartial == nullptr, false );
        BOOST_CHECK_EQUAL( massRatePartial->isMassRatePartialWrtMassNonZero( ), false );


        // Calculate analytical partials.
        massRatePartial->update( 0.0 );


        Eigen::MatrixXd partialWrtMass = Eigen::Vector1d::Zero( );
        massRatePartial->wrtMassOfBody( partialWrtMass.block( 0, 0, 1, 1 ) );

        double partialWrtEngine1Thrust = massRatePartial->wrtParameter(
                    constantThrustParameter1 )( 0 );
        double partialWrtSpecificImpulse1 = massRatePartial->wrtParameter(
                    constantSpecificImpulseParameter1 )( 0 );
        double partialWrtEngine2Thrust = massRatePartial->wrtParameter(
                    constantThrustParameter2 )( 0 );
        double partialWrtSpecificImpulse2 = massRatePartial->wrtParameter(
                    constantSpecificImpulseParameter2 )( 0 );


        double testPartialWrtEngine1Thrust = calculateMassRateWrtParameterPartials(
                    constantThrustParameter1, massRateModel, 1.0 );
        double testPartialWrtSpecificImpulse1 = calculateMassRateWrtParameterPartials(
                    constantSpecificImpulseParameter1, massRateModel, 0.001 );
        double testPartialWrtEngine2Thrust = calculateMassRateWrtParameterPartials(
                    constantThrustParameter2, massRateModel, 1.0 );
        double testPartialWrtSpecificImpulse2 = calculateMassRateWrtParameterPartials(
                    constantSpecificImpulseParameter2, massRateModel, 0.001 );

        BOOST_CHECK_EQUAL( partialWrtMass( 0 ), 0.0 );
        BOOST_CHECK_CLOSE_FRACTION( testPartialWrtEngine1Thrust, partialWrtEngine1Thrust, 1.0E-10 );
        BOOST_CHECK_CLOSE_FRACTION( testPartialWrtSpecificImpulse1, partialWrtSpecificImpulse1, 1.0E-8 );
        if( i == 0 )
        {
            BOOST_CHECK_EQUAL( partialWrtEngine2Thrust, 0.0 );
            BOOST_CHECK_EQUAL( partialWrtSpecificImpulse2, 0.0 );
        }
        else if( i == 1 )
        {
            BOOST_CHECK_CLOSE_FRACTION( testPartialWrtEngine2Thrust, partialWrtEngine2Thrust, 1.0E-10 );
            BOOST_CHECK_CLOSE_FRACTION( testPartialWrtSpecificImpulse2, partialWrtSpecificImpulse2, 1.0E-8 );
        }
    }
}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
