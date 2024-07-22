/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Noomen, R. AE2230-I Flight and Orbital Mechanics Lecture Notes, Ch. Perturbations (2),
 *             Delft University of Technology, 2022.
 *      Knocke, Philip et al. "Earth radiation pressure effects on satellites."
 *          Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <memory>
#include <tuple>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/math/basic/rotationRepresentations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"


namespace tudat
{
namespace unit_tests
{
using namespace tudat;
using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiation_pressure_acceleration)

//! Test acceleration for idealized GOCE spacecraft (Noomen 2022)
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_GOCE )
{
    const double expectedReceivedIrradiance = 1371;
    const Eigen::Vector3d expectedAcceleration = Eigen::Vector3d(1, 1, 0).normalized() * 5.2e-9;
    const double expectedReceivedFraction = 1.0;

    auto luminosityModel = std::make_shared<ConstantLuminosityModel>(
        computeLuminosityFromIrradiance( 1371.0, physical_constants::ASTRONOMICAL_UNIT ) );
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1.2);
    IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
            sourceModel,
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
            [] () { return Eigen::Vector3d::Zero(); },
            targetModel,
            []() { return (Eigen::Vector3d(1, 1, 0).normalized() * physical_constants::ASTRONOMICAL_UNIT).eval(); },
            []() { return Eigen::Quaterniond::Identity(); },
            []() { return 1050; },
            std::make_shared<NoOccultingBodyOccultationModel>());

    sourceModel->updateMembers(TUDAT_NAN);
    targetModel->updateMembers(TUDAT_NAN);
    accelerationModel.updateMembers(TUDAT_NAN);

    const auto actualReceivedIrradiance = accelerationModel.getReceivedIrradiance();
    const auto actualAcceleration = accelerationModel.getAcceleration();
    const auto actualReceivedFraction = accelerationModel.getSourceToTargetReceivedFraction();

    BOOST_CHECK_CLOSE_FRACTION(actualReceivedIrradiance, expectedReceivedIrradiance, 1e-15);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-2);
    BOOST_CHECK_CLOSE_FRACTION(actualReceivedFraction, expectedReceivedFraction, 1e-15);
}

//! Test that acceleration of cannonball target with isotropic point source is invariant under rotation
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_IsotropicPointSource_CannonballTarget_RotationInvariance )
{
    using mathematical_constants::PI;

    std::vector<Eigen::Vector3d> actualAccelerations;

    // Iterate over arbitrary values for Euler angles
    for (double x : {0., 0.984, 2.579, 2*PI})
    {
        for (double y : {0., 0.743, 1.903, PI})
        {
            for (double z : {0., 0.646, 5.634, 2*PI})
            {
                auto rotation = basic_mathematics::getQuaternionFrom313EulerAngles(Eigen::Vector3d(x, y, z));

                auto luminosityModel = std::make_shared<ConstantLuminosityModel>(1);
                auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
                auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(42, 1.42);
                IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
                        sourceModel,
                        std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                        [] () { return Eigen::Vector3d::Zero(); },
                        targetModel,
                        [=] () { return Eigen::Vector3d(0, 1, 1); },
                        [=] () { return rotation; }, // vary target rotation
                        []() { return 1; },
                        std::make_shared<NoOccultingBodyOccultationModel>());

                sourceModel->updateMembers(TUDAT_NAN);
                targetModel->updateMembers(TUDAT_NAN);
                accelerationModel.updateMembers(TUDAT_NAN);

                actualAccelerations.push_back(accelerationModel.getAcceleration());
            }
        }
    }

    // Check that all calculated accelerations are identical
    BOOST_CHECK(std::all_of(
            actualAccelerations.begin(),
            actualAccelerations.end(),
            [&] (Eigen::Vector3d const &e) { return e.isApprox(actualAccelerations.front()); }));
}

//! Test occultation for cannonball target acceleration with isotropic point source
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_IsotropicPointSource_CannonballTarget_Occultation )
{
    {
        // Occulting body not interfering
        const double expectedReceivedFraction = 1.0;

        const auto occultingBodyPosition = Eigen::Vector3d(5, 5, 0);
        const auto targetPosition = Eigen::Vector3d(10, 0, 0);

        auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
            std::make_shared<ConstantLuminosityModel>(1));
        const auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return occultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [] () { return Eigen::Vector3d::Zero(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                sourceToTargetOccultationModel);

        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        sourceToTargetOccultationModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualReceivedFraction = accelerationModel.getSourceToTargetReceivedFraction();

        BOOST_CHECK(actualAcceleration > 0);
        BOOST_CHECK_CLOSE_FRACTION(actualReceivedFraction, expectedReceivedFraction, 1e-15);
    }

    {
        // Occulting body interfering with source -> target (target in umbra)
        const double expectedReceivedFraction = 0.0;

        const auto occultingBodyPosition = Eigen::Vector3d(5, 0, 0);
        const auto targetPosition = Eigen::Vector3d(10, 0, 0);

        auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
            std::make_shared<ConstantLuminosityModel>(1));
        const auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return occultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [] () { return Eigen::Vector3d::Zero(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                sourceToTargetOccultationModel);

        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        sourceToTargetOccultationModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualReceivedFraction = accelerationModel.getSourceToTargetReceivedFraction();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK_CLOSE_FRACTION(actualReceivedFraction, expectedReceivedFraction, 1e-15);
    }

    {
        // Occulting body interfering with source -> target (target in penumbra)
        const double expectedReceivedFraction = 0.4547;

        const auto occultingBodyPosition = Eigen::Vector3d::Zero();
        const auto occultingBodyRadius = 6378.137e3;

        const Eigen::Vector3d targetDirection( 0.018, 1.0, 0.0 );
        const Eigen::Vector3d targetPosition = (occultingBodyRadius + 1.0e3) * targetDirection.normalized();
        const Eigen::Vector3d sourcePosition = -149598000.0e3 * Eigen::Vector3d::UnitX();

        auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
            std::make_shared<ConstantLuminosityModel>(
                computeLuminosityFromIrradiance( 1367.0, physical_constants::ASTRONOMICAL_UNIT ) ));
        const auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(5, 1.2);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return occultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(occultingBodyRadius));

        IsotropicPointSourceRadiationPressureAcceleration unoccultedAccelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(6.96e8),
                [=] () { return sourcePosition; },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 42; },
                std::make_shared<NoOccultingBodyOccultationModel>());

        IsotropicPointSourceRadiationPressureAcceleration occultedAccelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(6.96e8),
                [=] () { return sourcePosition; },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 42; },
                sourceToTargetOccultationModel);

        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        sourceToTargetOccultationModel->updateMembers(TUDAT_NAN);
        unoccultedAccelerationModel.updateMembers(TUDAT_NAN);
        occultedAccelerationModel.updateMembers(TUDAT_NAN);

        const auto unoccultedAcceleration = unoccultedAccelerationModel.getAcceleration().norm();
        const auto occultedAcceleration = occultedAccelerationModel.getAcceleration().norm();

        BOOST_CHECK_CLOSE_FRACTION(occultedAcceleration / unoccultedAcceleration, expectedReceivedFraction, 1e-4);
    }
}

//! Test basic cases for paneled target acceleration with isotropic point source
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_IsotropicPointSource_PaneledTarget_Basic )
{
    // Set distance to speed of light to cancel to unity radiation pressure
    auto luminosityModel = std::make_shared<ConstantLuminosityModel>(
        computeLuminosityFromIrradiance( physical_constants::SPEED_OF_LIGHT, 1.0 ) );
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panels{
            std::make_shared< system_models::VehicleExteriorPanel >(1, -Eigen::Vector3d::UnitX(),
                        reflectionLawFromAbsorptivityAndDiffuseReflectivity(1, 0)),
            std::make_shared< system_models::VehicleExteriorPanel >(2, -Eigen::Vector3d::UnitY(),
                        reflectionLawFromSpecularAndDiffuseReflectivity(0, 1)),
            std::make_shared< system_models::VehicleExteriorPanel >(3, Eigen::Vector3d::UnitX(), // never pointing towards source in these tests
                        reflectionLawFromAbsorptivityAndDiffuseReflectivity(0.3, 0.4))
    };
    auto targetModel = std::make_shared<PaneledRadiationPressureTargetModel>(panels);

    // Using the same target position but different target rotations
    {
        // Only panel 1 towards source
        // Magnitude 1 because area 1 and only absorption
        const Eigen::Vector3d expectedAcceleration = Eigen::Vector3d::UnitX();
        IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [] () { return Eigen::Vector3d::Zero(); },
                targetModel,
                [] () { return Eigen::Vector3d::UnitX(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                std::make_shared<NoOccultingBodyOccultationModel>());

        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration();

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-15);
    }

    {
        // Only panel 2 towards source
        // Magnitude 3.333 because area 2 and factor 1.667 from diffuse reflection
        const Eigen::Vector3d expectedAcceleration = 2 * (1 + 2./3) * Eigen::Vector3d::UnitX();
        IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [] () { return Eigen::Vector3d::Zero(); },
                targetModel,
                [] () { return Eigen::Vector3d::UnitX(); },
                [] () {
                    const auto angle = unit_conversions::convertDegreesToRadians(-90.);
                    return Eigen::Quaterniond(Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()));
                },
                [] () { return 1; },
                std::make_shared<NoOccultingBodyOccultationModel>());

        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration();

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( actualAcceleration( i ) - expectedAcceleration( i ) ), ( 1.0E-15 * expectedAcceleration.norm( ) ) );
        }
    }

    {
        // Panel 1 and 2 angled 45° towards source
        // Panel 1 gives magnitude 1/√2 away from source (effective area is 1/√2)
        const Eigen::Vector3d expectedAccelerationDueToPanel1 = Eigen::Vector3d::UnitX() * 1 / sqrt(2);
        // Panel 2 due to diffuse reflection (effective area is 2/√2)
        const Eigen::Vector3d expectedAccelerationDueToPanel2 = (
                    1 * Eigen::Vector3d::UnitX() // incident light
                    + 2./3 * Eigen::Vector3d(1, 1, 0).normalized() // diffuse reflection
                ) * 2 / sqrt(2);
        const Eigen::Vector3d expectedAcceleration = expectedAccelerationDueToPanel1 + expectedAccelerationDueToPanel2;
        IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [] () { return Eigen::Vector3d::Zero(); },
                targetModel,
                [] () { return Eigen::Vector3d::UnitX(); },
                [] () {
                    const auto angle = unit_conversions::convertDegreesToRadians(-45.);
                    return Eigen::Quaterniond(Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()));
                },
                [] () { return 1; },
                std::make_shared<NoOccultingBodyOccultationModel>());

        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration();

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-15);
    }
}

//! Test radiation acceleration model for a paneled spacecraft in various orbits with respect to the Sun.
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_IsotropicPointSource_PaneledTarget_Realistic )
{
    // Box-and-wings model is partially obtained from Oliver Montenbruck, et al.
    //     "Semi-analytical solar radiation pressure modeling for QZS-1 orbit-normal and yaw-steering attitude".
    //     Advances in Space Research 59. 8(2017): 2088–2100.

    using namespace tudat::basic_astrodynamics;
    using namespace tudat::simulation_setup;
    using namespace tudat::ephemerides;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    SystemOfBodies bodies = createSystemOfBodies(
            getDefaultBodySettings( std::vector<std::string>{"Sun"}, initialEphemerisTime, finalEphemerisTime, "Sun" ) );
    auto radiationSourceModel =
            std::dynamic_pointer_cast<IsotropicPointRadiationSourceModel>(bodies.at("Sun")->getRadiationSourceModel());

    // Create vehicle
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 2000.0 );

    for ( int testCase = 0 ; testCase < 4 ; testCase++)
    {
        double inclination{};
        if ( testCase == 0 || testCase == 3 )
        {
            // Put vehicle on circular orbit around the Sun with i = 0 deg
            inclination = 0.0;
        }
        else if ( testCase == 1 )
        {
            // Put vehicle in circular orbit around the Sun with i = 90 deg
            inclination = 90.0;
        }
        else if ( testCase == 2 )
        {
            // Put vehicle in circular orbit around the Sun with arbitrary chosen inclination
            inclination = 20.0;
        }

        Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
        initialStateInKeplerianElements[ orbital_element_conversions::semiMajorAxisIndex ] = physical_constants::ASTRONOMICAL_UNIT;
        initialStateInKeplerianElements[ orbital_element_conversions::inclinationIndex ] = unit_conversions::convertDegreesToRadians( inclination );
        bodies.at( "Vehicle" )->setEphemeris(
                std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                     spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000") );

        auto orbitalPeriod = computeKeplerOrbitalPeriod(
                physical_constants::ASTRONOMICAL_UNIT,
                spice_interface::getBodyGravitationalParameter( "Sun" ));

        // Set-up rotational ephemeris for vehicle.
        if ( testCase < 3 ){
            // Define constant rotational model.
            Eigen::Vector7d rotationalStateVehicle;
            rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(Eigen::Quaterniond::Identity());
            bodies.at( "Vehicle" )->setRotationalEphemeris(
                    std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000", "VehicleFixed" ));
        }
        else if ( testCase == 3 )
        {
            // Define simple rotational model.
            bodies.at( "Vehicle" )->setRotationalEphemeris(
                    std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                            0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
        }
        bodies.processBodyFrameDefinitions( );


        std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panels;
        if ( testCase == 0 || testCase == 1 || testCase == 2 )
        {
            // Define panels properties for test cases with constant rotational model.
            panels = {
                    std::make_shared< system_models::VehicleExteriorPanel >(2.0, Eigen::Vector3d::UnitZ(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(4.0, Eigen::Vector3d::UnitZ(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    std::make_shared< system_models::VehicleExteriorPanel >(6.0, -Eigen::Vector3d::UnitZ(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(9.9, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(2.3, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    std::make_shared< system_models::VehicleExteriorPanel >(9.9, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(2.3, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    std::make_shared< system_models::VehicleExteriorPanel >(4.6, Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(2.7, Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    std::make_shared< system_models::VehicleExteriorPanel >(5.8, -Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(2.7, -Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
            };
        }
        else if ( testCase == 3 )
        {
            // Define panels properties for test cases with constant rotational model (simpler boxes-and-wings model).
            panels = {
                    std::make_shared< system_models::VehicleExteriorPanel >(9.9, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(2.3, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    std::make_shared< system_models::VehicleExteriorPanel >(9.9, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(2.3, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
            };
        }
        bodies.at( "Vehicle" )->setRadiationPressureTargetModels(
            { std::make_shared<PaneledRadiationPressureTargetModel>(panels) } );

        std::vector< double > areas;
        std::vector< Eigen::Vector3d > panelSurfaceNormals;
        std::vector< double > specularReflectivities;
        std::vector< double > diffuseReflectivities;
        for (auto& panel : panels)
        {
            areas.push_back(panel->getPanelArea());
            panelSurfaceNormals.push_back(panel->getFrameFixedSurfaceNormal( )( ) );
            auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(panel->getReflectionLaw());
            specularReflectivities.push_back(reflectionLaw->getSpecularReflectivity());
            diffuseReflectivities.push_back(reflectionLaw->getDiffuseReflectivity());
        }


        SelectedAccelerationMap accelerationMap {
                {"Vehicle", {
                        {"Sun", {
                                radiationPressureAcceleration()
                        }},
                }}
        };
        auto accelerationModelMap = createAccelerationModelsMap(bodies,
                                                                accelerationMap,
                                                                std::vector<std::string>{"Vehicle"},
                                                                std::vector<std::string>{"Sun"});
        auto accelerationModel = accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        std::vector< double > testTimes{ 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 / 4.0 * orbitalPeriod };


        // Compute panelled radiation pressure for various relative Sun positions.
        Eigen::Vector3d calculatedAcceleration, expectedAcceleration;

        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodies.at( "Sun" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setStateFromEphemeris( testTimes[ i ] );
            // Round vehicle state such that position vector only has a single non-zero component
            // This should physically be the case for the given test times, but the Kepler ephemeris returns small
            // non-zero values for other position components as well, leading to discrepancies between the simplified
            // panel accelerations calculated here (only considering the two pointing towards the sun) and
            // those calculated considering all panels in the tested class
            bodies.at( "Vehicle" )->setState(bodies.at( "Vehicle" )->getState().array().round());
            bodies.at( "Vehicle" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );

            bodies.at( "Sun" )->getRadiationSourceModel()->updateMembers( testTimes[ i ] );
            bodies.at( "Vehicle" )->getRadiationPressureTargetModel()->updateMembers( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAcceleration = accelerationModel->getAcceleration( );

            Eigen::Vector3d expectedVehicleToSunVector =
                    bodies.at( "Sun" )->getPosition() - bodies.at( "Vehicle" )->getPosition();
            Eigen::Vector3d expectedVehicleToSunVectorNormalized = expectedVehicleToSunVector.normalized();

            auto sourceIrradiance =
                    radiationSourceModel->evaluateIrradianceAtPosition(bodies.at( "Vehicle" )->getPosition()).front().first;
            auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;

            // Calculated accelerations.

            // Calculated acceleration if Sun-Vehicle vector aligned with +X-axis (acceleration generated by panels whose normals are
            // along +X-axis only)
            Eigen::Vector3d accelerationPositiveXaxisOrientedPanels, accelerationNegativeXaxisOrientedPanels, accelerationPositiveYaxisOrientedPanels, accelerationNegativeYaxisOrientedPanels;
            if( testCase < 3 )
            {
                double cosPanelInclinationPositiveXaxis = expectedVehicleToSunVectorNormalized.dot(Eigen::Vector3d::UnitX( ) );
                accelerationPositiveXaxisOrientedPanels = - cosPanelInclinationPositiveXaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[5] * (( 1.0 - specularReflectivities[5] ) * expectedVehicleToSunVectorNormalized
                        + 2.0 / 3.0 * diffuseReflectivities[5] * Eigen::Vector3d::UnitX( ) )
                        + areas[6] * (( 1.0 - specularReflectivities[6] ) * expectedVehicleToSunVectorNormalized
                        + ( 2.0 / 3.0 * diffuseReflectivities[6] + 2.0 * cosPanelInclinationPositiveXaxis * specularReflectivities[6] )
                        * Eigen::Vector3d::UnitX( ) ) );


                // Calculated acceleration generated by the panels whose normals are along -X-axis.
                double cosPanelInclinationNegativeXaxis = expectedVehicleToSunVectorNormalized.dot(- Eigen::Vector3d::UnitX() );
                accelerationNegativeXaxisOrientedPanels = - cosPanelInclinationNegativeXaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[3] * (( 1.0 - specularReflectivities[3] ) * expectedVehicleToSunVectorNormalized
                        + 2.0 / 3.0 * diffuseReflectivities[3] * - Eigen::Vector3d::UnitX( ) )
                        + areas[4] * (( 1.0 - specularReflectivities[4] ) * expectedVehicleToSunVectorNormalized
                        + ( 2.0 / 3.0 * diffuseReflectivities[4] + 2.0 * specularReflectivities[4] * cosPanelInclinationNegativeXaxis )
                        * - Eigen::Vector3d::UnitX( ) ) );

                // Calculated acceleration generated by the panels whose normals are along +Y-axis.
                double cosPanelInclinationPositiveYaxis = expectedVehicleToSunVectorNormalized.dot(Eigen::Vector3d::UnitY() );
                accelerationPositiveYaxisOrientedPanels = - cosPanelInclinationPositiveYaxis
                        * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                        * ( areas[7] * (( 1.0 - specularReflectivities[7] ) * expectedVehicleToSunVectorNormalized
                            + 2.0 / 3.0 * diffuseReflectivities[7] * Eigen::Vector3d::UnitY( ) )
                            + areas[8] * (( 1.0 - specularReflectivities[8] ) * expectedVehicleToSunVectorNormalized
                            + ( 2.0 / 3.0 * diffuseReflectivities[8] + 2.0 * specularReflectivities[8] * cosPanelInclinationPositiveYaxis )
                            * Eigen::Vector3d::UnitY( ) ) );

                // Calculated acceleration generated by the panels whose normals are along -Y-axis.
                double cosPanelInclinationNegativeYaxis = expectedVehicleToSunVectorNormalized.dot(- Eigen::Vector3d::UnitY() );
                accelerationNegativeYaxisOrientedPanels = - cosPanelInclinationNegativeYaxis
                        * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                        * ( areas[9] * (( 1.0 - specularReflectivities[9] ) * expectedVehicleToSunVectorNormalized
                            + 2.0 / 3.0 * diffuseReflectivities[9] * - Eigen::Vector3d::UnitY( ) )
                            + areas[10] * (( 1.0 - specularReflectivities[10] ) * expectedVehicleToSunVectorNormalized
                            + ( 2.0 / 3.0 * diffuseReflectivities[10] + 2.0 * specularReflectivities[10] * cosPanelInclinationNegativeYaxis )
                            * - Eigen::Vector3d::UnitY( ) ) );
            }

            // Calculated acceleration generated by the panels whose normals are along +Z-axis.
            double cosPanelInclinationPositiveZaxis = expectedVehicleToSunVectorNormalized.dot(Eigen::Vector3d::UnitZ() );
            Eigen::Vector3d accelerationPositiveZaxisOrientedPanels = - cosPanelInclinationPositiveZaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[0] * (( 1.0 - specularReflectivities[0] ) * expectedVehicleToSunVectorNormalized
                        + 2.0 / 3.0 * diffuseReflectivities[0] * Eigen::Vector3d::UnitZ( ) )
                        + areas[1] * (( 1.0 - specularReflectivities[1] ) * expectedVehicleToSunVectorNormalized
                        + ( 2.0 / 3.0 * diffuseReflectivities[1] + 2.0 * cosPanelInclinationPositiveZaxis * specularReflectivities[1] )
                        * Eigen::Vector3d::UnitZ( ) ) );

            // Calculated acceleration generated by the panels whose normals are along -Z-axis.
            double cosPanelInclinationNegativeZaxis = expectedVehicleToSunVectorNormalized.dot(- Eigen::Vector3d::UnitZ() );
            Eigen::Vector3d accelerationNegativeZaxisOrientedPanels = - cosPanelInclinationNegativeZaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * areas[2] * (( 1.0 - specularReflectivities[2] ) * expectedVehicleToSunVectorNormalized
                        + 2.0 / 3.0 * diffuseReflectivities[2] * - Eigen::Vector3d::UnitZ( ) );




            // Equatorial orbit case: vehicle in circular orbit around the Sun with i = 0 deg.
            if ( testCase == 0 ){

                // At t = 0: Sun-Vehicle vector expected along +X-axis
                // The acceleration is expected to be generated by the panels whose normals are along -X-axis only.
                if ( i == 0 ){
                    expectedAcceleration = accelerationNegativeXaxisOrientedPanels;
                }

                    // At t = period/4: Sun-Vehicle vector expected along +Y-axis
                    // The acceleration is expected to be generated by the panels whose normals are along -Y-axis only.
                else if ( i == 1 ){
                    expectedAcceleration = accelerationNegativeYaxisOrientedPanels;
                }

                    // At t = period/2: Sun-Vehicle vector expected along -X-axis
                    // The acceleration is expected to be generated by the panels whose normals are along +X-axis only.
                else if ( i == 2 ){
                    expectedAcceleration = accelerationPositiveXaxisOrientedPanels;
                }

                    // At t = 3 * period/4: Sun-Vehicle vector expected along -Y-axis
                    // The acceleration is expected to be generated by the panels whose normals are along +Y-axis only.
                else if ( i == 3 ){
                    expectedAcceleration = accelerationPositiveYaxisOrientedPanels;
                }
            }


            // Polar orbit case: vehicle in circular orbit around the Sun with i = 90 deg
            if ( testCase == 1 ){

                // At t = 0: Sun-Vehicle vector expected along +X-axis.
                // The acceleration is expected to be generated by the panels whose normals are along -X-axis only.
                if ( i == 0 ){
                    expectedAcceleration = accelerationNegativeXaxisOrientedPanels;
                }

                    // At t = period/4: Sun-Vehicle vector expected along +Z-axis.
                    // The acceleration is expected to be generated by the panels whose normals are along -Z-axis only.
                else if ( i == 1 ){
                    expectedAcceleration = accelerationNegativeZaxisOrientedPanels;
                }

                    // At t = period/2: Sun-Vehicle vector expected along -X-axis.
                    // The acceleration is expected to be generated by the panels whose normals are along +X-axis only.
                else if ( i == 2 ){
                    expectedAcceleration = accelerationPositiveXaxisOrientedPanels;
                }

                    // At t = 3*period/4: Sun-Vehicle vector expected along -Z-axis.
                    // The acceleration is expected to be generated by the panels whose normals are along +Z-axis only.
                else if ( i == 3 ){
                    expectedAcceleration = accelerationPositiveZaxisOrientedPanels;
                }
            }


            // Circular orbit with arbitrary chosen inclination
            if ( testCase == 2 ){

                // At t = 0: Sun-Vehicle vector expected along +X-axis.
                // The acceleration is expected to be generated by the panels whose normals are along -X-axis only.
                if ( i == 0 ){
                    expectedAcceleration = accelerationNegativeXaxisOrientedPanels;
                }

                    // At t = period/4: Sun-vehicle vector expected to have components along +Z and +Y axes.
                    // The acceleration is expected to be generated by the panels whose normals are along -Y and -Z axes.
                else if ( i == 1 ){
                    expectedAcceleration = accelerationNegativeZaxisOrientedPanels + accelerationNegativeYaxisOrientedPanels;
                }

                    // At t = period/2: Sun-Vehicle vector expected along -X-axis.
                    // The acceleration is expected to be generated by the panels whose normals are along +X-axis only.
                else if ( i == 2 ){
                    expectedAcceleration = accelerationPositiveXaxisOrientedPanels;
                }

                    // At t = period/4: Sun-vehicle vector expected to have components along -Z and -Y axes.
                    // The acceleration is expected to be generated by the panels whose normals are along +Y and +Z axes.
                else if ( i == 3 ){
                    expectedAcceleration = accelerationPositiveZaxisOrientedPanels + accelerationPositiveYaxisOrientedPanels;
                }

            }


            // Case with non-constant rotational model for the spacecraft.
            if ( testCase == 3 )
            {
                Eigen::Quaterniond currentRotationToInertialFrame =
                        bodies.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTimes.at( i ) );
                Eigen::Vector3d expectedPanelNormalPositiveXaxis = currentRotationToInertialFrame * Eigen::Vector3d::UnitX();
                Eigen::Vector3d expectedPanelNormalNegativeXaxis = currentRotationToInertialFrame * -Eigen::Vector3d::UnitX();


                double cosPanelInclinationPositiveXaxis = expectedVehicleToSunVectorNormalized.dot(expectedPanelNormalPositiveXaxis );

                // Determine the normal of the panels that are actually generating a radiation pressure acceleration,
                // depending on their current inertial orientation.
                Eigen::Vector3d expectedPanelNormal;
                if ( cosPanelInclinationPositiveXaxis >= 0.0 ){
                    expectedPanelNormal = expectedPanelNormalPositiveXaxis;
                }
                else{
                    expectedPanelNormal = expectedPanelNormalNegativeXaxis;
                }

                // Calculate expected acceleration (same characteristics for the panels oriented along the -X and +X axes).
                expectedAcceleration = -fabs( cosPanelInclinationPositiveXaxis )
                        * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                        * ( areas[0] * (( 1.0 - specularReflectivities[0] ) * expectedVehicleToSunVectorNormalized
                                + 2.0 / 3.0 * diffuseReflectivities[0] * expectedPanelNormal )
                        + areas[1] * (( 1.0 - specularReflectivities[1] ) * expectedVehicleToSunVectorNormalized
                                + ( 2.0 / 3.0 * diffuseReflectivities[1] + 2.0 * fabs(cosPanelInclinationPositiveXaxis) * specularReflectivities[1] )
                                * expectedPanelNormal ) );
            }

            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration[j] - expectedAcceleration[j] ), 3.0e-23 );
            }
        }
    }
}

//! Test radiation acceleration model for a paneled spacecraft with time-varying solar panels orientation.
// This is done by checking consistency of (i) a rotating spacecraft with fixed panels and
// (ii) a constant-attitude spacecraft with rotating panels
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_IsotropicPointSource_PaneledTarget_TimeVaryingPanelOrientation )
{
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::simulation_setup;
    using namespace tudat::ephemerides;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2000.0;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );


    // Compute spacecraft orbital period, and compute test times
    double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                                                                         spice_interface::getBodyGravitationalParameter( "Sun" ) );
    std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };

    // Put vehicle on circular orbit around Sun.
    Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                 spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );

    // Define rotational model parameters.
    std::vector< double > rightAscensionPole;
    rightAscensionPole.push_back( 0.0 );
    rightAscensionPole.push_back( 0.2 );

    std::vector< double > declinationPole;
    declinationPole.push_back( mathematical_constants::PI / 2.0 );
    declinationPole.push_back( 0.4 );

    std::vector< double > primeMeridianLongitude;
    primeMeridianLongitude.push_back( - mathematical_constants::PI / 2.0 );
    primeMeridianLongitude.push_back( - 0.2 );

    std::vector< double > rotationalRate;
    rotationalRate.push_back( 1.0E-5 );
    rotationalRate.push_back( 1.0E-5 );

    std::vector< double > numberSecondsSinceEpoch;
    numberSecondsSinceEpoch.push_back( 0.0 );
    numberSecondsSinceEpoch.push_back( 0.0 );


    // Case 0: vehicle-fixed axes aligned with inertial ones at time t = 0.
    // Case 1: arbitrary chosen rotational model for the vehicle.

    for ( unsigned int testCase = 0 ; testCase < 2 ; testCase++)
    {
        /// First calculation with simple rotational ephemeris and constant panel orientation.
        std::vector<Eigen::Vector3d> calculatedAcceleration;
        {
            // Define simple rotational ephemeris.
            bodies.at("Vehicle")->setRotationalEphemeris(std::make_shared<tudat::ephemerides::SimpleRotationalEphemeris>(
                    rightAscensionPole[testCase], declinationPole[testCase], primeMeridianLongitude[testCase],
                    rotationalRate[testCase], numberSecondsSinceEpoch[testCase], "ECLIPJ2000", "VehicleFixed"));

            // Create panels
            std::vector<std::shared_ptr< system_models::VehicleExteriorPanel >> panels{
                    std::make_shared< system_models::VehicleExteriorPanel >(2.0, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >(4.0, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
            };
            bodies.at("Vehicle")->setRadiationPressureTargetModels(
                { std::make_shared<PaneledRadiationPressureTargetModel>(panels) } );

            SelectedAccelerationMap accelerationMap{
                    {"Vehicle", {
                            {"Sun", {
                                    radiationPressureAcceleration( paneled_target )
                            }},
                    }}
            };
            auto accelerationModelMap = createAccelerationModelsMap(bodies,
                                                                    accelerationMap,
                                                                    std::vector<std::string>{"Vehicle"},
                                                                    std::vector<std::string>{"Sun"});
            auto accelerationModel = accelerationModelMap.at("Vehicle").at("Sun").at(0);


            // Compute radiation pressure acceleration for different Sun positions.
            Eigen::Vector3d sunCenteredVehiclePosition;
            std::shared_ptr<Ephemeris> vehicleEphemeris = bodies.at("Vehicle")->getEphemeris();

            for (unsigned int i = 0; i < testTimes.size(); i++)
            {
                auto currentTime = testTimes[i];

                // Update environment and acceleration
                bodies.at("Sun")->setStateFromEphemeris(currentTime);
                bodies.at("Vehicle")->setStateFromEphemeris(currentTime);
                bodies.at("Vehicle")->setCurrentRotationToLocalFrameFromEphemeris(currentTime);
                bodies.at("Vehicle")->updateMass(currentTime);
                bodies.at("Sun")->getRadiationSourceModel()->updateMembers(currentTime);
                simulation_setup::getRadiationPressureTargetModelOfType(
                    bodies.at("Vehicle"), simulation_setup::paneled_target )->updateMembers(currentTime);
                accelerationModel->updateMembers(currentTime);

                // Retrieve acceleration.
                calculatedAcceleration.push_back(accelerationModel->getAcceleration());
            }
        }



        /// Second calculation with constant rotational ephemeris and time-varying panel orientation
        std::vector< Eigen::Vector3d > calculatedAccelerationTimeVaryingPanelOrientation;
        {
            // Define constant rotational ephemeris
            Eigen::Vector7d rotationalStateVehicle;
            rotationalStateVehicle.segment(0, 4) = linear_algebra::convertQuaternionToVectorFormat(
                    Eigen::Quaterniond(Eigen::Matrix3d::Identity()));
            rotationalStateVehicle.segment(4, 3) = Eigen::Vector3d::Zero();
            bodies.at("Vehicle")->setRotationalEphemeris(
                    std::make_shared<ConstantRotationalEphemeris>(rotationalStateVehicle, "ECLIPJ2000",
                                                                  "VehicleFixed"));

            // Define time-varying panel orientation.
            double currentTime;

            const auto rotationToInertialFrameFunction = [&]()
            {
                Eigen::Quaterniond rotationAngle =
                        reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                                basic_mathematics::computeModulo(
                                        (currentTime - numberSecondsSinceEpoch[testCase]) * rotationalRate[testCase],
                                        2.0 * mathematical_constants::PI));
                Eigen::Quaterniond rotationAxis =
                        reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                                declinationPole[testCase],
                                rightAscensionPole[testCase],
                                primeMeridianLongitude[testCase]);

                return (rotationAngle * rotationAxis).inverse();
            };

            // Create panels
            std::vector<std::shared_ptr< system_models::VehicleExteriorPanel >> panels{
                    std::make_shared< system_models::VehicleExteriorPanel >([=] () { return rotationToInertialFrameFunction() * Eigen::Vector3d::UnitX(); }, 2.0, "",
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    std::make_shared< system_models::VehicleExteriorPanel >([=] () { return rotationToInertialFrameFunction() * -Eigen::Vector3d::UnitX(); }, 4.0, "",
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
            };
            bodies.at("Vehicle")->setRadiationPressureTargetModels(
                { std::make_shared<PaneledRadiationPressureTargetModel>(panels), std::make_shared< CannonballRadiationPressureTargetModel >( 1000.0, 0.3 ) } );

            SelectedAccelerationMap accelerationMap{
                    {"Vehicle", {
                            {"Sun", {
                                    radiationPressureAcceleration( paneled_target )
                            }},
                    }}
            };
            auto accelerationModelMap = createAccelerationModelsMap(bodies,
                                                                    accelerationMap,
                                                                    std::vector<std::string>{"Vehicle"},
                                                                    std::vector<std::string>{"Sun"});
            auto accelerationModelTimeVaryingPanelSurfaceNormal = accelerationModelMap.at("Vehicle").at("Sun").at(0);


            // Compute radiation pressure acceleration for different Sun positions.
            Eigen::Vector3d sunCenteredVehiclePosition;
            std::shared_ptr<Ephemeris> vehicleEphemeris = bodies.at("Vehicle")->getEphemeris();

            for (unsigned int i = 0; i < testTimes.size(); i++)
            {
                currentTime = testTimes[i];

                // Update environment and acceleration
                bodies.at("Sun")->setStateFromEphemeris(currentTime);
                bodies.at("Vehicle")->setStateFromEphemeris(currentTime);
                bodies.at("Vehicle")->setCurrentRotationToLocalFrameFromEphemeris(currentTime);
                bodies.at("Vehicle")->updateMass(currentTime);
                bodies.at("Sun")->getRadiationSourceModel()->updateMembers(currentTime);
                simulation_setup::getRadiationPressureTargetModelOfType(
                    bodies.at("Vehicle"), simulation_setup::paneled_target )->updateMembers(currentTime);
                accelerationModelTimeVaryingPanelSurfaceNormal->updateMembers(currentTime);

                // Retrieve acceleration.
                calculatedAccelerationTimeVaryingPanelOrientation.push_back(
                        accelerationModelTimeVaryingPanelSurfaceNormal->getAcceleration());
            }
        }

        for( unsigned int j = 0; j < testTimes.size() ; j++ )
        {
            for ( unsigned int i = 0 ; i < 3 ; i++ ){
                BOOST_CHECK_SMALL( std::fabs(
                                       calculatedAcceleration[j][i] - calculatedAccelerationTimeVaryingPanelOrientation[j][i] ), 3.0e-23 );
            }
        }
    }
}


//! Test basic case for paneled target acceleration with paneled source
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_StaticallyPaneledSource_PaneledTarget_Basic )
{
    // Set distance to speed of light to cancel to unity radiation pressure
    auto luminosityModel = std::make_shared<ConstantLuminosityModel>(
        computeLuminosityFromIrradiance( physical_constants::SPEED_OF_LIGHT, 1.0 ) );
    auto originalSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);

    const std::map<std::string, std::shared_ptr<IsotropicPointRadiationSourceModel>>& originalSourceModels {
        {"OrigSource", originalSourceModel}};
    const std::map<std::string, std::shared_ptr<basic_astrodynamics::BodyShapeModel>>& originalSourceBodyShapeModels {
        {"OrigSource", std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1)}};
    const std::map<std::string, std::function<Eigen::Vector3d()>>& originalSourcePositionFunctions {
        {"OrigSource", [] { return Eigen::Vector3d::UnitX(); }}};
    const std::map<std::string, std::shared_ptr<OccultationModel>>& originalSourceToSourceOccultationModels {
        {"OrigSource", std::make_shared<NoOccultingBodyOccultationModel>()}};
    auto sourcePanelRadiosityModelUpdater = std::make_unique<SourcePanelRadiosityModelUpdater>(
            [] { return Eigen::Vector3d::Zero(); },
            [] { return Eigen::Quaterniond::Identity(); },
            originalSourceModels, originalSourceBodyShapeModels, originalSourcePositionFunctions, originalSourceToSourceOccultationModels);

    std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
    radiosityModels.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
            "OrigSource", std::make_shared<ConstantSurfacePropertyDistribution>(1)));

    // Source is a single panel pointing in +X with purely diffuse reflection
    std::vector<RadiationSourcePanel> panels;
    panels.emplace_back(
            2, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitX(), std::move(radiosityModels));
    auto sourceModel = std::make_shared<StaticallyPaneledRadiationSourceModel>(
            std::move(sourcePanelRadiosityModelUpdater), std::move(panels));

    std::vector<std::shared_ptr< system_models::VehicleExteriorPanel >> targetPanels{
            std::make_shared< system_models::VehicleExteriorPanel >(1, -Eigen::Vector3d::UnitX(),
                        reflectionLawFromAbsorptivityAndDiffuseReflectivity(1, 0)),
            std::make_shared< system_models::VehicleExteriorPanel >(3, Eigen::Vector3d::UnitX(), // never pointing towards source in these tests
                        reflectionLawFromAbsorptivityAndDiffuseReflectivity(0.3, 0.4))
    };
    auto targetModel = std::make_shared<PaneledRadiationPressureTargetModel>(targetPanels);

    // Only panel 1 towards source
    // Magnitude 2/pi because target area 1, source area 2, for target only absorption, pi due to diffuse albedo of source
    const Eigen::Vector3d expectedAcceleration = 2 * Eigen::Vector3d::UnitX() / mathematical_constants::PI;
    PaneledSourceRadiationPressureAcceleration accelerationModel(
            sourceModel,
            [] () { return Eigen::Vector3d::Zero(); },
            [] () { return Eigen::Quaterniond::Identity(); },
            targetModel,
            [] () { return Eigen::Vector3d::UnitX(); },
            [] () { return Eigen::Quaterniond::Identity(); },
            [] () { return 1; },
            std::make_shared<NoOccultingBodyOccultationModel>());

    originalSourceModel->updateMembers(TUDAT_NAN);
    sourceModel->updateMembers(TUDAT_NAN);
    targetModel->updateMembers(TUDAT_NAN);
    accelerationModel.updateMembers(TUDAT_NAN);

    const auto actualAcceleration = accelerationModel.getAcceleration();

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-15);
}

//! Test occultation for cannonball target acceleration with paneled source
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_StaticallyPaneledSource_CannonballTarget_Occultation )
{
    const auto originalSourcePosition = Eigen::Vector3d(10, 0, 0);
    auto originalSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
        std::make_shared<ConstantLuminosityModel>(1));

    Eigen::Vector3d originalSourceToSourceOccultingBodyPosition;

    std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
    radiosityModels.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
            "OrigSource", std::make_shared<ConstantSurfacePropertyDistribution>(1)));

    const std::map<std::string, std::shared_ptr<IsotropicPointRadiationSourceModel>>& originalSourceModels {
        {"OrigSource", originalSourceModel}};
    const std::map<std::string, std::shared_ptr<basic_astrodynamics::BodyShapeModel>>& originalSourceBodyShapeModels {
        {"OrigSource", std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1)}};
    const std::map<std::string, std::function<Eigen::Vector3d()>>& originalSourcePositionFunctions {
        {"OrigSource", [=] { return originalSourcePosition; }}};
    const std::map<std::string, std::shared_ptr<OccultationModel>>& originalSourceToSourceOccultationModels {
        {"OrigSource", std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [&] () { return originalSourceToSourceOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1))}};
    auto sourcePanelRadiosityModelUpdater = std::make_unique<SourcePanelRadiosityModelUpdater>(
            [] { return Eigen::Vector3d::Zero(); },
            [] { return Eigen::Quaterniond::Identity(); },
            originalSourceModels, originalSourceBodyShapeModels, originalSourcePositionFunctions, originalSourceToSourceOccultationModels);

    // Source is a single panel pointing in +X
    std::vector<RadiationSourcePanel> panels;
    panels.emplace_back(
            1, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitX(), std::move(radiosityModels));
    auto sourceModel = std::make_shared<StaticallyPaneledRadiationSourceModel>(
            std::move(sourcePanelRadiosityModelUpdater), std::move(panels));

    const auto targetPosition = Eigen::Vector3d(10, 10, 0);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1);

    // Original source, source and target keep their position, occulting bodies move between tests
    {
        // Two occulting bodies but not interfering with radiation
        const unsigned int expectedVisibleAndEmittingSourcePanelCount = 1;

        originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(-5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, -5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                sourceToTargetOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleAndEmittingSourcePanelCount = accelerationModel.getVisibleAndEmittingSourcePanelCount();

        BOOST_CHECK(actualAcceleration > 0);
        BOOST_CHECK(actualVisibleAndEmittingSourcePanelCount == expectedVisibleAndEmittingSourcePanelCount);
    }

    {
        // Occulting body only interfering with original source -> source
        const unsigned int expectedVisibleAndEmittingSourcePanelCount = 0;

        originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, -5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                sourceToTargetOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleAndEmittingSourcePanelCount = accelerationModel.getVisibleAndEmittingSourcePanelCount();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK(actualVisibleAndEmittingSourcePanelCount == expectedVisibleAndEmittingSourcePanelCount);
    }

    {
        // Occulting body only interfering with source -> target
        const unsigned int expectedVisibleAndEmittingSourcePanelCount = 0;

        originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(-5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, 5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                sourceToTargetOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleAndEmittingSourcePanelCount = accelerationModel.getVisibleAndEmittingSourcePanelCount();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK(actualVisibleAndEmittingSourcePanelCount == expectedVisibleAndEmittingSourcePanelCount);
    }

    {
        // Occulting bodies interfering with both original source -> source and source -> target
        const unsigned int expectedVisibleAndEmittingSourcePanelCount = 0;

        originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, 5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                sourceToTargetOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleAndEmittingSourcePanelCount = accelerationModel.getVisibleAndEmittingSourcePanelCount();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK(actualVisibleAndEmittingSourcePanelCount == expectedVisibleAndEmittingSourcePanelCount);
    }
}

//! Test radiation acceleration model for LAGEOS with albedo and thermal radiation from Earth (Knocke 1988)
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_DynamicallyPanelunitestradedSource_CannonballTarget_LAGEOS )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::electromagnetism;
    using namespace tudat::ephemerides;
    using namespace tudat::basic_astrodynamics;

    spice_interface::loadStandardSpiceKernels( );

    double area = 0.28483;
    double coefficient = 1.13;
    double bodyMass = 406.9;

    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "J2000";
    const auto bodySettings = getDefaultBodySettings({"Earth", "Sun"}, globalFrameOrigin, globalFrameOrientation);
    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    const auto targetModelSettings = cannonballRadiationPressureTargetModelSettings(area, coefficient);
    const auto targetModel = createRadiationPressureTargetModel(targetModelSettings, "Vehicle", bodies).at( 0 );

    const auto sourceModelSettings = extendedRadiationSourceModelSettings({
            albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun"),
            delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun")
    }, {6, 12});
    const auto sourceModel =
            std::dynamic_pointer_cast<PaneledRadiationSourceModel>(createRadiationSourceModel(sourceModelSettings, globalFrameOrigin, bodies));

    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    auto ephemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    const auto startTime = spice_interface::convertDateStringToEphemerisTime("1977-02-22");

    Eigen::Vector6d initialStateInKeplerianElements =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                   ephemeris->getCartesianState(startTime), celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER);

    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
        initialStateInKeplerianElements[ orbital_element_conversions::semiMajorAxisIndex ],
        celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER, 0 );

    const int n_rev = 3;
    const double n_steps_per_rev = 10;
    for (int i = 0; i <= n_rev * n_steps_per_rev; ++i)
    {
        double t = startTime + i / n_steps_per_rev * orbitalPeriod;

        bodies.at("Earth")->setStateFromEphemeris(t);
        bodies.at("Earth")->setCurrentRotationToLocalFrameFromEphemeris(t);
        bodies.at("Sun")->setStateFromEphemeris(t);

        const Eigen::Quaterniond earthToGlobalRotation = bodies.at("Earth")->getCurrentRotationToGlobalFrame();

        const auto position = ephemeris->getCartesianPosition(t);
        const auto velocity = ephemeris->getCartesianVelocity(t);

        PaneledSourceRadiationPressureAcceleration accelerationModel(
            sourceModel,
            [] () { return Eigen::Vector3d::Zero(); },
            [=] () { return earthToGlobalRotation; },
            targetModel,
            [=] () { return position; },
            [] () { return Eigen::Quaterniond::Identity(); },
            [=] () { return bodyMass; },
                std::make_shared<NoOccultingBodyOccultationModel>()
        );

        sourceModel->updateMembers(t);
        accelerationModel.updateMembers(t);
        const auto acceleration = accelerationModel.getAcceleration();

        Eigen::Vector3d radialUnit = position.normalized();
        Eigen::Vector3d alongTrackUnit = (velocity - radialUnit * velocity.dot(radialUnit)).normalized();
        Eigen::Vector3d crossTrackUnit = radialUnit.cross(alongTrackUnit);

        auto radialAcceleration = acceleration.dot(radialUnit);
        auto alongTrackAcceleration = acceleration.dot(alongTrackUnit);
        auto crossTrackAcceleration = acceleration.dot(crossTrackUnit);

        // Compare with Knocke (1988), Fig. 4
        BOOST_CHECK_GT(radialAcceleration, 1.5e-10);
        BOOST_CHECK_LT(radialAcceleration, 3.9e-10);

        BOOST_CHECK_GT(alongTrackAcceleration, -4.0e-11);
        BOOST_CHECK_LT(alongTrackAcceleration, 4.0e-11);

        // Cross-track/normal acceleration from figure is less, but ours are a bit larger and more negative, likely
        // since we have a slightly different arc
        BOOST_CHECK_GT(crossTrackAcceleration, -3.0e-11);
        BOOST_CHECK_LT(crossTrackAcceleration, 1.0e-11);

        accelerationModel.resetCurrentTime( );
        accelerationModel.setPerpendicularSourceDirectionScaling( 1.5 );
        accelerationModel.setSourceDirectionScaling( 0.8 );
        accelerationModel.updateMembers(t);

        Eigen::Vector3d scaledAcceleration = accelerationModel.getAcceleration();
        Eigen::Vector3d directionToEarth = position.normalized( );

        Eigen::Vector3d originalAccelerationInEarthDirection = directionToEarth.dot( acceleration ) * directionToEarth;
        Eigen::Vector3d scaledAccelerationInEarthDirection = directionToEarth.dot( scaledAcceleration ) * directionToEarth;

        Eigen::Vector3d originalAccelerationPerpendicularToEarthDirection = acceleration - originalAccelerationInEarthDirection;
        Eigen::Vector3d scaledAccelerationPerpendicularToEarthDirection = scaledAcceleration - scaledAccelerationInEarthDirection;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( 1.5 * originalAccelerationInEarthDirection ), scaledAccelerationInEarthDirection, 1.0E-12 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( 0.8 * originalAccelerationPerpendicularToEarthDirection ), scaledAccelerationPerpendicularToEarthDirection, 1.0E-12 );

    }
}


BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
