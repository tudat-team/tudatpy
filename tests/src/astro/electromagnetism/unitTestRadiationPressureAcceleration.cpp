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

using namespace tudat::electromagnetism;
using TargetPanel = PaneledRadiationPressureTargetModel::Panel;
using SourcePanel = PaneledRadiationSourceModel::Panel;

BOOST_AUTO_TEST_SUITE(test_radiation_pressure_acceleration)

//! Test acceleration with all unity values
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_Unity )
{
    const auto expectedAcceleration = Eigen::Vector3d::UnitX();

    // Set distance to speed of light to cancel to unity radiation pressure
    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [](double) { return physical_constants::SPEED_OF_LIGHT; }, 1);
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1);
    IsotropicPointSourceRadiationPressureAcceleration accelerationModel(
            sourceModel,
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
            []() { return Eigen::Vector3d(1, 0, 0); },
            targetModel,
            []() { return Eigen::Vector3d(2, 0, 0); },
            []() { return Eigen::Quaterniond::Identity(); },
            []() { return 1; },
            std::make_shared<NoOccultingBodyOccultationModel>());

    sourceModel->updateMembers(TUDAT_NAN);
    targetModel->updateMembers(TUDAT_NAN);
    accelerationModel.updateMembers(TUDAT_NAN);

    const auto actualAcceleration = accelerationModel.getAcceleration();

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-15);
}

//! Test acceleration for idealized GOCE spacecraft (Noomen 2022)
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_GOCE )
{
    const double expectedReceivedIrradiance = 1371;
    const Eigen::Vector3d expectedAcceleration = Eigen::Vector3d(1, 1, 0).normalized() * 5.2e-9;
    const double expectedReceivedFraction = 1.0;

    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [](double) { return 1371; }, physical_constants::ASTRONOMICAL_UNIT);
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
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
        std::make_shared<ConstantLuminosityModel>(1));

    const auto targetPosition = Eigen::Vector3d(10, 0, 0);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1);

    // Source and target keep their position, occulting body moves between tests
    {
        // Occulting body not interfering
        const double expectedReceivedFraction = 1.0;

        const auto occultingBodyPosition = Eigen::Vector3d(5, 5, 0);

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
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualReceivedFraction = accelerationModel.getSourceToTargetReceivedFraction();

        BOOST_CHECK(actualAcceleration > 0);
        BOOST_CHECK_CLOSE_FRACTION(actualReceivedFraction, expectedReceivedFraction, 1e-15);
    }

    {
        // Occulting body interfering with source -> target
        const double expectedReceivedFraction = 0.0;

        const auto occultingBodyPosition = Eigen::Vector3d(5, 0, 0);

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
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualReceivedFraction = accelerationModel.getSourceToTargetReceivedFraction();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK_CLOSE_FRACTION(actualReceivedFraction, expectedReceivedFraction, 1e-15);
    }
}

//! Test basic cases for paneled target acceleration with isotropic point source
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_IsotropicPointSource_PaneledTarget_Basic )
{
    // Set distance to speed of light to cancel to unity radiation pressure
    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [](double) { return physical_constants::SPEED_OF_LIGHT; }, 1);
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
    std::vector<PaneledRadiationPressureTargetModel::Panel> panels{
            TargetPanel(1, -Eigen::Vector3d::UnitX(),
                        reflectionLawFromAbsorptivityAndDiffuseReflectivity(1, 0)),
            TargetPanel(2, -Eigen::Vector3d::UnitY(),
                        reflectionLawFromSpecularAndDiffuseReflectivity(0, 1)),
            TargetPanel(3, Eigen::Vector3d::UnitX(), // never pointing towards source in these tests
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

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-15);
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


        std::vector<TargetPanel> panels;
        if ( testCase == 0 || testCase == 1 || testCase == 2 )
        {
            // Define panels properties for test cases with constant rotational model.
            panels = {
                    TargetPanel(2.0, Eigen::Vector3d::UnitZ(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(4.0, Eigen::Vector3d::UnitZ(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    TargetPanel(6.0, -Eigen::Vector3d::UnitZ(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(9.9, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(2.3, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    TargetPanel(9.9, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(2.3, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    TargetPanel(4.6, Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(2.7, Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    TargetPanel(5.8, -Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(2.7, -Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
            };
        }
        else if ( testCase == 3 )
        {
            // Define panels properties for test cases with constant rotational model (simpler boxes-and-wings model).
            panels = {
                    TargetPanel(9.9, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(2.3, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
                    TargetPanel(9.9, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.0, 0.06)),
                    TargetPanel(2.3, -Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.1, 0.46)),
            };
        }
        bodies.at( "Vehicle" )->setRadiationPressureTargetModel(
                std::make_shared<PaneledRadiationPressureTargetModel>(panels));

        std::vector< double > areas;
        std::vector< Eigen::Vector3d > panelSurfaceNormals;
        std::vector< double > specularReflectivities;
        std::vector< double > diffuseReflectivities;
        for (auto& panel : panels)
        {
            areas.push_back(panel.getArea());
            panelSurfaceNormals.push_back(panel.getSurfaceNormal());
            auto reflectionLaw = std::dynamic_pointer_cast<SpecularDiffuseMixReflectionLaw>(panel.getReflectionLaw());
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

            auto sourceIrradiance = radiationSourceModel->evaluateIrradianceAtPosition(bodies.at( "Vehicle" )->getPosition());
            auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;

            // Calculated accelerations.

            // Calculated acceleration if Sun-Vehicle vector aligned with +X-axis (acceleration generated by panels whose normals are
            // along +X-axis only)
            double cosPanelInclinationPositiveXaxis = expectedVehicleToSunVectorNormalized.dot(Eigen::Vector3d::UnitX( ) );
            Eigen::Vector3d accelerationPositiveXaxisOrientedPanels = - cosPanelInclinationPositiveXaxis
                * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                * ( areas[5] * (( 1.0 - specularReflectivities[5] ) * expectedVehicleToSunVectorNormalized
                    + 2.0 / 3.0 * diffuseReflectivities[5] * Eigen::Vector3d::UnitX( ) )
                    + areas[6] * (( 1.0 - specularReflectivities[6] ) * expectedVehicleToSunVectorNormalized
                    + ( 2.0 / 3.0 * diffuseReflectivities[6] + 2.0 * cosPanelInclinationPositiveXaxis * specularReflectivities[6] )
                    * Eigen::Vector3d::UnitX( ) ) );


            // Calculated acceleration generated by the panels whose normals are along -X-axis.
            double cosPanelInclinationNegativeXaxis = expectedVehicleToSunVectorNormalized.dot(- Eigen::Vector3d::UnitX() );
            Eigen::Vector3d accelerationNegativeXaxisOrientedPanels = - cosPanelInclinationNegativeXaxis
                * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                * ( areas[3] * (( 1.0 - specularReflectivities[3] ) * expectedVehicleToSunVectorNormalized
                    + 2.0 / 3.0 * diffuseReflectivities[3] * - Eigen::Vector3d::UnitX( ) )
                    + areas[4] * (( 1.0 - specularReflectivities[4] ) * expectedVehicleToSunVectorNormalized
                    + ( 2.0 / 3.0 * diffuseReflectivities[4] + 2.0 * specularReflectivities[4] * cosPanelInclinationNegativeXaxis )
                    * - Eigen::Vector3d::UnitX( ) ) );

            // Calculated acceleration generated by the panels whose normals are along +Y-axis.
            double cosPanelInclinationPositiveYaxis = expectedVehicleToSunVectorNormalized.dot(Eigen::Vector3d::UnitY() );
            Eigen::Vector3d accelerationPositiveYaxisOrientedPanels = - cosPanelInclinationPositiveYaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[7] * (( 1.0 - specularReflectivities[7] ) * expectedVehicleToSunVectorNormalized
                        + 2.0 / 3.0 * diffuseReflectivities[7] * Eigen::Vector3d::UnitY( ) )
                        + areas[8] * (( 1.0 - specularReflectivities[8] ) * expectedVehicleToSunVectorNormalized
                        + ( 2.0 / 3.0 * diffuseReflectivities[8] + 2.0 * specularReflectivities[8] * cosPanelInclinationPositiveYaxis )
                        * Eigen::Vector3d::UnitY( ) ) );

            // Calculated acceleration generated by the panels whose normals are along -Y-axis.
            double cosPanelInclinationNegativeYaxis = expectedVehicleToSunVectorNormalized.dot(- Eigen::Vector3d::UnitY() );
            Eigen::Vector3d accelerationNegativeYaxisOrientedPanels = - cosPanelInclinationNegativeYaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[9] * (( 1.0 - specularReflectivities[9] ) * expectedVehicleToSunVectorNormalized
                        + 2.0 / 3.0 * diffuseReflectivities[9] * - Eigen::Vector3d::UnitY( ) )
                        + areas[10] * (( 1.0 - specularReflectivities[10] ) * expectedVehicleToSunVectorNormalized
                        + ( 2.0 / 3.0 * diffuseReflectivities[10] + 2.0 * specularReflectivities[10] * cosPanelInclinationNegativeYaxis )
                        * - Eigen::Vector3d::UnitY( ) ) );

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

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(calculatedAcceleration, expectedAcceleration, 1e-10)
        }
    }

}

//! Test basic case for paneled target acceleration with paneled source
BOOST_AUTO_TEST_CASE( testRadiationPressureAcceleration_StaticallyPaneledSource_PaneledTarget_Basic )
{
    // Set distance to speed of light to cancel to unity radiation pressure
    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [](double) { return physical_constants::SPEED_OF_LIGHT; }, 1);
    auto originalSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);

    std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
    radiosityModels.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
            std::make_shared<ConstantSurfacePropertyDistribution>(1)));

    // Source is a single panel pointing in +X with purely diffuse reflection
    std::vector<SourcePanel> panels;
    panels.emplace_back(
            2, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitX(), std::move(radiosityModels));
    auto sourceModel = std::make_shared<StaticallyPaneledRadiationSourceModel>("", std::move(panels));

    std::vector<TargetPanel> targetPanels{
            TargetPanel(1, -Eigen::Vector3d::UnitX(),
                        reflectionLawFromAbsorptivityAndDiffuseReflectivity(1, 0)),
            TargetPanel(3, Eigen::Vector3d::UnitX(), // never pointing towards source in these tests
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
            originalSourceModel,
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
            [] () { return Eigen::Vector3d::UnitX(); },
            std::make_shared<NoOccultingBodyOccultationModel>(),
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

    std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
    radiosityModels.push_back(std::make_unique<AngleBasedThermalSourcePanelRadiosityModel>(
            1000, 1000, std::make_shared<ConstantSurfacePropertyDistribution>(1)));

    // Source is a single panel pointing in +X
    std::vector<SourcePanel> panels;
    panels.emplace_back(
            1, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitX(), std::move(radiosityModels));
    auto sourceModel = std::make_shared<StaticallyPaneledRadiationSourceModel>("", std::move(panels));

    const auto targetPosition = Eigen::Vector3d(10, 10, 0);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(1, 1);

    // Original source, source and target keep their position, occulting bodies move between tests
    {
        // Two occulting bodies but not interfering with radiation
        const unsigned int expectedVisibleSourcePanelCount = 1;
        const unsigned int expectedIlluminatedSourcePanelCount = 1;
        const unsigned int expectedVisibleAndIlluminatedSourcePanelCount = 1;

        const auto originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(-5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, -5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));
        auto originalSourceToSourceOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return originalSourceToSourceOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                originalSourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [=] () { return originalSourcePosition; },
                sourceToTargetOccultationModel, originalSourceToSourceOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleSourcePanelCount = accelerationModel.getVisibleSourcePanelCount();
        const auto actualIlluminatedSourcePanelCount = accelerationModel.getIlluminatedSourcePanelCount();
        const auto actualVisibleAndIlluminatedSourcePanelCount = accelerationModel.getVisibleAndIlluminatedSourcePanelCount();

        BOOST_CHECK(actualAcceleration > 0);
        BOOST_CHECK(actualVisibleSourcePanelCount == expectedVisibleSourcePanelCount);
        BOOST_CHECK(actualIlluminatedSourcePanelCount == expectedIlluminatedSourcePanelCount);
        BOOST_CHECK(actualVisibleAndIlluminatedSourcePanelCount == expectedVisibleAndIlluminatedSourcePanelCount);
    }

    {
        // Occulting body only interfering with original source -> source
        const unsigned int expectedVisibleSourcePanelCount = 1;
        const unsigned int expectedIlluminatedSourcePanelCount = 0;
        const unsigned int expectedVisibleAndIlluminatedSourcePanelCount = 0;

        const auto originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, -5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));
        auto originalSourceToSourceOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return originalSourceToSourceOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                originalSourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [=] () { return originalSourcePosition; },
                sourceToTargetOccultationModel, originalSourceToSourceOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleSourcePanelCount = accelerationModel.getVisibleSourcePanelCount();
        const auto actualIlluminatedSourcePanelCount = accelerationModel.getIlluminatedSourcePanelCount();
        const auto actualVisibleAndIlluminatedSourcePanelCount = accelerationModel.getVisibleAndIlluminatedSourcePanelCount();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK(actualVisibleSourcePanelCount == expectedVisibleSourcePanelCount);
        BOOST_CHECK(actualIlluminatedSourcePanelCount == expectedIlluminatedSourcePanelCount);
        BOOST_CHECK(actualVisibleAndIlluminatedSourcePanelCount == expectedVisibleAndIlluminatedSourcePanelCount);
    }

    {
        // Occulting body only interfering with source -> target
        const unsigned int expectedVisibleSourcePanelCount = 0;
        const unsigned int expectedIlluminatedSourcePanelCount = 1;
        const unsigned int expectedVisibleAndIlluminatedSourcePanelCount = 0;

        const auto originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(-5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, 5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));
        auto originalSourceToSourceOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return originalSourceToSourceOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                originalSourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [=] () { return originalSourcePosition; },
                sourceToTargetOccultationModel, originalSourceToSourceOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleSourcePanelCount = accelerationModel.getVisibleSourcePanelCount();
        const auto actualIlluminatedSourcePanelCount = accelerationModel.getIlluminatedSourcePanelCount();
        const auto actualVisibleAndIlluminatedSourcePanelCount = accelerationModel.getVisibleAndIlluminatedSourcePanelCount();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK(actualVisibleSourcePanelCount == expectedVisibleSourcePanelCount);
        BOOST_CHECK(actualIlluminatedSourcePanelCount == expectedIlluminatedSourcePanelCount);
        BOOST_CHECK(actualVisibleAndIlluminatedSourcePanelCount == expectedVisibleAndIlluminatedSourcePanelCount);
    }

    {
        // Occulting bodies interfering with both original source -> source and source -> target
        const unsigned int expectedVisibleSourcePanelCount = 0;
        const unsigned int expectedIlluminatedSourcePanelCount = 0;
        const unsigned int expectedVisibleAndIlluminatedSourcePanelCount = 0;

        const auto originalSourceToSourceOccultingBodyPosition = Eigen::Vector3d(5, 0, 0);
        const auto sourceToTargetOccultingBodyPosition = Eigen::Vector3d(5, 5, 0);

        auto sourceToTargetOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return sourceToTargetOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));
        auto originalSourceToSourceOccultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
            "", [=] () { return originalSourceToSourceOccultingBodyPosition; },
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1));

        PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                [] () { return Eigen::Vector3d::Zero(); },
                [] () { return Eigen::Quaterniond::Identity(); },
                targetModel,
                [=] () { return targetPosition; },
                [] () { return Eigen::Quaterniond::Identity(); },
                [] () { return 1; },
                originalSourceModel,
                std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(1),
                [=] () { return originalSourcePosition; },
                sourceToTargetOccultationModel, originalSourceToSourceOccultationModel);

        originalSourceModel->updateMembers(TUDAT_NAN);
        sourceModel->updateMembers(TUDAT_NAN);
        targetModel->updateMembers(TUDAT_NAN);
        accelerationModel.updateMembers(TUDAT_NAN);

        const auto actualAcceleration = accelerationModel.getAcceleration().norm();
        const auto actualVisibleSourcePanelCount = accelerationModel.getVisibleSourcePanelCount();
        const auto actualIlluminatedSourcePanelCount = accelerationModel.getIlluminatedSourcePanelCount();
        const auto actualVisibleAndIlluminatedSourcePanelCount = accelerationModel.getVisibleAndIlluminatedSourcePanelCount();

        BOOST_CHECK_CLOSE_FRACTION(actualAcceleration, 0, 1e-15);
        BOOST_CHECK(actualVisibleSourcePanelCount == expectedVisibleSourcePanelCount);
        BOOST_CHECK(actualIlluminatedSourcePanelCount == expectedIlluminatedSourcePanelCount);
        BOOST_CHECK(actualVisibleAndIlluminatedSourcePanelCount == expectedVisibleAndIlluminatedSourcePanelCount);
    }
}


BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
