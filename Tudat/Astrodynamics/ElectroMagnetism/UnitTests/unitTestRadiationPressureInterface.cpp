/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN


#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_radiation_pressure_interface )

//! Test implementation of radiation pressure calculation
BOOST_AUTO_TEST_CASE( testRadiationPressureCalculation )
{
    // Calculate total solar power (Montenbruck & Gill, 2000, p.77)
    double totalSolarPower = 1367.0 * 4.0 * mathematical_constants::PI *
            physical_constants::ASTRONOMICAL_UNIT * physical_constants::ASTRONOMICAL_UNIT;

    // Calculate radiation pressure at 1 AU.
    double calculatedRadiationPressure = electro_magnetism::calculateRadiationPressure(
                totalSolarPower, physical_constants::ASTRONOMICAL_UNIT );

    // Set literature radiation pressure (Montenbruck & Gill, 2000, p.77)
    double expectedRadiationPressure = 4.56E-6;

    BOOST_CHECK_CLOSE_FRACTION( calculatedRadiationPressure, expectedRadiationPressure, 1.0E-4 );

    // Test calculation of radiation pressure from class interface.
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface =
            boost::make_shared< electro_magnetism::RadiationPressureInterface >(
                boost::lambda::constant( totalSolarPower ),
                boost::lambda::constant(
                    ( Eigen::Vector3d( ) <<
                      0.0, physical_constants::ASTRONOMICAL_UNIT / std::sqrt( 2.0 ), 0.0 ).finished( ) ),
                boost::lambda::constant(
                    ( Eigen::Vector3d( ) <<
                      physical_constants::ASTRONOMICAL_UNIT / std::sqrt( 2.0 ), 0.0, 0.0 ).finished( ) ),
                1.0, 1.0 );
    radiationPressureInterface->updateInterface( );
    double classCalculatedRadiationPressure =
            radiationPressureInterface->getCurrentRadiationPressure( );

    BOOST_CHECK_CLOSE_FRACTION( classCalculatedRadiationPressure, expectedRadiationPressure,
                                1.0E-4 );
    BOOST_CHECK_CLOSE_FRACTION( classCalculatedRadiationPressure, calculatedRadiationPressure,
                                2.0 * std::numeric_limits< double >::epsilon( ) );

}

//! Test application of shadow function.
BOOST_AUTO_TEST_CASE( testShadowFunctionLink )
{
    // Set test geometry (see testShadowFunctionForPartialShadow for details)
    double totalSolarPower = 1367.0 * 4.0 * mathematical_constants::PI *
            physical_constants::ASTRONOMICAL_UNIT * physical_constants::ASTRONOMICAL_UNIT;

    const Eigen::Vector3d occultingBodyPosition = Eigen::Vector3d::Zero( );
    const double occultedBodyRadius = 6.96e8; // Siedelmann 1992.
    const double occultingBodyRadius = 6378.137e3; // WGS-84.

    std::vector< boost::function< Eigen::Vector3d( ) > > occultingBodyPositionFunctions;
    occultingBodyPositionFunctions.push_back( boost::lambda::constant( occultingBodyPosition ) );
    std::vector< double > occultingBodyRadii;
    occultingBodyRadii.push_back( occultingBodyRadius );

    Eigen::Vector3d satelliteDirection( 0.018, 1.0, 0.0 );
    satelliteDirection.normalize( );
    const Eigen::Vector3d satellitePosition = ( occultingBodyRadius + 1.0e3 ) * satelliteDirection;
    const Eigen::Vector3d occultedBodyPosition = -149598000.0e3 * Eigen::Vector3d( 1.0, 0.0, 0.0 );

    // Create radiation pressure interface with occultation.
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface >
            occultedRadiationPressureInterface =
            boost::make_shared< electro_magnetism::RadiationPressureInterface >(
                boost::lambda::constant( totalSolarPower ),
                boost::lambda::constant( occultedBodyPosition ),
                boost::lambda::constant( satellitePosition ),
                1.0, 1.0, occultingBodyPositionFunctions, occultingBodyRadii,
                occultedBodyRadius );
    occultedRadiationPressureInterface->updateInterface( );

    // Create radiation pressure interface without occultation.
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface >
            unoccultedRadiationPressureInterface =
            boost::make_shared< electro_magnetism::RadiationPressureInterface >(
                boost::lambda::constant( totalSolarPower ),
                boost::lambda::constant( occultedBodyPosition ),
                boost::lambda::constant( satellitePosition ),
                1.0, 1.0 );
    unoccultedRadiationPressureInterface->updateInterface( );

    // Test application of shadow function (see testShadowFunctionForPartialShadow)
    // Satellite partially visible. Satellite is located between the locations of the previous
    // tests in penumbra. According to analytical derivations in Matlab the shadow function should
    // be around 0.4547. Altitude of satellite = 1000 km (from unitTestMissionGeometry.cpp).
    BOOST_CHECK_CLOSE_FRACTION(
                occultedRadiationPressureInterface->getCurrentRadiationPressure( ) /
                unoccultedRadiationPressureInterface->getCurrentRadiationPressure( ) , 0.4547,
                                1.0E-4 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
