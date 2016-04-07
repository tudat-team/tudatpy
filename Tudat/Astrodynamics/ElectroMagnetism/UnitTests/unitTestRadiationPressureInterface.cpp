/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      150408    D. Dirkx          File created.
 *
 *    References
 *
 *    Notes
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
                    ( Eigen::Vector3d( )<<
                      0.0, physical_constants::ASTRONOMICAL_UNIT / std::sqrt( 2.0 ), 0.0 ).finished( ) ),
                boost::lambda::constant(
                    ( Eigen::Vector3d( )<<
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
