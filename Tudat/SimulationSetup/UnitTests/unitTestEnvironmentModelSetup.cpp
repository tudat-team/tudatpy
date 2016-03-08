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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120202    K. Kumar          Separated from unitTestBasicMathematics.cpp into new
 *                                  Interpolators sub-directory.
 *      120529    E.A.G. Heeren     Boostified unit test.
 *      120615    T. Secretin       Minor layout changes.
 *      120716    D. Dirkx          Updated with interpolator architecture.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"
#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/SimulationSetup/createAtmosphereModel.h"
#include "Tudat/SimulationSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/createGravityField.h"
#include "Tudat/SimulationSetup/createRotationModel.h"
#include "Tudat/SimulationSetup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace spice_interface;
using namespace input_output;
using namespace reference_frames;

BOOST_AUTO_TEST_SUITE( test_environment_model_setup )

//! Test set up of atmosphere environment models.
BOOST_AUTO_TEST_CASE( test_atmosphereModelSetup )
{

    // Create settings for tabulated atmosphere.
    boost::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
            boost::make_shared< TabulatedAtmosphereSettings >(
                input_output::getTudatRootPath( ) + "/External/AtmosphereTables/" +
                "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );

    // Create settings for exponential atmosphere
    double densityScaleHeight = 8.0E3;
    double constantTemperature = 270.0;
    double densityAtZeroAltitude = 1.225;
    double specificGasConstant = 287.1;
    boost::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
            boost::make_shared< ExponentialAtmosphereSettings >(
                densityScaleHeight, constantTemperature,
                densityAtZeroAltitude, specificGasConstant );

    // Create atmpshere models using setup function
    boost::shared_ptr< aerodynamics::AtmosphereModel > exponentialAtmosphere =
            createAtmosphereModel( exponentialAtmosphereSettings, "Earth" );
    boost::shared_ptr< aerodynamics::AtmosphereModel > tabulatedAtmosphere =
            createAtmosphereModel( tabulatedAtmosphereSettings, "Earth" );

    // Create atmosphere models manually.
    aerodynamics::TabulatedAtmosphere manualTabulatedAtmosphere(
                input_output::getTudatRootPath( ) + "/External/AtmosphereTables/" +
                "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    aerodynamics::ExponentialAtmosphere manualExponentialAtmosphere(
                densityScaleHeight, constantTemperature, densityAtZeroAltitude,
                specificGasConstant );

    // Verify equivalence of automatically set up and manual models.
    BOOST_CHECK_EQUAL( manualTabulatedAtmosphere.getDensity( 32.0, 0.0, 0.0, 0.0 ),
                       tabulatedAtmosphere->getDensity( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualTabulatedAtmosphere.getPressure( 32.0, 0.0, 0.0, 0.0 ),
                       tabulatedAtmosphere->getPressure( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualTabulatedAtmosphere.getTemperature( 32.0, 0.0, 0.0, 0.0 ),
                       tabulatedAtmosphere->getTemperature( 32.0, 0.0, 0.0, 0.0 ) );

    BOOST_CHECK_EQUAL( manualExponentialAtmosphere.getDensity( 32.0, 0.0, 0.0, 0.0 ),
                       exponentialAtmosphere->getDensity( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualExponentialAtmosphere.getPressure( 32.0, 0.0, 0.0, 0.0 ),
                       exponentialAtmosphere->getPressure( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualExponentialAtmosphere.getTemperature( 32.0, 0.0, 0.0, 0.0 ),
                       exponentialAtmosphere->getTemperature( 32.0, 0.0, 0.0, 0.0 ) );

}

//! Test set up of ephemeris environment models.
BOOST_AUTO_TEST_CASE( test_ephemerisSetup )
{
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls" );

    {
        // Create settings for approximate planet positions.
        ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData bodyIdentifier =
                ephemerides::ApproximatePlanetPositionsBase::mars;
        bool useCircularCoplanarApproximation = 0;
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximateEphemerisSettings =
                boost::make_shared< ApproximatePlanetPositionSettings >(
                    bodyIdentifier, useCircularCoplanarApproximation );

        // Create ephemeris using setup function.
        boost::shared_ptr< ephemerides::Ephemeris > approximateEphemeris =
                createBodyEphemeris( approximateEphemerisSettings, "Earth" );

        // Create manual ephemeris.
        ephemerides::ApproximatePlanetPositions manualApproximateEphemeris(
                    bodyIdentifier );

        // Verify equivalence of automatically set up and manual models.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( manualApproximateEphemeris.getCartesianStateFromEphemeris( 1.0E7 ) ),
                    ( approximateEphemeris->getCartesianStateFromEphemeris( 1.0E7 ) ),
                    std::numeric_limits< double >::epsilon( ) );
    }

    {
        // Create spice ephemeris.
        boost::shared_ptr< EphemerisSettings > spiceEphemerisSettings =
                boost::make_shared< DirectSpiceEphemerisSettings >( "Earth", "J2000" );
        boost::shared_ptr< ephemerides::Ephemeris > spiceEphemeris =
                createBodyEphemeris( spiceEphemerisSettings, "Moon" );

        // Compare spice ephemeris against direct spice state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spice_interface::getBodyCartesianStateAtEpoch(
                          "Moon", "Earth", "J2000", "None", 1.0E7 ) ),
                    ( spiceEphemeris->getCartesianStateFromEphemeris( 1.0E7 ) ),
                    std::numeric_limits< double >::epsilon( ) );
    }

    {
        // Create tabulated spice ephemeris
        boost::shared_ptr< EphemerisSettings > spiceEphemerisSettings =
                boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                    1.0E7 - 50.0 * 600.0, 1.0E7 + 50.0 * 600.0, 600.0, "Earth", "J2000" );
        boost::shared_ptr< ephemerides::Ephemeris > spiceEphemeris =
                createBodyEphemeris( spiceEphemerisSettings, "Moon" );

        // Compare tabulated spice ephemeris against direct spice state on node point.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spice_interface::getBodyCartesianStateAtEpoch(
                          "Moon", "Earth", "J2000", "None", 1.0E7 ) ),
                    ( spiceEphemeris->getCartesianStateFromEphemeris( 1.0E7 ) ),
                    std::numeric_limits< double >::epsilon( ) );

        // Manually create table of states from spice
        std::map< double, basic_mathematics::Vector6d > tabulatedStates;
        double currentTime = 1.0E7 - 50.0 * 600.0;
        while( currentTime <= 1.0E7 + 50.0 * 600.0 )
        {
            tabulatedStates[ currentTime ] = spice_interface::getBodyCartesianStateAtEpoch(
                        "Moon", "Earth", "J2000", "None", currentTime );
            currentTime += 600.0;
        }

        // Create tabulated ephemeris.
        boost::shared_ptr< EphemerisSettings > tabulatedEphemerisSettings =
                boost::make_shared< TabulatedEphemerisSettings >(
                    tabulatedStates, "Earth", "J2000" );
        boost::shared_ptr< ephemerides::Ephemeris > tabulatedEphemeris =
                createBodyEphemeris( tabulatedEphemerisSettings, "Moon" );

        // Manually create tabulated ephemeris.
        boost::shared_ptr< ephemerides::Ephemeris > manualTabulatedEphemeris =
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::make_shared< interpolators::LagrangeInterpolator
                    < double, basic_mathematics::Vector6d > >( tabulatedStates, 6 ),
                    "Earth", "J2000" );

        // Compare ephemerides away from node point.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spiceEphemeris->getCartesianStateFromEphemeris( 1.0E7 + 110.0 ) ),
                    ( tabulatedEphemeris->getCartesianStateFromEphemeris( 1.0E7 + 110.0 ) ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spiceEphemeris->getCartesianStateFromEphemeris( 1.0E7 + 110.0 ) ),
                    ( manualTabulatedEphemeris->getCartesianStateFromEphemeris( 1.0E7 + 110.0 ) ),
                    std::numeric_limits< double >::epsilon( ) );
    }


}

//! Test set up of gravity field model environment models.
BOOST_AUTO_TEST_CASE( test_gravityFieldSetup )
{
    // Load Spice kernel with gravitational parameters.
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc" );

    // Create settings for spice central gravity field model.
    boost::shared_ptr< GravityFieldSettings > spiceCentralGravityFieldSettings =
            boost::make_shared< GravityFieldSettings >( central_spice );

    // Create spice central gravity field model from setup function.
    boost::shared_ptr< gravitation::GravityFieldModel > spiceCentralGravityField =
            createGravityFieldModel( spiceCentralGravityFieldSettings, "Venus" );

    // Check correct creation of gravity field.
    BOOST_CHECK_EQUAL(
                ( spice_interface::getBodyGravitationalParameter( "Venus" ) ),
                ( spiceCentralGravityField->getGravitationalParameter( ) ) );

    // Settings for spherical harmonic acceleration.
    double gravitationalParameter = 398600.4418E9;
    Eigen::Vector3d testPosition( 7.0e6, 8.0e6, 9.0e6 );
    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    // Create settings for central gravity field.
    boost::shared_ptr< CentralGravityFieldSettings > centralGravityFieldSettings =
            boost::make_shared< CentralGravityFieldSettings >( gravitationalParameter );

    // Create central gravity field with setup function.
    boost::shared_ptr< gravitation::GravityFieldModel > centralGravityField =
            createGravityFieldModel( centralGravityFieldSettings, "Earth" );

    // Create central gravity field manually.
    gravitation::GravityFieldModel manualCentralGravityField(
                gravitationalParameter );

    // Verify equivalence of automatically set up and manual models.
    BOOST_CHECK_EQUAL(
                ( manualCentralGravityField.getGravitationalParameter( ) ),
                ( centralGravityField->getGravitationalParameter( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualCentralGravityField.getGradientOfPotential( testPosition ) ),
                ( centralGravityField->getGradientOfPotential( testPosition ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Create settings for sh gravity field.
    boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > shGravityFieldSettings =
            boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients,
                "Earth_fixed" );

    // Create sh gravity field with setup function.
    boost::shared_ptr< gravitation::GravityFieldModel > shGravityField =
            createGravityFieldModel( shGravityFieldSettings, "Earth" );

    // Create sh gravity field manually.
    gravitation::SphericalHarmonicsGravityField manualShGravityField(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients );

    // Verify equivalence of automatically set up and manual models.
    BOOST_CHECK_EQUAL(
                ( manualShGravityField.getGravitationalParameter( ) ),
                ( shGravityField->getGravitationalParameter( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualShGravityField.getGradientOfPotential( testPosition ) ),
                ( shGravityField->getGradientOfPotential( testPosition ) ),
                std::numeric_limits< double >::epsilon( ) );

}

//! Test set up of rotation model environment models.
BOOST_AUTO_TEST_CASE( test_rotationModelSetup )
{

    // Create settings for simple rotation model.
    Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
    spiceInitialRotationToTargetFrameMatrix
            << -0.9548214974296336, 0.2665104385944917, 0.1314841974018291,
            -0.296591573568662, -0.882413772579987, -0.3652114078848295,
            0.01869081416890206, -0.3877088083617987, 0.9215923900425707;
    double venusRotationRate = unit_conversions::convertDegreesToRadians( -1.4813688 ) /
            physical_constants::JULIAN_DAY;
    boost::shared_ptr< SimpleRotationModelSettings > simpleRotationSettings =
            boost::make_shared< SimpleRotationModelSettings >
            ( "IAU_VENUS", "J2000", Eigen::Quaterniond( spiceInitialRotationToTargetFrameMatrix ),
              1.0E7, venusRotationRate );

    // Create rotation model using setup function
    boost::shared_ptr< ephemerides::RotationalEphemeris > approximateEphemeris =
            createRotationModel( simpleRotationSettings, "Earth" );

    // Create rotation model manually.
    ephemerides::SimpleRotationalEphemeris manualApproximateEphemeris(
                Eigen::Quaterniond( spiceInitialRotationToTargetFrameMatrix ),
                venusRotationRate, 1.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000,
                "J2000", "IAU_VENUS" );

    // Verify equivalence of automatically set up and manual models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::Matrix3d( manualApproximateEphemeris.getRotationToBaseFrame(
                                       4.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) ),
                ( Eigen::Matrix3d( approximateEphemeris->getRotationToBaseFrame(
                                       4.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::Matrix3d( manualApproximateEphemeris.getRotationToTargetFrame(
                                       4.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) ),
                ( Eigen::Matrix3d( approximateEphemeris->getRotationToTargetFrame(
                                       4.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) ),
                std::numeric_limits< double >::epsilon( ) );

}

//! Test set up of rotation model environment models.
BOOST_AUTO_TEST_CASE( test_RadiationPressureInterfaceSetup )
{

}

//! Test set up of rotation model environment models.
BOOST_AUTO_TEST_CASE( test_ShapeModelSetup )
{

}


//! Test set up of rotation model environment models.
BOOST_AUTO_TEST_CASE( test_flightConditionsSetup )
{

    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "naif0009.tls" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings(
                "Earth", 0.0, 1.0E7 );

    bodySettings[ "Vehicle" ] = boost::make_shared< BodySettings >( );

    bodySettings[ "Vehicle" ] ->aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                1.0, 2.0, 3.0, Eigen::Vector3d::Zero( ),
                ( Eigen::Vector3d( )<<-1.1, 0.1, 2.3 ).finished( ),
                Eigen::Vector3d::Zero( ), 1, 1 );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    double testHeadingAngle = 1.229357188236127;
    double testFlightPathAngle = -0.024894033070522;
    double testLatitude = -0.385027359562548;
    double testLongitude = -1.849449608688977;

    double angleOfAttack = 1.232;
    double angleOfSideslip = -0.00322;
    double bankAngle = 2.323432;

    boost::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
            createFlightConditions( bodyMap.at( "Vehicle" ), bodyMap.at( "Earth" ),
                                    "Vehicle", "Earth",
                                    boost::lambda::constant( angleOfAttack ),
                                    boost::lambda::constant( angleOfSideslip ),
                                    boost::lambda::constant( bankAngle ) );

    basic_mathematics::Vector6d vehicleBodyFixedState =
            ( basic_mathematics::Vector6d( )<< -1656517.23153109, -5790058.28764025, -2440584.88186829,
              6526.30784888051, -2661.34558272018, 2377.09572383163 ).finished( );

    double testTime = 0.5E7;

    basic_mathematics::Vector6d vehicleInertialState =
            ephemerides::transformStateToFrame(
                vehicleBodyFixedState,
                bodyMap[ "Earth" ]->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ),
            bodyMap[ "Earth" ]->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ) );



    bodyMap[ "Earth" ]->setState( basic_mathematics::Vector6d::Zero( ) );
    bodyMap[ "Vehicle" ]->setState( vehicleInertialState );
    bodyMap[ "Earth" ]->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );


    vehicleFlightConditions->updateConditions( );

    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( latitude_angle ) -
                           testLatitude), 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( longitude_angle ) -
                           testLongitude), 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( heading_angle ) -
                           testHeadingAngle), 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( flight_path_angle ) -
                           testFlightPathAngle), 10.0 * std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( angle_of_attack ) -
                           angleOfAttack ), 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( angle_of_sideslip ) -
                           angleOfSideslip), 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( bank_angle ) -
                           bankAngle), 10.0 * std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( vehicleFlightConditions->getCurrentBodyCenteredBodyFixedState( ), vehicleBodyFixedState,
                                       ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

