/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Basics/testMacros.h"

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#endif

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
using namespace input_output;
using namespace reference_frames;
using namespace basic_astrodynamics;

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

#if USE_CSPICE
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
#endif

#if USE_CSPICE
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

    boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > defaultEarthField =
            boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                createGravityFieldModel( getDefaultGravityFieldSettings(
                                         "Earth", TUDAT_NAN, TUDAT_NAN ), "Earth" ) );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getGravitationalParameter( ) ), ( 0.3986004418E15 ) );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getReferenceRadius( ) ), ( 6378137.0 ) );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( ).rows( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( ).cols( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getSineCoefficients( ).rows( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getSineCoefficients( ).cols( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( )( 2, 0 ) ), -0.484165371736E-03 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( )( 5, 3 ) ), -0.451955406071E-06 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getSineCoefficients( )( 7, 1 ) ), 0.954336911867E-07 );

    boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > defaultMoonField =
            boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                createGravityFieldModel( getDefaultGravityFieldSettings(
                                         "Moon", TUDAT_NAN, TUDAT_NAN ), "Moon" ) );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getGravitationalParameter( ) ), ( 0.4902800238000000E+13 ) );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getReferenceRadius( ) ), ( 0.17380E+07 ) );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getCosineCoefficients( ).rows( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getCosineCoefficients( ).cols( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getSineCoefficients( ).rows( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getSineCoefficients( ).cols( ) ), 51 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getCosineCoefficients( )( 5, 3 ) ), 0.5493176535439800E-06 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getSineCoefficients( )( 7, 1 ) ), -0.1744763377093700E-06 );

}
#endif

#if USE_CSPICE
//! Test set up of gravity field model variations environment models.
BOOST_AUTO_TEST_CASE( test_gravityFieldVariationSetup )
{
    // Load Spice kernel with gravitational parameters.
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "naif0009.tls" );


    // Settings for spherical harmonic acceleration.
    double gravitationalParameter = 398600.4418E9;
    double referenceRadius = 6378.0E3;
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


    // Define body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Moon" ] = getDefaultSingleBodySettings( "Moon", 0.0, 1.0E7 );
    bodySettings[ "Sun" ] = getDefaultSingleBodySettings( "Sun", 0.0, 1.0E7 );
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", 0.0, 1.0E7 );

    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, "IAU_Earth" );


    // Define settings for tidal variations.
    double loveNumber = 0.25;
    std::vector< std::vector< std::complex< double > > > fullLoveNumberVector =
            gravitation::getFullLoveNumbersVector( loveNumber, 3, 2 );
    double testTime = 0.5E7;


    Eigen::MatrixXd cosineCorrections1, cosineCorrections2, cosineCorrections3;
    Eigen::MatrixXd sineCorrections1, sineCorrections2, sineCorrections3;

    // Calculate non-interpolated corrections from two bodies in single object.
    {
        // Define correction settings
        std::vector< std::string > deformingBodies;
        deformingBodies.push_back( "Moon" );
        deformingBodies.push_back( "Sun" );

        bodySettings[ "Earth" ]->gravityFieldVariationSettings.push_back(
                    boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector, referenceRadius ) );

        // Create bodies
        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Update states.
        bodyMap[ "Earth" ]->setStateFromEphemeris( testTime );
        bodyMap[ "Earth" ]->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
        bodyMap[ "Sun" ]->setStateFromEphemeris( testTime );
        bodyMap[ "Moon" ]->setStateFromEphemeris( testTime );

        // Update gravity field
        boost::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                boost::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodyMap[ "Earth" ]->getGravityFieldModel( ) );
        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCorrections1 = earthGravityField->getCosineCoefficients( ) - cosineCoefficients;
        sineCorrections1 = earthGravityField->getSineCoefficients( ) - sineCoefficients;


    }

    // Calculate non-interpolated corrections from two bodies in separate objects.
    {
        bodySettings[ "Earth" ]->gravityFieldVariationSettings.clear( );

        // Define correction settings
        std::vector< std::string > deformingBodies;
        deformingBodies.push_back( "Moon" );
        bodySettings[ "Earth" ]->gravityFieldVariationSettings.push_back(
                    boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector, referenceRadius ) );

        deformingBodies.clear( );
        deformingBodies.push_back( "Sun" );
        bodySettings[ "Earth" ]->gravityFieldVariationSettings.push_back(
                    boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector, referenceRadius ) );

        // Create bodies
        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Update states.
        bodyMap[ "Earth" ]->setStateFromEphemeris( testTime );
        bodyMap[ "Earth" ]->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
        bodyMap[ "Sun" ]->setStateFromEphemeris( testTime );
        bodyMap[ "Moon" ]->setStateFromEphemeris( testTime );

        // Update gravity field
        boost::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                boost::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodyMap[ "Earth" ]->getGravityFieldModel( ) );
        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCorrections2 = earthGravityField->getCosineCoefficients( ) - cosineCoefficients;
        sineCorrections2 = earthGravityField->getSineCoefficients( ) - sineCoefficients;
    }

    {
        bodySettings[ "Earth" ]->gravityFieldVariationSettings.clear( );

        // Define correction settings
        std::vector< std::string > deformingBodies;
        deformingBodies.push_back( "Moon" );
        deformingBodies.push_back( "Sun" );
        bodySettings[ "Earth" ]->gravityFieldVariationSettings.push_back(
                    boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector, referenceRadius,
                        boost::make_shared< ModelInterpolationSettings >(
                            0.25E7, 0.75E7, 600.0,
                            boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ) ) );

        // Create bodies
        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Update gravity field
        boost::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                boost::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodyMap[ "Earth" ]->getGravityFieldModel( ) );
        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCorrections3 = earthGravityField->getCosineCoefficients( ) - cosineCoefficients;
        sineCorrections3= earthGravityField->getSineCoefficients( ) - sineCoefficients;

    }

    // Mutually compare results of three methods.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cosineCorrections1.block( 2, 0, 2, 3 ),  cosineCorrections2.block( 2, 0, 2, 3 ), 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cosineCorrections1.block( 2, 0, 2, 3 ),  cosineCorrections3.block( 2, 0, 2, 3 ), 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sineCorrections1.block( 2, 1, 2, 3 ),  sineCorrections2.block( 2, 1, 2, 3 ), 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sineCorrections1.block( 2, 1, 2, 3 ),  sineCorrections2.block( 2, 1, 2, 3 ), 1.0E-10 );

    // Check whether 0 coefficients are actually 0.
    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; j < 6; j++ )
        {
            if( i < 2 || i > 3 || j > 2 )
            {
                BOOST_CHECK_EQUAL( cosineCorrections1( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( cosineCorrections2( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( cosineCorrections3( i, j ), 0.0 );

                BOOST_CHECK_EQUAL( sineCorrections1( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections2( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections3( i, j ), 0.0 );
            }
            else if( j == 0 )
            {
                BOOST_CHECK_EQUAL( sineCorrections1( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections2( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections3( i, j ), 0.0 );
            }

        }
    }

    // Calculate corrections manually and compare against created results.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directMoonTide =
            gravitation::calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                fullLoveNumberVector, 0.4902800238000000E+13 / gravitationalParameter,
                referenceRadius,  spice_interface::getBodyCartesianPositionAtEpoch(
                    "Moon", "Earth", "IAU_Earth", "None", testTime ), 3, 2 );
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directSunTide =
            gravitation::calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                fullLoveNumberVector, spice_interface::getBodyGravitationalParameter( "Sun" ) / gravitationalParameter,
                referenceRadius,  spice_interface::getBodyCartesianPositionAtEpoch(
                    "Sun", "Earth", "IAU_Earth", "None", testTime ), 3, 2 );

    for( unsigned int n = 2; n <= 3; n++ )
    {
        for( unsigned m = 0; m <=2; m++ )
        {
            BOOST_CHECK_SMALL( directMoonTide.first( n, m ) + directSunTide.first( n, m ) - cosineCorrections1( n, m ), 1.0E-20 );
            BOOST_CHECK_SMALL( directMoonTide.second( n, m ) + directSunTide.second( n, m ) - sineCorrections1( n, m ), 1.0E-20 );
        }
    }
}
#endif

#if USE_CSPICE
//! Test set up of rotation models.
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
#endif

#if USE_CSPICE
//! Test set up of radiation pressure interfacel environment models.
BOOST_AUTO_TEST_CASE( test_radiationPressureInterfaceSetup )
{
    // Load Spice kernels
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "naif0009.tls" );

    // Define body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", 0.0, 1.0E7 );
    bodySettings[ "Sun" ] = getDefaultSingleBodySettings( "Sun", 0.0, 1.0E7 );

    // Get settings for vehicle
    double area = 2.34;
    double coefficient = 1.2;
    basic_mathematics::Vector6d initialKeplerElements =
            ( basic_mathematics::Vector6d( ) << 12000.0E3, 0.13, 0.3, 0.0, 0.0, 0.0 ).finished( );
    bodySettings[ "Vehicle" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Vehicle" ]->radiationPressureSettings[ "Sun" ] =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", area, coefficient );
    bodySettings[ "Vehicle" ]->ephemerisSettings =
            boost::make_shared< KeplerEphemerisSettings >(
                initialKeplerElements,
                0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    BOOST_CHECK_EQUAL( bodyMap[ "Vehicle" ]->getRadiationPressureInterfaces( ).size( ), 1 );
    BOOST_CHECK_EQUAL( bodyMap[ "Vehicle" ]->getRadiationPressureInterfaces( ).count( "Sun" ), 1 );

    double testTime = 0.5E7;

    // Update environment to current time.
    bodyMap[ "Sun" ]->setTemplatedStateFromEphemeris< double, double >( testTime );
    bodyMap[ "Earth" ]->setTemplatedStateFromEphemeris< double, double >( testTime );
    bodyMap[ "Vehicle" ]->setTemplatedStateFromEphemeris< double, double >( testTime );

    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > vehicleRadiationPressureInterface =
            bodyMap[ "Vehicle" ]->getRadiationPressureInterfaces( ).at( "Sun" );

    vehicleRadiationPressureInterface->updateInterface( testTime );
    double sourceDistance = ( ( bodyMap[ "Vehicle" ]->getState( ) -  bodyMap[ "Sun" ]->getState( ) ).
            segment( 0, 3 ) ).norm( );
    double expectedRadiationPressure = electro_magnetism::calculateRadiationPressure(
                defaultRadiatedPowerValues.at( "Sun" ), sourceDistance );

    BOOST_CHECK_CLOSE_FRACTION( expectedRadiationPressure,
                                vehicleRadiationPressureInterface->getCurrentRadiationPressure( ),
                                std::numeric_limits< double >::epsilon( ) );

}
#endif

#if USE_CSPICE
//! Test set up of body shape environment models (see testShapeModels).
BOOST_AUTO_TEST_CASE( test_shapeModelSetup )
{
    // Define test values.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    // Test spherical body setup
    {
        boost::shared_ptr< BodyShapeSettings > shapeSettings = boost::make_shared< SphericalBodyShapeSettings >(
                    equatorialRadius );
        boost::shared_ptr< BodyShapeModel > shapeModel = createBodyShapeModel(
                    shapeSettings, "Earth" );


        double calculatedAltitude = shapeModel->getAltitude(
                    testCartesianPosition );
        BOOST_CHECK_CLOSE_FRACTION( calculatedAltitude, testCartesianPosition.norm( ) - equatorialRadius,
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( shapeModel->getAverageRadius( ), equatorialRadius,
                                    std::numeric_limits< double >::epsilon( ) );
    }
    // Test oblate spheroid setup
    {
        boost::shared_ptr< BodyShapeSettings > shapeSettings = boost::make_shared< OblateSphericalBodyShapeSettings >(
                    equatorialRadius, flattening );
        boost::shared_ptr< BodyShapeModel > shapeModel = createBodyShapeModel(
                    shapeSettings, "Earth" );


        double calculatedAltitude = shapeModel->getAltitude(
                    testCartesianPosition );
        double manualAltitude = coordinate_conversions::calculateAltitudeOverOblateSpheroid(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );
        BOOST_CHECK_CLOSE_FRACTION( calculatedAltitude, manualAltitude,
                                    std::numeric_limits< double >::epsilon( ) );
    }

}
#endif

#if USE_CSPICE
//! Test set up of flight conditions object.
BOOST_AUTO_TEST_CASE( test_flightConditionsSetup )
{
    // Load Spice kernels
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "naif0009.tls" );

    // Define body settings/
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings(
                "Earth", 0.0, 1.0E7 );
    bodySettings[ "Vehicle" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Vehicle" ] ->aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                1.0, 2.0, 3.0, Eigen::Vector3d::Zero( ),
                ( Eigen::Vector3d( )<<-1.1, 0.1, 2.3 ).finished( ),
                Eigen::Vector3d::Zero( ), 1, 1 );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define expected aerodynamic angles (see testAerodynamicAngleCalculator)
    double testHeadingAngle = 1.229357188236127;
    double testFlightPathAngle = -0.024894033070522;
    double testLatitude = -0.385027359562548;
    double testLongitude = -1.849449608688977;

    double angleOfAttack = 1.232;
    double angleOfSideslip = -0.00322;
    double bankAngle = 2.323432;

    // Create flight conditions object.
    boost::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
            createFlightConditions( bodyMap.at( "Vehicle" ), bodyMap.at( "Earth" ),
                                    "Vehicle", "Earth",
                                    boost::lambda::constant( angleOfAttack ),
                                    boost::lambda::constant( angleOfSideslip ),
                                    boost::lambda::constant( bankAngle ) );

    // Set vehicle body-fixed state (see testAerodynamicAngleCalculator)
    basic_mathematics::Vector6d vehicleBodyFixedState =
            ( basic_mathematics::Vector6d( )<< -1656517.23153109, -5790058.28764025, -2440584.88186829,
              6526.30784888051, -2661.34558272018, 2377.09572383163 ).finished( );

    // Set states in environment.
    double testTime = 0.5E7;
    basic_mathematics::Vector6d vehicleInertialState =
            ephemerides::transformStateToFrame(
                vehicleBodyFixedState,
                bodyMap[ "Earth" ]->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ),
            bodyMap[ "Earth" ]->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ) );
    bodyMap[ "Earth" ]->setState( basic_mathematics::Vector6d::Zero( ) );
    bodyMap[ "Vehicle" ]->setState( vehicleInertialState );
    bodyMap[ "Earth" ]->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );

    // Update flight conditions
    vehicleFlightConditions->updateConditions( testTime );

    // Check whether calulcate angles and body-fixed state correspond to expected results
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( latitude_angle ) - testLatitude),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( longitude_angle ) - testLongitude),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( heading_angle ) - testHeadingAngle),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( flight_path_angle) - testFlightPathAngle),
                10.0 * std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( angle_of_attack ) - angleOfAttack ),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( angle_of_sideslip) - angleOfSideslip),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( bank_angle ) - bankAngle),
                10.0 * std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                vehicleFlightConditions->getCurrentBodyCenteredBodyFixedState( ),vehicleBodyFixedState,
                ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );


}
#endif

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

