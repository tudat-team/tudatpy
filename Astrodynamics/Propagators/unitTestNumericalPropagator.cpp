/*! \file unitTestNumericalPropagator.cpp
 *    Source file that defines a unit test that tests the numerical propagator
 *    included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 16 February, 2011
 *    Last modified     : 16 Feburary, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      The unit test at present only checks that the code is internally
 *      consistent and doesn't check the result against benchmark data.
 *      Validation of the code against benchmark data should be added to
 *      ensure that the output of both simulation cases tested is correct.
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110216    K. Kumar          First creation of code.
 *      110216    E. Iorfida        Moved "Full Earth gravity case" as first
 *                                  case to test. Set the value of the
 *                                  gravitational parameter of "half Earth
 *                                  gravity case" with (pointerToEarth->
 *                                  getGravitationalParameter( ) / 2.0).
 */

// Include statements.
#include "unitTestNumericalPropagator.h"

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of implementation of numerical propagator class.
bool testNumericalPropagator( )
{
    // Test to see if the orbit of a satellite around the Earth is correctly
    // reproduced with 1x the gravity of the Earth and 2x half the gravity of
    // the Earth.

    // Test result initialised to false.
    bool isNumericalPropagatorErroneous = false;

    // Run full Earth gravity case.

    // Create pointer to the state of Asterix given in Cartesian elements.
    CartesianElements* pointerToStateOfAsterixForFullEarthGravity
            = new CartesianElements;

    // Fill initial state vector with position and
    // velocity given for Asterix.
    // Position is given in kilometers and
    // velocity is given in kilometers per second.
    pointerToStateOfAsterixForFullEarthGravity->setCartesianElementX( 7000.0 );
    pointerToStateOfAsterixForFullEarthGravity->setCartesianElementY( 0.0 );
    pointerToStateOfAsterixForFullEarthGravity->setCartesianElementZ( 0.0 );
    pointerToStateOfAsterixForFullEarthGravity->setCartesianElementXDot( 2.0 );
    pointerToStateOfAsterixForFullEarthGravity->setCartesianElementYDot( 5.0 );
    pointerToStateOfAsterixForFullEarthGravity->setCartesianElementZDot( 7.0 );

    // Convert initial state vector to meters from
    // kilometers for both cases.
    pointerToStateOfAsterixForFullEarthGravity->state =
            unit_conversions::convertKilometersToMeters(
                    pointerToStateOfAsterixForFullEarthGravity->state );

    // Create map of propagation history of Asterix for full gravity and half
    // gravity, and an iterator that works for both.
    std::map < double, State* > asterixPropagationHistoryFullGravity;

    // Create a pointer to new vehicle for Asterix.
    Vehicle* pointerToAsterixForFullEarthGravity = new Vehicle;

    // Create pre-defined Earth object for full gravity.
    CelestialBody* pointerToEarth
            = predefined_celestial_bodies::createPredefinedCelestialBody(
                    predefined_celestial_bodies::earth );

    // Create a pointer to a new Gravity object for full Earth gravity.
    Gravity* pointerToFullEarthGravity = new Gravity;

    // Set Earth as central body for full Earth gravity.
    pointerToFullEarthGravity->setBody( pointerToEarth );

    // Create a pointer to a new RK4 integrator.
    RungeKutta4thOrderFixedStepsize* pointerToRK4ForFullEarthGravity
            = new RungeKutta4thOrderFixedStepsize;

    // Set an initial stepsize for the integrators.
    pointerToRK4ForFullEarthGravity->setInitialStepsize( 30.0 );

    // Create numerical propagator object for full Earth.
    NumericalPropagator numericalPropagatorForFullEarthGravity;

    // Set fixed output interval for output in numerical
    // propagator object.
    numericalPropagatorForFullEarthGravity.setFixedOutputInterval( 3600.0 );

    // Set the propagation start time.
    numericalPropagatorForFullEarthGravity.setPropagationIntervalStart( 0.0 );

    // Set the propagation end time.
    numericalPropagatorForFullEarthGravity
            .setPropagationIntervalEnd( 86400.0 );

    // Set the integrator to use RK4.
    numericalPropagatorForFullEarthGravity
            .setIntegrator( pointerToRK4ForFullEarthGravity );

    // Add Asterix as the body that has to be propagated.
    numericalPropagatorForFullEarthGravity
            .addBody( pointerToAsterixForFullEarthGravity );

    // Add full Earth gravity as force acting on Asterix.
    numericalPropagatorForFullEarthGravity
            .addForceModel( pointerToAsterixForFullEarthGravity,
                            pointerToFullEarthGravity );

    // Set initial state of Asterix.
    numericalPropagatorForFullEarthGravity.setInitialState(
            pointerToAsterixForFullEarthGravity,
            pointerToStateOfAsterixForFullEarthGravity );

    // Run simulation for full Earth.
    numericalPropagatorForFullEarthGravity.propagate( );

    // Get propagation history of Asterix.
    asterixPropagationHistoryFullGravity
            = numericalPropagatorForFullEarthGravity
              .getPropagationHistoryAtFixedOutputIntervals(
                      pointerToAsterixForFullEarthGravity );

    // Run half Earth gravity case.

    // Create pointer to the state of Asterix given in Cartesian elements.
    CartesianElements* pointerToStateOfAsterixForHalfEarthGravity
            = new CartesianElements;

    // Fill initial state vector with position and
    // velocity given for Asterix.
    // Position is given in kilometers and
    // velocity is given in kilometers per second.
    pointerToStateOfAsterixForHalfEarthGravity->setCartesianElementX( 7000.0 );
    pointerToStateOfAsterixForHalfEarthGravity->setCartesianElementY( 0.0 );
    pointerToStateOfAsterixForHalfEarthGravity->setCartesianElementZ( 0.0 );
    pointerToStateOfAsterixForHalfEarthGravity->setCartesianElementXDot( 2.0 );
    pointerToStateOfAsterixForHalfEarthGravity->setCartesianElementYDot( 5.0 );
    pointerToStateOfAsterixForHalfEarthGravity->setCartesianElementZDot( 7.0 );

    // Convert initial state vector to meters from
    // kilometers for both cases.
    pointerToStateOfAsterixForHalfEarthGravity->state =
            unit_conversions::convertKilometersToMeters(
                    pointerToStateOfAsterixForHalfEarthGravity->state );

    // Create map of propagation history of Asterix for full gravity and half
    // gravity, and an iterator that works for both.
    std::map < double, State* > asterixPropagationHistoryHalfGravity;

    // Create a pointer to new vehicle for Asterix.
    Vehicle* pointerToAsterixForHalfEarthGravity = new Vehicle;

    // Create half Earth object for half gravity.
    CelestialBody* pointerToHalfEarth = new CelestialBody;

    // Create gravity field model for half Earth gravity.
    SphericalHarmonicsGravityField* pointerToHalfEarthGravityField
            = new SphericalHarmonicsGravityField;

    // Set half Earth gravity field to central gravity.
    pointerToHalfEarthGravityField->setDegreeOfExpansion( 0 );
    pointerToHalfEarthGravityField->setOrderOfExpansion( 0 );

    // Set origin of half Earth gravity field.
    CartesianPositionElements* pointerToOriginForHalfEarthGravityField
            = new CartesianPositionElements;
    pointerToOriginForHalfEarthGravityField->setCartesianElementX( 0.0 );
    pointerToOriginForHalfEarthGravityField->setCartesianElementY( 0.0 );
    pointerToOriginForHalfEarthGravityField->setCartesianElementZ( 0.0 );
    pointerToHalfEarthGravityField
            ->setOrigin( pointerToOriginForHalfEarthGravityField );

    // Set gravity parameter for half Earth gravity field.
    pointerToHalfEarthGravityField->setGravitationalParameter(
            pointerToEarth->getGravitationalParameter( ) / 2.0 );

    // Set gravity field model for half Earth.
    pointerToHalfEarth->setGravityFieldModel( pointerToHalfEarthGravityField );

    // Create a pointer to a new Gravity object for half Earth gravity.
    Gravity* pointerToHalfEarthGravity = new Gravity;

    // Set half Earth as central body for half Earth gravity.
    pointerToHalfEarthGravity->setBody( pointerToHalfEarth );

    // Create a pointer to a new RK4 integrator.
    RungeKutta4thOrderFixedStepsize* pointerToRK4ForHalfEarthGravity
            = new RungeKutta4thOrderFixedStepsize;

    // Set an initial stepsize for the integrators.
    pointerToRK4ForHalfEarthGravity->setInitialStepsize( 30.0 );

    // Create numerical propagator object for full Earth.
    NumericalPropagator numericalPropagatorForHalfEarthGravity;

    // Set fixed output interval for output in numerical
    // propagator object.
    numericalPropagatorForHalfEarthGravity.setFixedOutputInterval( 3600.0 );

    // Set the propagation start time.
    numericalPropagatorForHalfEarthGravity.setPropagationIntervalStart( 0.0 );

    // Set the propagation end time.
    numericalPropagatorForHalfEarthGravity
            .setPropagationIntervalEnd( 86400.0 );

    // Set the integrator to use RK4.
    numericalPropagatorForHalfEarthGravity
            .setIntegrator( pointerToRK4ForHalfEarthGravity );

    // Add Asterix as the body that has to be propagated.
    numericalPropagatorForHalfEarthGravity
            .addBody( pointerToAsterixForHalfEarthGravity );

    // Add half Earth gravity as force acting on Asterix twice.
    numericalPropagatorForHalfEarthGravity
            .addForceModel( pointerToAsterixForHalfEarthGravity,
                            pointerToHalfEarthGravity );
    numericalPropagatorForHalfEarthGravity
            .addForceModel( pointerToAsterixForHalfEarthGravity,
                            pointerToHalfEarthGravity );

    // Set initial state of Asterix.
    numericalPropagatorForHalfEarthGravity.setInitialState(
            pointerToAsterixForHalfEarthGravity,
            pointerToStateOfAsterixForHalfEarthGravity );

    // Run simulation for half Earth.
    numericalPropagatorForHalfEarthGravity.propagate( );

    // Get propagation history of Asterix for half Earth gravity case.
    asterixPropagationHistoryHalfGravity
            = numericalPropagatorForHalfEarthGravity
              .getPropagationHistoryAtFixedOutputIntervals(
                      pointerToAsterixForHalfEarthGravity );

    // Check if results of full and half Earth gravity cases match.
    for ( unsigned int i = 0;
          i < numericalPropagatorForHalfEarthGravity
              .getPropagationIntervalEnd( )
              / numericalPropagatorForHalfEarthGravity
              .getFixedOutputInterval( ); i++ )
    {
        if ( asterixPropagationHistoryFullGravity[
                i * numericalPropagatorForFullEarthGravity
                .getFixedOutputInterval( ) ]->state
             != asterixPropagationHistoryHalfGravity[
                     i * numericalPropagatorForHalfEarthGravity
                     .getFixedOutputInterval( ) ]->state )
        {
            isNumericalPropagatorErroneous = true;

            std::cerr << "The numerical propagator does not produce "
                      << "consistent results, as running a simulation with "
                      << "full Earth gravity does not yield the same results "
                      << "as running a simulation with twice half Earth "
                      << "gravity." << std::endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isNumericalPropagatorErroneous;
}

}

// End of file.
