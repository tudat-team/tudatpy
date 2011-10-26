/*! \file unitTestNumericalPropagator.cpp
 *    Source file that defines a unit test that tests the numerical propagator included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 5
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
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 16 February, 2011
 *    Last modified     : 20 September, 2011
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
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110216    E. Iorfida        Moved "Full Earth gravity case" as first case to test. Set the
 *                                  value of the gravitational parameter of "half-Earth gravity
 *                                  case" with pointerToEarth->getGravitationalParameter( )/2.0.
 *      110602    K. Kumar          Updated code to not use dynamic memory allocation and to work
 *                                  with other code updates.
 *      110815    K. Kumar          Included the mass of Asterix.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

// Include statements.
#include <iostream>
#include <map>
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/Bodies/planet.h"
#include "Astrodynamics/Bodies/vehicle.h"
#include "Astrodynamics/EnvironmentModels/sphericalHarmonicsGravityField.h"
#include "Astrodynamics/ForceModels/gravitationalForceModel.h"
#include "Astrodynamics/Propagators/cartesianStateNumericalPropagator.h"
#include "Astrodynamics/Propagators/numericalPropagator.h"
#include "Astrodynamics/Propagators/seriesPropagator.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Astrodynamics/States/cartesianPositionElements.h"
#include "Astrodynamics/States/state.h"
#include "Mathematics/unitConversions.h"
#include "Mathematics/NumericalIntegrators/rungeKutta4thOrderFixedStepsize.h"

//! Test implementation of numerical propagator class.
int main( )
{
    // Using declarations.
    using namespace tudat;

    // Test to see if the orbit of a satellite around the Earth is reproduced
    // the same with 1x the gravity of the Earth and 2x half the gravity of
    // the Earth. In other words, the results from the two cases are compared
    // directly.

    // Using declarations.
    using tudat::unit_conversions::convertKilometersToMeters;

    // Test result initialised to false.
    bool isNumericalPropagatorErroneous = false;

    // Run full Earth gravity case.

    // Create pointer to the state of Asterix given in Cartesian elements.
    CartesianElements stateOfAsterixForFullEarthGravity;

    // Fill initial state vector with position and
    // velocity given for Asterix.
    // Position is given in kilometers and
    // velocity is given in kilometers per second.
    stateOfAsterixForFullEarthGravity.setCartesianElementX( 7000.0 );
    stateOfAsterixForFullEarthGravity.setCartesianElementY( 0.0 );
    stateOfAsterixForFullEarthGravity.setCartesianElementZ( 0.0 );
    stateOfAsterixForFullEarthGravity.setCartesianElementXDot( 2.0 );
    stateOfAsterixForFullEarthGravity.setCartesianElementYDot( 5.0 );
    stateOfAsterixForFullEarthGravity.setCartesianElementZDot( 7.0 );

    // Convert initial state vector to meters from
    // kilometers for both cases.
    stateOfAsterixForFullEarthGravity.state = convertKilometersToMeters(
                stateOfAsterixForFullEarthGravity.state );

    // Create map of propagation history of Asterix for full gravity and half
    // gravity.
    std::map< double, State > asterixPropagationHistoryFullGravity;

    // Create a vehicle object for Asterix.
    Vehicle asterixForFullEarthGravity;

    // Set mass of Asterix.
    asterixForFullEarthGravity.setMass( 1000.0 );

    // Create predefined Earth object for full gravity.
    Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( Planet::earth );

    // Create a pointer to a gravitational force model for full Earth gravity.
    GravitationalForceModel fullEarthGravitationalForceModel;

    // Set Asterix as body subject to the full Earth gravitational force.
    fullEarthGravitationalForceModel.setBodySubjectToForce( &asterixForFullEarthGravity );

    // Set Earth as central body for full Earth gravity.
    fullEarthGravitationalForceModel.setGravitationalBody( &predefinedEarth );

    // Create a pointer to a new RK4 integrator.
    RungeKutta4thOrderFixedStepsize rk4ForFullEarthGravity;

    // Set an initial stepsize for the integrators.
    rk4ForFullEarthGravity.setInitialStepsize( 30.0 );

    // Create numerical propagator object for full Earth gravity.
    CartesianStateNumericalPropagator cartesianStateNumericalPropagatorForFullEarthGravity;

    // Set Cartesian state numerical propagator as state derivative class for
    // RK4 integrator.
    rk4ForFullEarthGravity.setObjectContainingStateDerivative(
                &cartesianStateNumericalPropagatorForFullEarthGravity );

    // Set the integrator to use RK4.
    cartesianStateNumericalPropagatorForFullEarthGravity.setIntegrator( &rk4ForFullEarthGravity );

    // Add Asterix as the body that has to be propagated.
    cartesianStateNumericalPropagatorForFullEarthGravity.addBody( &asterixForFullEarthGravity );

    // Add full Earth gravity as force acting on Asterix.
    cartesianStateNumericalPropagatorForFullEarthGravity.addForceModel(
                &asterixForFullEarthGravity, &fullEarthGravitationalForceModel );

    // Create series propagator for full Earth gravity.
    SeriesPropagator seriesPropagatorForFullEarthGravity;

    // Set series propagation start time.
    seriesPropagatorForFullEarthGravity.setSeriesPropagationStart( 0.0 );

    // Set the propagation end time.
    seriesPropagatorForFullEarthGravity.setSeriesPropagationEnd( 86400.0 );

    // Set fixed output interval for series propagation
    seriesPropagatorForFullEarthGravity.setFixedOutputInterval( 3600.0 );

    // Set propagator for series propagation.
    seriesPropagatorForFullEarthGravity.setPropagator(
            &cartesianStateNumericalPropagatorForFullEarthGravity );

    // Set initial state of Asterix for series propagation.
    seriesPropagatorForFullEarthGravity.setInitialState( &asterixForFullEarthGravity,
                                                         &stateOfAsterixForFullEarthGravity );

    // Run simulation for full Earth.
    seriesPropagatorForFullEarthGravity.execute( );

    // Get propagation history of Asterix.
    asterixPropagationHistoryFullGravity = seriesPropagatorForFullEarthGravity
              .getPropagationHistoryAtFixedOutputIntervals( &asterixForFullEarthGravity );

    // Run half Earth gravity case.

    // Create pointer to the state of Asterix given in Cartesian elements.
    CartesianElements stateOfAsterixForHalfEarthGravity;

    // Fill initial state vector with position and
    // velocity given for Asterix.
    // Position is given in kilometers and
    // velocity is given in kilometers per second.
    stateOfAsterixForHalfEarthGravity.setCartesianElementX( 7000.0 );
    stateOfAsterixForHalfEarthGravity.setCartesianElementY( 0.0 );
    stateOfAsterixForHalfEarthGravity.setCartesianElementZ( 0.0 );
    stateOfAsterixForHalfEarthGravity.setCartesianElementXDot( 2.0 );
    stateOfAsterixForHalfEarthGravity.setCartesianElementYDot( 5.0 );
    stateOfAsterixForHalfEarthGravity.setCartesianElementZDot( 7.0 );

    // Convert initial state vector to meters from
    // kilometers for both cases.
    stateOfAsterixForHalfEarthGravity.state = convertKilometersToMeters(
                stateOfAsterixForHalfEarthGravity.state );

    // Create map of propagation history of Asterix for full gravity and half
    // gravity.
    std::map < double, State > asterixPropagationHistoryHalfGravity;

    // Create a pointer to new vehicle for Asterix.
    Vehicle asterixForHalfEarthGravity;

    // Set mass of Asterix.
    asterixForHalfEarthGravity.setMass( 1000.0 );

    // Create half Earth object for half gravity.
    CelestialBody halfEarth;

    // Create gravity field model for half Earth gravity.
    SphericalHarmonicsGravityField halfEarthGravityField;

    // Set half Earth gravity field to central gravity.
    halfEarthGravityField.setDegreeOfExpansion( 0 );
    halfEarthGravityField.setOrderOfExpansion( 0 );

    // Set origin of half Earth gravity field.
    CartesianPositionElements originForHalfEarthGravityField;
    originForHalfEarthGravityField.setCartesianElementX( 0.0 );
    originForHalfEarthGravityField.setCartesianElementY( 0.0 );
    originForHalfEarthGravityField.setCartesianElementZ( 0.0 );
    halfEarthGravityField.setOrigin( &originForHalfEarthGravityField );

    // Set gravity parameter for half Earth gravity field.
    halfEarthGravityField.setGravitationalParameter( predefinedEarth
                                                     .getGravitationalParameter( ) / 2.0 );

    // Set gravity field model for half Earth.
    halfEarth.setGravityFieldModel( &halfEarthGravityField );

    // Create a pointer to a gravitational force model for half Earth gravity.
    GravitationalForceModel halfEarthGravitationalForceModel;

    // Set Asterix as body subject to half of the Earth gravitational force.
    halfEarthGravitationalForceModel.setBodySubjectToForce( &asterixForHalfEarthGravity );

    // Set half Earth as central body for half Earth gravity.
    halfEarthGravitationalForceModel.setGravitationalBody( &halfEarth );

    // Create a pointer to a new RK4 integrator.
    RungeKutta4thOrderFixedStepsize rk4ForHalfEarthGravity;

    // Set an initial stepsize for the integrators.
    rk4ForHalfEarthGravity.setInitialStepsize( 30.0 );

    // Create numerical propagator object for full Earth gravity.
    CartesianStateNumericalPropagator cartesianStateNumericalPropagatorForHalfEarthGravity;

    // Set Cartesian state numerical propagator as state derivative class for
    // RK4 integrator.
    rk4ForHalfEarthGravity.setObjectContainingStateDerivative(
            &cartesianStateNumericalPropagatorForHalfEarthGravity );

    // Set the integrator to use RK4.
    cartesianStateNumericalPropagatorForHalfEarthGravity.setIntegrator( &rk4ForHalfEarthGravity );

    // Add Asterix as the body that has to be propagated.
    cartesianStateNumericalPropagatorForHalfEarthGravity.addBody( &asterixForHalfEarthGravity );

    // Add half Earth gravitational force models as forces acting on Asterix.
    cartesianStateNumericalPropagatorForHalfEarthGravity.addForceModel(
                &asterixForHalfEarthGravity, &halfEarthGravitationalForceModel );
    cartesianStateNumericalPropagatorForHalfEarthGravity.addForceModel(
                &asterixForHalfEarthGravity, &halfEarthGravitationalForceModel );

    // Create series propagator for half Earth gravity.
    SeriesPropagator seriesPropagatorForHalfEarthGravity;

    // Set series propagation start time.
    seriesPropagatorForHalfEarthGravity.setSeriesPropagationStart( 0.0 );

    // Set the propagation end time.
    seriesPropagatorForHalfEarthGravity.setSeriesPropagationEnd( 86400.0 );

    // Set fixed output interval for series propagation
    seriesPropagatorForHalfEarthGravity.setFixedOutputInterval( 3600.0 );

    // Set propagator for series propagation.
    seriesPropagatorForHalfEarthGravity.setPropagator(
                &cartesianStateNumericalPropagatorForHalfEarthGravity );

    // Set initial state of Asterix for series propagation.
    seriesPropagatorForHalfEarthGravity.setInitialState( &asterixForHalfEarthGravity,
                                                         &stateOfAsterixForHalfEarthGravity );

    // Run simulation for full Earth.
    seriesPropagatorForHalfEarthGravity.execute( );

    // Get propagation history of Asterix for half Earth gravity case.
    asterixPropagationHistoryHalfGravity = seriesPropagatorForHalfEarthGravity
              .getPropagationHistoryAtFixedOutputIntervals( &asterixForHalfEarthGravity );

    // Check if results of full and half Earth gravity cases match.
    for ( unsigned int i = 0;
          i < ( seriesPropagatorForHalfEarthGravity.getSeriesPropagationEnd( )
                / seriesPropagatorForHalfEarthGravity.getFixedOutputInterval( ) ); i++ )
    {
        if ( asterixPropagationHistoryFullGravity[
                i * seriesPropagatorForFullEarthGravity.getFixedOutputInterval( ) ].state
             != asterixPropagationHistoryHalfGravity[
                i * seriesPropagatorForHalfEarthGravity.getFixedOutputInterval( ) ].state )
        {
            isNumericalPropagatorErroneous = true;

            std::cerr << "The numerical propagator does not produce consistent results, as "
                      << "running a simulation with full Earth gravity does not yield the "
                      << "same results as running a simulation with twice half Earth gravity."
                      << std::endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isNumericalPropagatorErroneous )
    {
        std::cerr << "testNumericalPropagator failed!" << std::endl;
    }

    return isNumericalPropagatorErroneous;
}

// End of file.
