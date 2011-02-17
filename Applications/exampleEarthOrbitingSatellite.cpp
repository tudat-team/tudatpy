/*! \file exampleEarthOrbitingSatellite.cpp
 *    Source file that contains an example application of an Earth-orbiting
 *    satellite mission scenario using Tudat.
 *
 *    Path              : /Applications/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Date created      : 11 November, 2010
 *    Last modified     : 17 February, 2011
 *
 *    References
 *
 *    Notes
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
 *      101111    K. Kumar          File created.
 *      110113    K. Kumar          Scenario updated to use latest version of
 *                                  code; added file header and footer.
 *      110202    K. Kumar          Scenario updated to use latest version of
 *                                  code.
 *      110216    K. Kumar          Migrated to applications namespace.
 *      110217    K. Kumar          Function name changed.
 */

// Include statements.
#include "exampleEarthOrbitingSatellite.h"

// Using declarations.
using std::cout;
using std::endl;

//! Namespace for all applications.
namespace applications
{

//! Execute example of an Earth-orbiting satellite.
void executeEarthOrbitingSatelliteExample( )
{
    // Create pointer to the state of Asterix given in Cartesian elements.
    CartesianElements* pointerToStateOfAsterix = new CartesianElements;

    // Fill initial state vector with position and
    // velocity given for Asterix.
    // Position is given in kilometers and
    // velocity is given in kilometers per second.
    pointerToStateOfAsterix->setCartesianElementX( 7000.0 );
    pointerToStateOfAsterix->setCartesianElementY( 0.0 );
    pointerToStateOfAsterix->setCartesianElementZ( 0.0 );
    pointerToStateOfAsterix->setCartesianElementXDot( 0.0 );
    pointerToStateOfAsterix->setCartesianElementYDot( 5.0 );
    pointerToStateOfAsterix->setCartesianElementZDot( 7.0 );

    // Convert initial state vector to meters from
    // kilometers.
    pointerToStateOfAsterix->state =
            unit_conversions::convertKilometersToMeters(
                    pointerToStateOfAsterix->state );

    // Create map of propagation history of Asterix.
    std::map < double, State* > asterixPropagationHistory;

    // Create a pointer to new vehicle for Asterix.
    Vehicle* pointerToAsterix = new Vehicle;

    // Create pre-defined Earth object.
    CelestialBody* pointerToEarth = predefined_planets::
                                    createPredefinedPlanet(
                                            predefined_planets::earth );

    // Create a pointer to a new Gravity object for Earth.
    Gravity* pointerToEarthGravity = new Gravity;

    // Set Earth as central body for gravity.
    pointerToEarthGravity->setBody( pointerToEarth );

    // Create a pointer to a new RK4 integrator.
    RungeKutta4thOrderFixedStepsize* pointerToRK4
            = new RungeKutta4thOrderFixedStepsize;

    // Set an initial stepsize for our integrator.
    pointerToRK4->setInitialStepsize( 30.0 );

    // Create numerical propagator object.
    NumericalPropagator numericalPropagator;

    // Set fixed output interval for output in numerical
    // propagator object.
    numericalPropagator.setFixedOutputInterval( 60.0 );

    // Set the propagation start time.
    numericalPropagator.setPropagationIntervalStart( 0.0 );

    // Set the propagation end time.
    numericalPropagator.setPropagationIntervalEnd( 86400.0 );

    // Set the integrator to use RK4.
    numericalPropagator.setIntegrator( pointerToRK4 );

    // Add Asterix as the body that has to be propagated.
    numericalPropagator.addBody( pointerToAsterix );

    // Add Earth gravity as force acting on Asterix.
    numericalPropagator.addForceModel( pointerToAsterix,
                                       pointerToEarthGravity );

    // Set initial state of Asterix.
    numericalPropagator.setInitialState(
            pointerToAsterix, pointerToStateOfAsterix );

    // Run simulation.86400
    numericalPropagator.propagate( );

    // Get propagation history of Asterix.
    asterixPropagationHistory = numericalPropagator.
                                getPropagationHistoryAtFixedOutputIntervals(
                                        pointerToAsterix );

    // Output final state vector of Asterix to screen.
    cout << "Asterix final state in km(/s):" << endl;
    cout << unit_conversions::convertMetersToKilometers(
            numericalPropagator.getFinalState( pointerToAsterix )->state )
         << endl;

    // Write propagation history of Asterix to file.
    WritingOutputToFile::
            writePropagationHistoryToFile(
                    asterixPropagationHistory,
                    "AsterixExampleEarthOrbitingSatellite.dat" );

}

}

// End of file.
