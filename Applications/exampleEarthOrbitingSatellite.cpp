/*! \file exampleEarthOrbitingSatellite.cpp
 *    Source file that contains an example application of an Earth-orbiting
 *    satellite mission scenario using Tudat.
 *
 *    Path              : /Applications/
 *    Version           : 6
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
 *    Last modified     : 15 August, 2011
 *
 *    References
 *
 *    Notes
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
 *      101111    K. Kumar          File created.
 *      110113    K. Kumar          Scenario updated to use latest version of code; added file
 *                                  header and footer.
 *      110202    K. Kumar          Scenario updated to use latest version of code.
 *      110216    K. Kumar          Migrated to applications namespace.
 *      110217    K. Kumar          Function name changed.
 *      110815    K. Kumar          Updated with mass of Asterix.
 *      111024    K. Kumar          Modified to be executable program with main-function as
 *                                  suggested by M. Persson.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/Bodies/planet.h"
#include "Astrodynamics/Bodies/vehicle.h"
#include "Astrodynamics/ForceModels/gravitationalForceModel.h"
#include "Astrodynamics/Propagators/cartesianStateNumericalPropagator.h"
#include "Astrodynamics/Propagators/seriesPropagator.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Mathematics/unitConversions.h"
#include "Mathematics/NumericalIntegrators/rungeKutta4thOrderFixedStepsize.h"
#include "Output/writingOutputToFile.h"

#include "Astrodynamics/astrodynamicsFunctions.h"

//! Execute example of an Earth-orbiting satellite.
int main( )
{
    // Using declarations.
    using std::cout;
    using std::endl;
    using namespace tudat;

    // Create the state of Asterix given in Cartesian elements.
    CartesianElements stateOfAsterix;

    // Fill initial state vector with position and
    // velocity given for Asterix.
    // Position is given in kilometers and
    // velocity is given in kilometers per second.
    stateOfAsterix.setCartesianElementX( 7000.0 );
    stateOfAsterix.setCartesianElementY( 0.0 );
    stateOfAsterix.setCartesianElementZ( 0.0 );
    stateOfAsterix.setCartesianElementXDot( 0.0 );
    stateOfAsterix.setCartesianElementYDot( 5.0 );
    stateOfAsterix.setCartesianElementZDot( 7.0 );

    // Convert initial state vector to meters from
    // kilometers.
    stateOfAsterix.state = unit_conversions::convertKilometersToMeters( stateOfAsterix.state );

    // Create map of propagation history of Asterix.
    std::map < double, State > asterixPropagationHistory;

    // Create a new vehicle object for Asterix.
    Vehicle asterix;

    // Set mass of Asterix.
    asterix.setMass( 1000.0 );

    // Create pre-defined Earth object.
    Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( Planet::earth );

    // Create a new gravitational force model for Earth.
    GravitationalForceModel earthGravitiationalForceModel;

    // Set Asterix as body subject to Earth gravitational force.
    earthGravitiationalForceModel.setBodySubjectToForce( &asterix );

    // Set Earth as central body for gravity.
    earthGravitiationalForceModel.setGravitationalBody( &predefinedEarth );

    // Create a new RK4 integrator object.
    RungeKutta4thOrderFixedStepsize rungeKutta4;

    // Set an initial stepsize for our integrator.
    rungeKutta4.setInitialStepsize( 30.0 );

    // Create Cartesian state numerical propagator object.
    CartesianStateNumericalPropagator cartesianStateNumericalPropagator;

    // Set object containing state derivative function for integrator. In this
    // case the Cartesian state numerical propagator contains the state
    // derivative function.
    rungeKutta4.setObjectContainingStateDerivative( &cartesianStateNumericalPropagator );

    // Set the integrator to use RK4.
    cartesianStateNumericalPropagator.setIntegrator( &rungeKutta4 );

    // Add Asterix as the body that has to be propagated.
    cartesianStateNumericalPropagator.addBody( &asterix );

    // Add Earth gravity as force acting on Asterix.
    cartesianStateNumericalPropagator.addForceModel( &asterix, &earthGravitiationalForceModel );

    // Create series propagator object to propagate timeseries.
    SeriesPropagator seriesPropagator;

    // Set Cartesian state numerical propagator for timeseries propagation.
    seriesPropagator.setPropagator( &cartesianStateNumericalPropagator );

    // Set start of the timeseries for propagation.
    seriesPropagator.setSeriesPropagationStart( 0.0 );

    // Set end of the timeseries for propagation.
    seriesPropagator.setSeriesPropagationEnd( 86400.0 );

    // Set fixed output interval for timeseries propagation
    seriesPropagator.setFixedOutputInterval( 60.0 );

    // Set initial state of Asterix for timeseries propagation.
    seriesPropagator.setInitialState( &asterix, &stateOfAsterix );

    // Run timeseries propagation.
    seriesPropagator.execute( );

    // Get propagation history of Asterix.
    asterixPropagationHistory
            = seriesPropagator.getPropagationHistoryAtFixedOutputIntervals( &asterix );

    // Output final state vector of Asterix to screen.
    cout << "Asterix final state in km(/s):" << endl;
    cout << unit_conversions::convertMetersToKilometers(
            asterixPropagationHistory[ seriesPropagator.getSeriesPropagationEnd( ) ].state )
         << endl;

    // Write propagation history of Asterix to file.
    WritingOutputToFile fileWriter;
    fileWriter.writePropagationHistoryToFile( asterixPropagationHistory,
                                              "asterixExampleEarthOrbitingSatellite.dat" );

    return 0;
}

// End of file.
