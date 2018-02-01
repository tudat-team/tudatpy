/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_KEPLERPROPAGATORTESTDATA_H
#define TUDAT_KEPLERPROPAGATORTESTDATA_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

using namespace orbital_element_conversions;

//! Typedef for propagation history.
typedef std::map < double, Eigen::Vector6d > PropagationHistory;

//! Get Earth gravitational parameter for benchmark data from (Melman, 2010).
double getMelmanEarthGravitationalParameter( )
{
    // Return Earth gravitational parameter [m^3 s^-2].
    return 3.986004415e14;
}

//! Get benchmark data from (Melman, 2010).
PropagationHistory getMelmanBenchmarkData( )
{
    // Declare benchmark pragation history.
    PropagationHistory benchmarkPropagationHistory;

    // Populate benchmark propagation history.
    Eigen::Vector6d stateInCartesianElements;
    Eigen::Vector6d stateInKeplerianElements;

    stateInCartesianElements << 6.75e6, 0.0, 0.0, 0.0, 8.0595973215e3, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getMelmanEarthGravitationalParameter( ) );
    benchmarkPropagationHistory[ 0.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -6.1318272067e6, 5.1974105627e6, 0.0,
            -4.7375063953e3, -4.8565484865e3, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getMelmanEarthGravitationalParameter( ) );
    benchmarkPropagationHistory[ 86400.0 ] = stateInKeplerianElements;

    // Return benchmark data.
    return benchmarkPropagationHistory;
}

//! Get ODTBX benchmark data.
PropagationHistory getODTBXBenchmarkData( )
{
    // Declare benchmark pragation history.
    PropagationHistory benchmarkPropagationHistory;

    // Populate benchmark propagation history.
    Eigen::Vector6d stateInKeplerianElements;

    stateInKeplerianElements << 42165.3431351313e3, 0.26248354351331, 0.30281462522101,
            4.71463172847351, 4.85569272927819, 2.37248926702153;

    // Set time step.
    double timeStep = 8640.0;

    for ( unsigned int i = 0; i < 11; i++ )
    {
        benchmarkPropagationHistory[ static_cast< double >( i ) * timeStep ]
                = stateInKeplerianElements;
    }

    benchmarkPropagationHistory[  1.0 * timeStep ]( 5 ) = 2.79722436211144;
    benchmarkPropagationHistory[  2.0 * timeStep ]( 5 ) = 3.18337407409023;
    benchmarkPropagationHistory[  3.0 * timeStep ]( 5 ) = 3.57400974200765;
    benchmarkPropagationHistory[  4.0 * timeStep ]( 5 ) = 4.01425565759545;
    benchmarkPropagationHistory[  5.0 * timeStep ]( 5 ) = 4.57232665706546;
    benchmarkPropagationHistory[  6.0 * timeStep ]( 5 ) = 5.35956850972672;
    benchmarkPropagationHistory[  7.0 * timeStep ]( 5 ) = 0.137251905665217;
    benchmarkPropagationHistory[  8.0 * timeStep ]( 5 ) = 1.14521863765007;
    benchmarkPropagationHistory[  9.0 * timeStep ]( 5 ) = 1.86433634881636;
    benchmarkPropagationHistory[ 10.0 * timeStep ]( 5 ) = 2.38486787064101;

    return benchmarkPropagationHistory;
}

//! Get GTOP gravitational parameter for benchmark data from GTOP.
double getGTOPGravitationalParameter( )
{
    // Return inaccurate gravitational parameter of the Sun [m^3 s^-2].
    return 1.327e20;
}

//! Get benchmark data from GTOP.
PropagationHistory getGTOPBenchmarkData( )
{
    // Declare benchmark pragation history.
    PropagationHistory benchmarkPropagationHistory;

    // Populate benchmark propagation history. Obtained by propagating an orbit starting at
    // x = 1.5e11, V_y = 6.0e4 for a period of 100 days four times consecutively. Since GTOP does
    // not work with similar Keplerian orbital elements, the cartesian elements resulting from that
    // are converted to Keplerian elements first. These initial starting coordinates correspond to
    // a semi major axis of -7.24873e+010 meters and an eccentricity of 3.06933.
    Eigen::Vector6d stateInCartesianElements, stateInKeplerianElements;

    stateInCartesianElements << 1.5e11, 0.0, 0.0, 0.0, 6.0e4, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 0.0 ] = stateInKeplerianElements;

    stateInCartesianElements << 50369576778.98602, 453006898372.5074, 2.156946592732799e-005,
            -14654.13750690802, 46884.94068619227, 4.665334803219454e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -76810236076.38216, 842661848023.4473, 6.100268297443444e-005,
            -14683.57015580225, 43917.12010513522, 4.48721854707566e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 * 2.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -203052258817.9893, 1216808019495.603, 9.937145023346651e-005,
            -14543.34378775917, 42828.66589049961, 4.403399939593385e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 * 3.0 ] = stateInKeplerianElements;

    stateInCartesianElements << -328225472457.8796, 1584186440047.591, 0.0001371949389119038,
            -14437.813524927732, 42264.20425643964, 4.355914471377053e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
                stateInCartesianElements, getGTOPGravitationalParameter( ) );
    benchmarkPropagationHistory[ 8640000.0 * 4.0 ] = stateInKeplerianElements;

    // Return benchmark data.
    return benchmarkPropagationHistory;
}

} // namespace unit_tests
} // namespace tudat
#endif // TUDAT_KEPLERPROPAGATORTESTDATA_H
