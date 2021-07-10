/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Izzo, D. and Vinko, T. ACT - Informatics - GTOP Database, ESA Advanced Concept Team, last
 *          accessed on 2012-01-12. http://www.esa.int/gsp/ACT/inf/op/globopt.htm.
 *      Musegaas, P. Gravity Assist calculation Verification.xlsx, last accessed: 3 December 2012,
 *          http://tudat.tudelft.nl/projects/tudat/wiki/Unit_tests, 2012.
 *
 *    Notes
 *      Three main functions are tested in these unit tests.
 *        Regarding the deltaV calculation gravity assist method:
 *          There is a complicated if-statement in this method. Hence many unit test are performed
 *          to test the functionality. Also various limit cases failed previously, hence many tests
 *          for this are also included:
 *              Case 1: required bending angle > maximum bending angle:
 *                  Two tests were written. In the first one no velocity effect is needed. This
 *                  test has a low accuracy, which should be replaced one day (it still relies on
 *                  hand calculator calculations done in 2011). In the second one a combination of
 *                  bending-effect deltaV and velocity-effect deltaV is calculated. This test has
 *                  been calculated using Tudat, and was verified using Excel.
 *                  Could definitely be improved.
 *              Case 2: no assist is required:
 *                  One test was written.
 *              Case 3: velocity effect deltaV only, using eccentricity iteration scheme:
 *                  Four tests were written. The first one calculates a case from Cassini-1 of GTOP
 *                  with high precision. The other three test limit cases: low incoming, high
 *                  outgoing velocity; high incoming, low outgoing velocity; low incoming, low
 *                  outgoing velocity. These tests were calculated using Tudat, but verified in
 *                  Excel to be exactly correct.
 *              Case 4: velocity effect deltaV only, using pericenter radius iteration scheme:
 *                  The same four tests as for case 3 were used.
 *        Regarding the unpowered gravity assist propagator:
 *          One test was written, based on GTOP. This should be a satisfactory test.
 *        Regarding the powered gravity assist propagator:
 *          Two tests were written. The first one is similar to the unpowered gravity assist
 *          propagator. The second one is reverse engineered from the Cassini-1 test, similar to
 *          the one in the deltaV calculation test.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>

#include "tudat/astro/mission_segments/createTransferTrajectory.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

//! Test of gravity assist code.
BOOST_AUTO_TEST_SUITE( test_patched_conic )

//! Test bending angle Delta-V computation.
BOOST_AUTO_TEST_CASE( TestPatchedConic )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );

   // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
