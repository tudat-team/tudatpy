/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      120611    P. Musegaas       First creation of code.
 *      120827    P. Musegaas       Adaptation to own ephemeris type.
 *
 *    References
 *      Musegaas, P. (2012). Optimization of Space Trajectories Including Multiple Gravity Assists
 *          and Deep Space Maneuvers. MSc Thesis, Delft University of Technology, Delft,
 *          The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <vector>

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>


#include <Tudat/InputOutput/basicInputOutput.h>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"

namespace tudat
{
namespace unit_tests
{

// Additional namespaces to be used.
using namespace tudat::input_output;
using namespace tudat::input_output::parsed_data_vector_utilities;
using namespace tudat::transfer_trajectories;

//! Test implementation of trajectory class
BOOST_AUTO_TEST_SUITE( test_trajectory )

//! Test delta-V computation for simple MGA trajectory model.
BOOST_AUTO_TEST_CASE( testMGATrajectory )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 15.0e-2;

    // Expected test result based on the ideal Cassini 1 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 4930.72686847243;

    // Specify required parameters
    // Specify the number of legs and type of legs.
    const int numberOfLegs = 6;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[ 0 ] = mga_Departure; legTypeVector[ 1 ] = mga_Swingby;
    legTypeVector[ 2 ] = mga_Swingby; legTypeVector[ 3 ] = mga_Swingby; legTypeVector[ 4 ] = mga_Swingby;
    legTypeVector[ 5 ] = capture;

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer >
            ephemerisVector( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
    ephemerisVector[ 5 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );

    // Create gravitational parameter vector
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.24860e14, 3.24860e14, 3.9860119e14, 1.267e17,
            3.79e16;

    // Create variable vector.
    Eigen::VectorXd variableVector ( numberOfLegs + 1 );
    variableVector << -789.8117, 158.302027105278, 449.385873819743, 54.7489684339665,
                      1024.36205846918, 4552.30796805542, 1/*dummy*/;
    variableVector *= physical_constants::JULIAN_DAY;

    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes ( 2 ), eccentricities ( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities << 0., 0.98;

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector
    Eigen::VectorXd minimumPericenterRadii ( numberOfLegs );
    minimumPericenterRadii << 6778000., 6351800., 6351800., 6778000., 600000000., 600000000.;

    // Create the trajectory problem.
    Trajectory Cassini1 ( numberOfLegs, legTypeVector, ephemerisVector,
                          gravitationalParameterVector, variableVector, sunGravitationalParameter,
                          minimumPericenterRadii, semiMajorAxes, eccentricities );

    // Start the deltaV vector.
    double resultingDeltaV;
    Cassini1.calculateTrajectory( resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test delta-V computation for MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory1 )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 2.0e-3;

    // Expected test result based on the ideal Messenger trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8630.83256199051;

    // Specify required parameters
    // Specify the number of legs and type of legs.
    const int numberOfLegs = 5;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[0] = mga1DsmVelocity_Departure; legTypeVector[1] = mga1DsmVelocity_Swingby;
    legTypeVector[2] = mga1DsmVelocity_Swingby; legTypeVector[3] = mga1DsmVelocity_Swingby;
    legTypeVector[4] = capture;

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer >
            ephemerisVector( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );

    // Create gravitational parameter vector
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.9860119e14, 3.24860e14, 3.24860e14, 2.2321e13;

    // Create variable vector.
    Eigen::VectorXd variableVector ( numberOfLegs /*time of flight*/ + 1 /*start epoch*/ +
                                     4 * ( numberOfLegs - 1 ) /*additional variables for model,
                                                                except the final capture leg*/ );

    // Add the time of flight and start epoch, which are in JD.
    variableVector << 1171.64503236 * physical_constants::JULIAN_DAY,
                      399.999999715 * physical_constants::JULIAN_DAY,
                      178.372255301 * physical_constants::JULIAN_DAY,
                      299.223139512 * physical_constants::JULIAN_DAY,
                      180.510754824 * physical_constants::JULIAN_DAY,
                      1, // The capture time is irrelevant for the final leg.
    // Add the additional variables.
                      0.234594654679, 1408.99421278, 0.37992647165 * 2 * 3.14159265358979,
                      std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2, // 1st leg.
                      0.0964769387134, 1.35077257078, 1.80629232251 * 6.378e6, 0.0, // 2nd leg.
                      0.829948744508, 1.09554368115, 3.04129845698 * 6.052e6, 0.0, // 3rd leg.
                      0.317174785637, 1.34317576594, 1.10000000891 * 6.052e6, 0.0; //4th leg.

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector
    Eigen::VectorXd minimumPericenterRadii ( numberOfLegs );
    minimumPericenterRadii << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN;

    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes ( 2 ), eccentricities ( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ),
                     std::numeric_limits< double >::infinity( );
    eccentricities << 0., 0.;

    // Create the trajectory problem.
    Trajectory Messenger ( numberOfLegs, legTypeVector, ephemerisVector,
                           gravitationalParameterVector, variableVector, sunGravitationalParameter,
                           minimumPericenterRadii, semiMajorAxes, eccentricities );

    // Start the deltaV vector.
    double resultingDeltaV;
    Messenger.calculateTrajectory( resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test delta-V computation for another MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory2 )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 4.0e-2;

    // Expected test result based on the ideal Cassini 2 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8385.15784516116;

    // Specify required parameters
    // Specify the number of legs and type of legs.
    const int numberOfLegs = 6;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[0] = mga1DsmVelocity_Departure; legTypeVector[1] = mga1DsmVelocity_Swingby;
    legTypeVector[2] = mga1DsmVelocity_Swingby; legTypeVector[3] = mga1DsmVelocity_Swingby;
    legTypeVector[4] = mga1DsmVelocity_Swingby; legTypeVector[5] = capture;

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer >
            ephemerisVector( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
    ephemerisVector[ 5 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );
    // Create gravitational parameter vector
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.24860e14, 3.24860e14, 3.9860119e14, 1.267e17,
                                    3.79e16;

    // Create variable vector.
    Eigen::VectorXd variableVector ( numberOfLegs /*time of flight*/ + 1 /*start epoch*/ +
                                     4 * ( numberOfLegs - 1 ) /*additional variables for model,
                                                                except the final capture leg*/ );

    // Add the time of flight and start epoch, which are in JD.
    variableVector << -779.046753814506 * physical_constants::JULIAN_DAY,
                      167.378952534645 * physical_constants::JULIAN_DAY,
                      424.028254165204 * physical_constants::JULIAN_DAY,
                      53.2897409769205 * physical_constants::JULIAN_DAY,
                      589.766954923325 * physical_constants::JULIAN_DAY,
                      2200.00000000000 * physical_constants::JULIAN_DAY,
                      1, // The capture time is irrelevant for the final leg.
    // Add the additional variables.
            0.769483451363201, 3259.11446832345, 0.525976214695235 * 2 * 3.14159265358979,
            std::acos(  2 * 0.38086496458657 - 1 ) - 3.14159265358979 / 2, // 1st leg.
            0.513289529822621, -1.5937371121191, 1.34877968657176 * 6.052e6, 0.0, // 2nd leg.
            0.0274175362264024, -1.95952512232447, 1.05 * 6.052e6, 0.0, // 3rd leg.
            0.263985256705873, -1.55498859283059, 1.30730278372017 * 6.378e6, 0.0, //4th leg.
            0.599984695281461, -1.5134625299674, 69.8090142993495 * 7.1492e7, 0.0; //5th leg.

    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes ( 2 ), eccentricities ( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ),
                     std::numeric_limits< double >::infinity( );
    eccentricities << 0., 0.;

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector
    Eigen::VectorXd minimumPericenterRadii ( numberOfLegs );
    minimumPericenterRadii << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN;

    // Create the trajectory problem.
    Trajectory Cassini2 ( numberOfLegs, legTypeVector, ephemerisVector,
                          gravitationalParameterVector, variableVector, sunGravitationalParameter,
                          minimumPericenterRadii, semiMajorAxes, eccentricities );

    // Start the deltaV vector.
    double resultingDeltaV;
    Cassini2.calculateTrajectory( resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test update functionality for simple MGA trajectory model.
BOOST_AUTO_TEST_CASE( testUpdateMGATrajectory )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 15.0e-2;

    // Expected test result based on the ideal Cassini 1 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 4930.72686847243;

    // Specify required parameters
    // Specify the number of legs and type of legs.
    const int numberOfLegs = 6;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[0] = mga_Departure; legTypeVector[1] = mga_Swingby;
    legTypeVector[2] = mga_Swingby; legTypeVector[3] = mga_Swingby; legTypeVector[4] = mga_Swingby;
    legTypeVector[5] = capture;

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer >
            ephemerisVector( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
    ephemerisVector[ 5 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );

    // Create gravitational parameter vector
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.24860e14, 3.24860e14, 3.9860119e14, 1.267e17,
            3.79e16;

    // Create a dummy variable vector. This dummy variable vector is finally multiplied by 1.01 to
    // ensure that the code breaks if the update function is not complete. This is done in favour
    // of passing the NaNs, which was done in the mission legs tests, because this will generate
    // error messages caused by an impossible conversion betweeen mean and eccentric anomaly.
    Eigen::VectorXd dummyVariableVector ( numberOfLegs + 1 );
    dummyVariableVector << -789.8117, 158.302027105278, 449.385873819743, 54.7489684339665,
                           1024.36205846918, 4552.30796805542, 1/*dummy*/;
    dummyVariableVector *= physical_constants::JULIAN_DAY;
    dummyVariableVector *= 1.01;


    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes ( 2 ), eccentricities ( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities << 0., 0.98;

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector
    Eigen::VectorXd minimumPericenterRadii ( numberOfLegs );
    minimumPericenterRadii << 6778000., 6351800., 6351800., 6778000., 600000000., 600000000.;

    // Create the trajectory problem.
    Trajectory cassini1 ( numberOfLegs, legTypeVector, ephemerisVector,
                          gravitationalParameterVector, dummyVariableVector,
                          sunGravitationalParameter, minimumPericenterRadii, semiMajorAxes,
                          eccentricities );

    // Start the deltaV vector.
    double resultingDeltaV;
    cassini1.calculateTrajectory( resultingDeltaV );

    // Create variable vector containing the updated variables.
    Eigen::VectorXd variableVector ( numberOfLegs + 1 );
    variableVector << -789.8117, 158.302027105278, 449.385873819743, 54.7489684339665,
                      1024.36205846918, 4552.30796805542, 1/*dummy*/;
    variableVector *= physical_constants::JULIAN_DAY;

    // Pass the variable to the class.
    cassini1.updateVariableVector( variableVector );
    cassini1.updateEphemeris( );

    // Recalculate the trajectory.
    cassini1.calculateTrajectory( resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test update functionality for MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testUpdateMGA1DSMVFTrajectory )
{
    // Set tolerance. Due to ephemeris differences, higher accuracy cannot be achieved.
    const double tolerance = 2.0e-3;

    // Expected test result based on the ideal Messenger trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8630.83256199051;

    // Specify required parameters
    // Specify the number of legs and type of legs.
    const int numberOfLegs = 5;
    std::vector< int > legTypeVector;
    legTypeVector.resize( numberOfLegs );
    legTypeVector[0] = mga1DsmVelocity_Departure; legTypeVector[1] = mga1DsmVelocity_Swingby;
    legTypeVector[2] = mga1DsmVelocity_Swingby; legTypeVector[3] = mga1DsmVelocity_Swingby;
    legTypeVector[4] = capture;

    // Create the ephemeris vector.
    std::vector< ephemerides::EphemerisPointer >
            ephemerisVector( numberOfLegs );
    ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
    ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );

    // Create gravitational parameter vector
    Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
    gravitationalParameterVector << 3.9860119e14, 3.9860119e14, 3.24860e14, 3.24860e14, 2.2321e13;

    // Create a dummy variable vector. This dummy variable vector is finally multiplied by 1.01 to
    // ensure that the code breaks if the update function is not complete. This is done in favour
    // of passing the NaNs, which was done in the mission legs tests, because this will generate
    // error messages caused by an impossible conversion betweeen mean and eccentric anomaly.
    Eigen::VectorXd dummyVariableVector ( numberOfLegs /*time of flight*/ + 1 /*start epoch*/ +
                                          4 * ( numberOfLegs - 1 ) /*additional variables for
                                          model, except the final capture leg*/ );

    // Add the time of flight and start epoch, which are in JD.
    dummyVariableVector << 1171.64503236 * physical_constants::JULIAN_DAY,
                           399.999999715 * physical_constants::JULIAN_DAY,
                           178.372255301 * physical_constants::JULIAN_DAY,
                           299.223139512 * physical_constants::JULIAN_DAY,
                           180.510754824 * physical_constants::JULIAN_DAY,
                           1, // The capture time is irrelevant for the final leg.
    // Add the additional variables.
                           0.234594654679, 1408.99421278, 0.37992647165 * 2 * 3.14159265358979,
                           std::acos(  2 * 0.498004040298 - 1 ) - 3.14159265358979 / 2, // 1st leg.
                           0.0964769387134, 1.35077257078, 1.80629232251 * 6.378e6, 0.0,// 2nd leg.
                           0.829948744508, 1.09554368115, 3.04129845698 * 6.052e6, 0.0, // 3rd leg.
                           0.317174785637, 1.34317576594, 1.10000000891 * 6.052e6, 0.0; //4th leg.
    dummyVariableVector *= 1.01;

    // Sun gravitational parameter.
    const double sunGravitationalParameter = 1.32712428e20;

    // Create minimum pericenter radii vector.
    Eigen::VectorXd minimumPericenterRadii ( numberOfLegs );
    minimumPericenterRadii << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN;

    // Create departure and capture variables.
    Eigen::VectorXd semiMajorAxes ( 2 ), eccentricities ( 2 );
    semiMajorAxes << std::numeric_limits< double >::infinity( ),
                     std::numeric_limits< double >::infinity( );
    eccentricities << 0., 0.;

    // Create the trajectory problem.
    Trajectory messenger ( numberOfLegs, legTypeVector, ephemerisVector,
                           gravitationalParameterVector, dummyVariableVector,
                           sunGravitationalParameter, minimumPericenterRadii, semiMajorAxes,
                           eccentricities );

    // Start the deltaV vector.
    double resultingDeltaV;
    messenger.calculateTrajectory( resultingDeltaV );

    // Create variable vector containing the updated variables.
    Eigen::VectorXd variableVector ( numberOfLegs /*time of flight*/ + 1 /*start epoch*/ +
                                     4 * ( numberOfLegs - 1 ) /*additional variables for model,
                                                                except the final capture leg*/ );

    // Add the time of flight and start epoch, which are in JD.
    variableVector << 1171.64503236 * physical_constants::JULIAN_DAY,
                      399.999999715 * physical_constants::JULIAN_DAY,
                      178.372255301 * physical_constants::JULIAN_DAY,
                      299.223139512 * physical_constants::JULIAN_DAY,
                      180.510754824 * physical_constants::JULIAN_DAY,
                      1, // The capture time is irrelevant for the final leg.
    // Add the additional variables.
                      0.234594654679, 1408.99421278, 0.37992647165 * 2 * 3.14159265358979,
                      std::acos(  2 * 0.498004040298 - 1 ) - 3.14159265358979 / 2, // 1st leg.
                      0.0964769387134, 1.35077257078, 1.80629232251 * 6.378e6, 0.0, // 2nd leg.
                      0.829948744508, 1.09554368115, 3.04129845698 * 6.052e6, 0.0, // 3rd leg.
                      0.317174785637, 1.34317576594, 1.10000000891 * 6.052e6, 0.0; //4th leg.

    // Pass the variable to the class.
    messenger.updateVariableVector( variableVector );
    messenger.updateEphemeris( );

    // Recalculate the trajectory.
    messenger.calculateTrajectory( resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
