/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *    Hackworth, M. Magnetic Fields and Magnetic Forces,
 *          http://www2.cose.isu.edu/~hackmart/Magnetic%20Fields_engphys.pdf, last accessed: 20th
 *          January, 2014.
 *    Khan, S. Magnetism 3 Whats happens when a speeding proton goes through a magnetic field, 
 *          https://www.khanacademy.org/science/physics/electricity-and-magnetism/v/magnetism-3,
 *          2013, last accessed: 16th January, 2014.
 *    Learning About Electronics Magnetic Force Calculator, 
 *          http://www.learningaboutelectronics.com/Articles/Magnetic-force-calculator.php#answer1,
 *          last accessed: 20th January, 2014.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/lorentzStaticMagneticForce.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/lorentzStaticMagneticAcceleration.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_lorentz_static_magnetic_acceleration_and_force_models )

//! Test implementation of static magnetic force model.
BOOST_AUTO_TEST_CASE( testStaticMagneticForceKhanAcademy )
{
    // Benchmark data is obtained using online video (Khan, 2013).

    // Tests approximate force on a proton due to a static uniform magnetic field. 

    // Set expected static magnetic force [N].
    const Eigen::Vector3d expectedStaticMagneticForce 
                = Eigen::Vector3d( 0.0, -4.80653199e-12, 0.0 );

    // Set velocity of accelerated body [m/s].
    const Eigen::Vector3d velocityOfBodySubjectToAcceleration = Eigen::Vector3d( 6.0e7, 0.0, 0.0 );

    // Set strength of local magnetic field [T].
    const Eigen::Vector3d localMagneticField = Eigen::Vector3d( 0.0, 0.0, 0.5 );

    // Set charge of accelerated body [C].
    const double chargeOfBodySubjectToAcceleration = 1.60217733e-19;

    // Compute static magnetic force [N].
    const Eigen::Vector3d computedStaticMagneticForce
                = electro_magnetism::computeLorentzForceDueToStaticMagneticField(
                    velocityOfBodySubjectToAcceleration,
                    localMagneticField,
                    chargeOfBodySubjectToAcceleration );

    // Compare computed and expected static magnetic force vectors.
    BOOST_CHECK_SMALL( computedStaticMagneticForce.x( ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_CLOSE_FRACTION( computedStaticMagneticForce.y( ),
                                expectedStaticMagneticForce.y( ),
                                std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_SMALL( computedStaticMagneticForce.z( ),
                       std::numeric_limits< double >::min( ) );
}


//! Test implementation of static magnetic force model.
BOOST_AUTO_TEST_CASE( testStaticMagneticForceArbitraryCharge )
{
    // Benchmark data is obtained using Learning About Electronics.

    // Tests force using arbitrary values for parameters.

    // Set expected static magnetic force [N].
    const Eigen::Vector3d expectedStaticMagneticForce = Eigen::Vector3d( 218.0, -172.0, 4.0 );

    // Set velocity of accelerated body [m/s].
    const Eigen::Vector3d velocityOfBodySubjectToAcceleration = Eigen::Vector3d( 4.0, 5.0, -3.0 );

    // Set strength of local magnetic field [T].
    const Eigen::Vector3d localMagneticField = Eigen::Vector3d( 2.0, 3.0, 20.0 );

    // Set charge of accelerated body [C].
    const double chargeOfBodySubjectToAcceleration = 2.0;

    // Compute static magnetic force [N].
    const Eigen::Vector3d computedStaticMagneticForce
                = electro_magnetism::computeLorentzForceDueToStaticMagneticField(
                    velocityOfBodySubjectToAcceleration,
                    localMagneticField,
                    chargeOfBodySubjectToAcceleration );

    // Compare computed and expected static magnetic force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedStaticMagneticForce,
                                       expectedStaticMagneticForce,
                                       std::numeric_limits< double >::epsilon( ) );
}

//! Test implementation of static magnetic force model.
BOOST_AUTO_TEST_CASE( testStaticMagneticForceHackworth )
{
    // Benchmark data is obtained from pdf (Hackworth, M.).
    // Tests force on charge in uniform magnetic field.

    // Set expected static magnetic force [N].
    const Eigen::Vector3d expectedStaticMagneticForce = Eigen::Vector3d( 77.0, 28.0, -10.0 );

    // Set velocity of accelerated body [m/s].
    const Eigen::Vector3d velocityOfBodySubjectToAcceleration = Eigen::Vector3d( -2.0, 3.0, -7.0 );

    // Set strength of local magnetic field [T].
    const Eigen::Vector3d localMagneticField = Eigen::Vector3d( 4.0, -11.0, 0.0 );

    // Set charge of accelerated body [C].
    const double chargeOfBodySubjectToAcceleration = -1.0;

    // Compute static magnetic force [N].
    const Eigen::Vector3d computedStaticMagneticForce
                = electro_magnetism::computeLorentzForceDueToStaticMagneticField(
                    velocityOfBodySubjectToAcceleration,
                    localMagneticField,
                    chargeOfBodySubjectToAcceleration );

    // Compare computed and expected static magnetic force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedStaticMagneticForce,
                                       expectedStaticMagneticForce,
                                       std::numeric_limits< double >::epsilon( ) );
}

//! Test implementation of static magnetic acceleration model.
BOOST_AUTO_TEST_CASE( testStaticMagneticAccelerationKhanAcademy )
{
    // Benchmark data is obtained using online video (Khan, 2013).

    // Tests approximate force on a proton due to a static uniform magnetic field. 

    // Set expected static magnetic force [N].
    const Eigen::Vector3d expectedStaticMagneticAcceleration =
            Eigen::Vector3d( 0.0, -2.8736514e15, 0.0 );

    // Set velocity of accelerated body [m/s].
    const Eigen::Vector3d velocityOfBodySubjectToAcceleration = Eigen::Vector3d( 6.0e7, 0.0, 0.0 );

    // Set strength of local magnetic field [T].
    const Eigen::Vector3d localMagneticField = Eigen::Vector3d( 0.0, 0.0, 0.5 );

    // Set charge of accelerated body [C].
    const double chargeOfBodySubjectToAcceleration = 1.60217733e-19;

    // Set mass of accelerated body [kg].
    const double massOfBodySubjectToAcceleration = 1.67262178e-27;

    // Compute static magnetic acceleration [N].
    const Eigen::Vector3d computedStaticMagneticAcceleration
                = electro_magnetism::computeLorentzAccelerationDueToStaticMagneticField(
                    velocityOfBodySubjectToAcceleration,
                    localMagneticField,
                    chargeOfBodySubjectToAcceleration,
                    massOfBodySubjectToAcceleration );

    // Compare computed and expected static magnetic acceleration vectors.
    BOOST_CHECK_SMALL( computedStaticMagneticAcceleration.x( ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_CLOSE_FRACTION( computedStaticMagneticAcceleration.y( ),
                                expectedStaticMagneticAcceleration.y( ),
                                1.0e-7 );

    BOOST_CHECK_SMALL( computedStaticMagneticAcceleration.z( ),
                       std::numeric_limits< double >::min( ) );
}

//! Test implementation of static magnetic acceleration model.
BOOST_AUTO_TEST_CASE( testStaticMagneticAccelerationArbitraryCharge )
{
    // Benchmark data is obtained using Learning About Electronics.

    // Tests force using arbitrary values for parameters.

    // Set expected static magnetic force [N].
    const Eigen::Vector3d expectedStaticMagneticAcceleration 
                = Eigen::Vector3d( 43.6, -34.4, 0.8 );

    // Set velocity of accelerated body [m/s].
    const Eigen::Vector3d velocityOfBodySubjectToAcceleration = Eigen::Vector3d( 4.0, 5.0, -3.0 );

    // Set strength of local magnetic field [T].
    const Eigen::Vector3d localMagneticField = Eigen::Vector3d( 2.0, 3.0, 20.0 );

    // Set charge of accelerated body [C].
    const double chargeOfBodySubjectToAcceleration = 2.0;

    // Set mass of accelerated body [kg].
    const double massOfBodySubjectToAcceleration = 5.0;

    // Compute static magnetic acceleration [N].
    const Eigen::Vector3d computedStaticMagneticAcceleration
                = electro_magnetism::computeLorentzAccelerationDueToStaticMagneticField(
                    velocityOfBodySubjectToAcceleration,
                    localMagneticField,
                    chargeOfBodySubjectToAcceleration,
                    massOfBodySubjectToAcceleration );

    // Compare computed and expected static magnetic acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedStaticMagneticAcceleration,
                                       expectedStaticMagneticAcceleration,
                                       std::numeric_limits< double >::epsilon( ) );
}




BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
