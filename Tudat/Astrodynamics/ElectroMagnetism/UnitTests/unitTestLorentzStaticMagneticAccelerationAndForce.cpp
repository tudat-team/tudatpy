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
 *      130124    S. Hirsh          File created.
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
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <TudatCore/Basics/testMacros.h>

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
