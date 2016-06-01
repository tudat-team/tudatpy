/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      110622    F.M. Engelen      File created.
 *      110822    D. Dirkx          Removed no longer necessary unit tests.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      120328    D. Dirkx          Boostified unit tests.
 *      120405    K. Kumar          Ensured no interference between unit tests by placing them in
 *                                  local scope.
 *      121020    D. Dirkx          Update to new acceleration model architecture.
 *
 *    References
 *
 *    Notes
 *      The unit tests here are based off of expected values that are internally computed. Ideally,
 *      these should be based off of published values that can be found in literature. This will
 *      have to be updated in future for the code to be considered completely tested. In addition,
 *      more test values are required, as more the unit tests are benchmarked off of one set of
 *      data.
 *
 *      The class objects aerodynamicCoefficientInterface, aerodynamicForce, and aerodynamicMoment
 *      are declared multiple times in local scope currently since their member variables aren't
 *      set at construction, but rather through set-functions. Once these are adapted to be set
 *      through the constructor, const class objects can be declared that can be shared between the
 *      unit tests.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicRotationalAcceleration.h"
#include "Tudat/SimulationSetup/createFlightConditions.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

using namespace aerodynamics;
using namespace simulation_setup;

//! Test implementation of aerodynamic force and acceleration models.
BOOST_AUTO_TEST_CASE( testAerodynamicForceAndAcceleration )
{
    // Set force coefficients.
    const Eigen::Vector3d forceCoefficients( 1.1, 1.2, 1.3 );

    // Set dynamical model parameters.
    const double density = 3.5e-5;
    const double airSpeed = 3.491e3;
    const double dynamicPressure = 0.5 * density * airSpeed * airSpeed;
    const double referenceArea = 2.2;
    const double referenceLength = 3.2;
    const double mass = 1.93;

    // Compute expected force.
    const Eigen::Vector3d expectedForce = forceCoefficients * dynamicPressure * referenceArea;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Test 1: test the force model implemented as free function with primitive arguments.
    {
        // Compute force.
        Eigen::Vector3d force = computeAerodynamicForce( dynamicPressure,
                                                         referenceArea, forceCoefficients );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 2: test the force model implemented as free function, with coefficient interface
    //         argument.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
                createConstantCoefficientAerodynamicCoefficientInterface(
                    forceCoefficients, Eigen::Vector3d::Zero( ),
                    referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic force using free function with coefficient interface argument.
        Eigen::Vector3d force = computeAerodynamicForce( dynamicPressure,
                                                         aerodynamicCoefficientInterface );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );

        // Test aerodynamic coefficient interface properties
        BOOST_CHECK_EQUAL(
                    aerodynamicCoefficientInterface->getIndependentVariableNames( ).size( ), 0 );

        bool isVariableIndexTooHigh = 0;
        try
        {
            aerodynamicCoefficientInterface->getIndependentVariableName( 0 );
        }
        catch ( std::runtime_error )
        {
            isVariableIndexTooHigh = 1;
        }
        BOOST_CHECK_EQUAL( isVariableIndexTooHigh, 1 );
    }

    // Test 3: test the acceleration model implemented as free function with primitive arguments,
    //         based on the force that can be derived from the computed acceleration.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
                createConstantCoefficientAerodynamicCoefficientInterface(
                    forceCoefficients, Eigen::Vector3d::Zero( ),
                    referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic force from aerodynamic acceleration free function with primitive
        // arguments.
        Eigen::Vector3d force = computeAerodynamicAcceleration(
                    dynamicPressure, aerodynamicCoefficientInterface, mass ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 4: test the acceleration model implemented as free function with coefficient interface
    //         argument, based on the force that can be derived from the computed acceleration.
    {
        // Compute aerodynamic force from aerodynamic acceleration free function with
        // coefficient interface argument.
        Eigen::Vector3d force = computeAerodynamicAcceleration(
                    dynamicPressure, referenceArea, forceCoefficients, mass ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 5: Test the acceleration model class without inverted coefficients.
    {
        // Create aaerodynamic acceleration model class, no inverted coefficients, direct mass
        // and reference area.
        AerodynamicAccelerationPointer accelerationClass
                = boost::make_shared< AerodynamicAcceleration >(
                    boost::lambda::constant( forceCoefficients ),
                    boost::lambda::constant( density ),
                    boost::lambda::constant( airSpeed ),
                    mass, referenceArea, false );
        accelerationClass->updateMembers( );
        Eigen::Vector3d force = accelerationClass->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );

        // Create aerodynamic acceleration model class, no inverted coefficients, mass and
        // reference area set through boost::functions.
        AerodynamicAccelerationPointer accelerationClass2 =
                boost::make_shared< AerodynamicAcceleration >(
                    boost::lambda::constant( forceCoefficients ),
                    boost::lambda::constant( density ),
                    boost::lambda::constant( airSpeed ),
                    boost::lambda::constant( mass ),
                    boost::lambda::constant( referenceArea ),
                    false );
        accelerationClass2->updateMembers( );
        force = accelerationClass2->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 6: Test the acceleration model class with inverted coefficients
    {
        // Create aaerodynamic acceleration model class, inverted coefficients, direct mass
        // and reference area.
        AerodynamicAccelerationPointer accelerationClass =
                boost::make_shared< AerodynamicAcceleration >(
                    boost::lambda::constant( -forceCoefficients ),
                    boost::lambda::constant( density ),
                    boost::lambda::constant( airSpeed ),
                    mass, referenceArea, true );
        accelerationClass->updateMembers( );
        Eigen::Vector3d force = accelerationClass->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );

        // Create aerodynamic acceleration model class, inverted coefficients, mass and
        // reference area set through boost::functions.
        AerodynamicAccelerationPointer accelerationClass2 =
                boost::make_shared< AerodynamicAcceleration >(
                    boost::lambda::constant( -forceCoefficients ),
                    boost::lambda::constant( density ),
                    boost::lambda::constant( airSpeed ),
                    boost::lambda::constant( mass ),
                    boost::lambda::constant( referenceArea ),
                    true );
        accelerationClass2->updateMembers( );
        force = accelerationClass2->getAcceleration( ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }
}

//! Test implementation of aerodynamic moment and rotational acceleration models.
BOOST_AUTO_TEST_CASE( testAerodynamicMomentAndRotationalAcceleration )
{
    // Set moment coefficients.
    const Eigen::Vector3d momentCoefficients( -3.2, 1.0, 8.4 );

    // Set dynamical model parameters.
    const double dynamicPressure = 123.6;
    const double referenceArea = 1.7;
    const double referenceLength = 2.6;
    const double mass = 12.46;

    // Calculate expected moment.
    const Eigen::Vector3d expectedMoment = dynamicPressure * referenceArea *
            referenceLength * momentCoefficients;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Test 1: test the moment model implemented as free function with primitive arguments.
    {
        // Compute aerodynamic moment using free function with primitive arguments.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure, referenceArea,
                                                           referenceLength, momentCoefficients );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 2: test the moment moment implemented as free function with coefficient interface
    //         argument.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
        createConstantCoefficientAerodynamicCoefficientInterface(
            Eigen::Vector3d::Zero( ), momentCoefficients,
            referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic moment using free function with coefficient interface argument.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure,
                                                           aerodynamicCoefficientInterface );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 3: test the rotational acceleration model implemented as free function with coefficient
    //         interface argument, based on the force that can be derived from the computed
    //         acceleration.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterfacePointer aerodynamicCoefficientInterface =
                createConstantCoefficientAerodynamicCoefficientInterface(
                    Eigen::Vector3d::Zero( ), momentCoefficients,
                    referenceLength, referenceArea, referenceLength, Eigen::Vector3d::Zero( ) );

        // Compute aerodynamic moment from aerodynamic rotational acceleration free function with
        // primitive arguments.
        Eigen::Vector3d moment = computeAerodynamicRotationalAcceleration(
                    dynamicPressure, aerodynamicCoefficientInterface, mass ) * mass;

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 4: test the rotational acceleration model implemented as free function with primitive
    //         arguments, based on the force that can be derived from the computed acceleration.
    {
        // Compute aerodynamic moment from aerodynamic rotational acceleration free function with
        // primitive arguments.
        Eigen::Vector3d moment = computeAerodynamicRotationalAcceleration(
                    dynamicPressure, referenceArea, referenceLength,
                    momentCoefficients, mass ) * mass;

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
