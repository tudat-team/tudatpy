/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110622    F.M. Engelen      First creation of code.
 *      110822    D. Dirkx          Removed no longer necessary unit tests.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      120328    D. Dirkx          Boostified unit tests.
 *      120405    K. Kumar          Ensured no interference between unit tests by placing them in
 *                                  local scope.
 *
 *    References
 *
 *    The unit tests here are based off of expected values that are internally computed. Ideally,
 *    these should be based off of published values that can be found in literature. This will have
 *    to be updated in future for the code to be considered completely tested. In addition, more
 *    test values are required, as more the unit tests are benchmarked off of one set of data.
 *
 *    The class objects aerodynamicCoefficientInterface, aerodynamicForce, and aerodynamicMoment
 *    are declared multiple times in local scope currently since their member variables aren't set
 *    at construction, but rather through set-functions. Once these are adapted to be set through
 *    the constructor, const class objects can be declared that can be shared between the unit
 *    tests.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicForce.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicMoment.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicRotationalAcceleration.h"

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

//! Test implementation of aerodynamic force and acceleration models.
BOOST_AUTO_TEST_CASE( testAerodynamicForceAndAcceleration )
{
    // Using declarations.
    using tudat::AerodynamicCoefficientInterface;
    using tudat::astrodynamics::force_models::AerodynamicForce;
    using tudat::astrodynamics::states::State;
    using tudat::astrodynamics::force_models::computeAerodynamicForce;
    using tudat::astrodynamics::acceleration_models::computeAerodynamicAcceleration;

    // Set force coefficients.
    const Eigen::Vector3d forceCoefficients( 1.1, 1.2, 1.3 );

    // Set dynamical model parameters.
    const double dynamicPressure = 50.0;
    const double referenceArea = 2.2;
    const double referenceLength = 3.2;
    const double mass = 1.93;

    // Compute expected force.
    const Eigen::Vector3d expectedForce = forceCoefficients * dynamicPressure * referenceArea;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Test 1: test the force model implemented as class.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface
                = boost::make_shared< AerodynamicCoefficientInterface >( );
        aerodynamicCoefficientInterface->setCurrentForceCoefficients( forceCoefficients );
        aerodynamicCoefficientInterface->setReferenceArea( referenceArea );
        aerodynamicCoefficientInterface->setReferenceLength( referenceLength );

        // Declare aerodynamic force model object with coefficient interface and dynamic pressure.
        AerodynamicForce aerodynamicForce( aerodynamicCoefficientInterface );
        aerodynamicForce.setDynamicPressure( dynamicPressure );

        // Compute force from aerodynamic force class.
        aerodynamicForce.computeForce( boost::make_shared< State >( ), 0.0 );
        Eigen::Vector3d force = aerodynamicForce.getForce( );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 2: test the force model implemented as free function with primitive arguments.
    {
        // Compute force.
        Eigen::Vector3d force = computeAerodynamicForce( dynamicPressure,
                                                         referenceArea, forceCoefficients );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 3: test the force model implemented as free function, with coefficient interface
    //         argument.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterface aerodynamicCoefficientInterface;
        aerodynamicCoefficientInterface.setCurrentForceCoefficients( forceCoefficients );
        aerodynamicCoefficientInterface.setReferenceArea( referenceArea );
        aerodynamicCoefficientInterface.setReferenceLength( referenceLength );

        // Compute aerodynamic force using free function with coefficient interface argument.
        Eigen::Vector3d force = computeAerodynamicForce( dynamicPressure,
                                                         aerodynamicCoefficientInterface );

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 4: test the acceleration model implemented as free function with primitive arguments,
    //         based on the force that can be derived from the computed acceleration.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterface aerodynamicCoefficientInterface;
        aerodynamicCoefficientInterface.setCurrentForceCoefficients( forceCoefficients );
        aerodynamicCoefficientInterface.setReferenceArea( referenceArea );
        aerodynamicCoefficientInterface.setReferenceLength( referenceLength );

        // Compute aerodynamic force from aerodynamic acceleration free function with primitive
        // arguments.
        Eigen::Vector3d force = computeAerodynamicAcceleration(
                    dynamicPressure, aerodynamicCoefficientInterface, mass ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }

    // Test 5: test the acceleration model implemented as free function with coefficient interface
    //         argument, based on the force that can be derived from the computed acceleration.
    {
        // Compute aerodynamic force from aerodynamic acceleration free function with
        // coefficient interface argument.
        Eigen::Vector3d force = computeAerodynamicAcceleration(
                    dynamicPressure, referenceArea, forceCoefficients, mass ) * mass;

        // Check if computed force matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedForce, force, tolerance );
    }
}

//! Test implementation of aerodynamic moment and rotational acceleration models.
BOOST_AUTO_TEST_CASE( testAerodynamicMomentAndRotationalAcceleration )
{
    // Using declarations.
    using tudat::AerodynamicCoefficientInterface;
    using tudat::astrodynamics::moment_models::AerodynamicMoment;
    using tudat::astrodynamics::moment_models::computeAerodynamicMoment;
    using tudat::astrodynamics::rotational_acceleration_models::
    computeAerodynamicRotationalAcceleration;
    using tudat::astrodynamics::states::State;

    // Set force coefficients.
    const Eigen::Vector3d forceCoefficients( 2.6, 6.7, 0.3 );

    // Set moment coefficients.
    const Eigen::Vector3d momentCoefficients( 0.0, 1.0, 0.0 );

    // Set dynamical model parameters.
    const double dynamicPressure = 123.6;
    const double referenceArea = 1.7;
    const double referenceLength = 2.6;
    const double mass = 12.46;

    // Set moment arm used to compute moment due to aerodynamic force.
    const Eigen::Vector3d momentArm( 12.1, 0.0, 0.0 );

    // Compute expected force.
    const Eigen::Vector3d expectedForce = forceCoefficients * dynamicPressure * referenceArea;

    // Calculate expected moment.
    const Eigen::Vector3d expectedMoment = dynamicPressure * referenceArea *
            referenceLength * momentCoefficients;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Test 1: test the moment model implemented as class.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface
                = boost::make_shared< AerodynamicCoefficientInterface >( );
        aerodynamicCoefficientInterface->setCurrentForceCoefficients( forceCoefficients );
        aerodynamicCoefficientInterface->setCurrentMomentCoefficients( momentCoefficients );
        aerodynamicCoefficientInterface->setReferenceArea( referenceArea );
        aerodynamicCoefficientInterface->setReferenceLength( referenceLength );

        // Declare aerodynamic moment model object with coefficient interface and dynamic pressure.
        AerodynamicMoment aerodynamicMoment( aerodynamicCoefficientInterface );
        aerodynamicMoment.setDynamicPressure( dynamicPressure );

        // Compute moment from aerodynamic moment class.
        State dummyState;
        aerodynamicMoment.computeMoment( boost::make_shared< State >( ), 0.0 );
        Eigen::Vector3d moment = aerodynamicMoment.getMoment( );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 2: test the moment model implemented as free function with primitive arguments.
    {
        // Compute aerodynamic moment using free function with primitive arguments.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure, referenceArea,
                                                           referenceLength, momentCoefficients );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 3: test the moment moment implemented as free function with coefficient interface
    //         argument.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterface aerodynamicCoefficientInterface;
        aerodynamicCoefficientInterface.setCurrentForceCoefficients( forceCoefficients );
        aerodynamicCoefficientInterface.setCurrentMomentCoefficients( momentCoefficients );
        aerodynamicCoefficientInterface.setReferenceArea( referenceArea );
        aerodynamicCoefficientInterface.setReferenceLength( referenceLength );

        // Compute aerodynamic moment using free function with coefficient interface argument.
        Eigen::Vector3d moment = computeAerodynamicMoment( dynamicPressure,
                                                           aerodynamicCoefficientInterface );

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 4: test the rotational acceleration model implemented as free function with coefficient
    //         interface argument, based on the force that can be derived from the computed
    //         acceleration.
    {
        // Set coefficients and model parameters in aerodynamics coefficient interface object.
        AerodynamicCoefficientInterface aerodynamicCoefficientInterface;
        aerodynamicCoefficientInterface.setCurrentForceCoefficients( forceCoefficients );
        aerodynamicCoefficientInterface.setCurrentMomentCoefficients( momentCoefficients );
        aerodynamicCoefficientInterface.setReferenceArea( referenceArea );
        aerodynamicCoefficientInterface.setReferenceLength( referenceLength );

        // Compute aerodynamic moment from aerodynamic rotational acceleration free function with
        // primitive arguments.
        Eigen::Vector3d moment = computeAerodynamicRotationalAcceleration(
                    dynamicPressure, aerodynamicCoefficientInterface, mass ) * mass;

        // Check if computed moment matches expected.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMoment, moment, tolerance );
    }

    // Test 5: test the rotational acceleration model implemented as free function with primitive
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
