/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"
#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_yarkovsky_acceleration )

const double yarkovskyParameter = -2.88e-14 * 6.839e-12;
const double AU = physical_constants::ASTRONOMICAL_UNIT;


BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationVerySimple )
{
    const Eigen::Vector6d state = { AU, 0.0, 0.0, 0.0, AU, 0.0 };

    const Eigen::Vector3d expectedYarkovskyAcceleration = Eigen::Vector3d{ 0.0, yarkovskyParameter, 0.0 };
    const Eigen::Vector3d computedYarkovskyAcceleration = electromagnetism::computeYarkovskyAcceleration(
            yarkovskyParameter, state );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.e-10 );
}


BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationSimple )
{
    const Eigen::Vector6d state = { 2.0 * AU, AU, 1.0e10, 3.0e5, 4.0e5, 12.0e5 };

    const Eigen::Vector3d expectedYarkovskyAcceleration =
            yarkovskyParameter * 0.19982142476807888 * Eigen::Vector3d{ 3.0, 4.0, 12.0 } / 13.0;
    const Eigen::Vector3d computedYarkovskyAcceleration = electromagnetism::computeYarkovskyAcceleration(
            yarkovskyParameter, state );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.0e-10 );
}

Eigen::Vector6d bodyState = Eigen::Vector6d{ AU, 0, 0, 12.0e5, 3.0e5, 4.0e5 };

Eigen::Vector6d getBodyState( ) { return bodyState; }

Eigen::Vector6d centralBodyState = Eigen::Vector6d::Zero( );

Eigen::Vector6d getCentralBodyState( ) { return centralBodyState; }

BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationModelClassUpdateMembers )
{
    electromagnetism::YarkovskyAccelerationPointer yarkovskyAccelerationModel = std::make_shared< electromagnetism::YarkovskyAcceleration >(
            yarkovskyParameter, &getBodyState, &getCentralBodyState );

    yarkovskyAccelerationModel->updateMembers( 0.0 );
    const Eigen::Vector3d expectedYarkovskyAcceleration = yarkovskyParameter * Eigen::Vector3d{ 12.0, 3.0, 4.0 } / 13.0;
    const Eigen::Vector3d computedYarkovskyAcceleration = yarkovskyAccelerationModel->getAcceleration( );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.0e-10 )
}

BOOST_AUTO_TEST_SUITE_END( )

}

}