/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <type_traits>

#include <boost/test/unit_test.hpp>

#include <Eigen/Geometry>

#include <Tudat/Basics/tudatTypeTraits.h>
#include <Tudat/Basics/timeType.h>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_tudat_type_traits )

//! Test if Time objects cast to the expected precision
BOOST_AUTO_TEST_CASE( testTypeTraits )
{
    BOOST_CHECK_EQUAL( is_eigen_matrix< Eigen::MatrixXd >::value, true );
    BOOST_CHECK_EQUAL( is_eigen_matrix< Eigen::VectorXd >::value, true );
    BOOST_CHECK_EQUAL( is_eigen_matrix< Eigen::Quaterniond >::value, false );
    BOOST_CHECK_EQUAL( is_eigen_matrix< int >::value, false );
    BOOST_CHECK_EQUAL( is_eigen_matrix< double >::value, false );

    bool isTimeType;
    isTimeType = is_time_type< int >::value;
    BOOST_CHECK_EQUAL( isTimeType, false );

    isTimeType = is_time_type< float >::value;
    BOOST_CHECK_EQUAL( isTimeType, false );

    isTimeType = is_time_type< double >::value;
    BOOST_CHECK_EQUAL( isTimeType, true );

    isTimeType = is_time_type< long double >::value;
    BOOST_CHECK_EQUAL( isTimeType, false );

    isTimeType = is_time_type< Time >::value;
    BOOST_CHECK_EQUAL( isTimeType, true );



    bool isStateType;
    isStateType = is_state_scalar< int >::value;
    BOOST_CHECK_EQUAL( isStateType, false );

    isStateType = is_state_scalar< float >::value;
    BOOST_CHECK_EQUAL( isStateType, false );

    isStateType = is_state_scalar< double >::value;
    BOOST_CHECK_EQUAL( isStateType, true );

    isStateType = is_state_scalar< long double >::value;
    BOOST_CHECK_EQUAL( isStateType, true );

    isStateType = is_state_scalar< Time >::value;
    BOOST_CHECK_EQUAL( isStateType, false );



    bool isStateAndTimeType;
    isStateAndTimeType = is_state_scalar_and_time_type< int, int >::value;
    BOOST_CHECK_EQUAL( isStateAndTimeType, false );

    isStateAndTimeType = is_state_scalar_and_time_type< int, double >::value;
    BOOST_CHECK_EQUAL( isStateAndTimeType, false );

    isStateAndTimeType = is_state_scalar_and_time_type< double, float >::value;
    BOOST_CHECK_EQUAL( isStateAndTimeType, false );

    isStateAndTimeType = is_state_scalar_and_time_type< double, double >::value;
    BOOST_CHECK_EQUAL( isStateAndTimeType, true );

    isStateAndTimeType = is_state_scalar_and_time_type< long double, long double >::value;
    BOOST_CHECK_EQUAL( isStateAndTimeType, false );

    isStateAndTimeType = is_state_scalar_and_time_type< long double, Time >::value;
    BOOST_CHECK_EQUAL( isStateAndTimeType, true );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

