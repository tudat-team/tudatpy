/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TEST_LIGHT_TIME_CORRECTIONS_H
#define TUDAT_TEST_LIGHT_TIME_CORRECTIONS_H

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

//! Test correction function for light-time solution unit test.
/*!
 *  Test correction function for light-time solution unit test.
 *  Return value has no physical meaning, but is used to check the correct implementation
 *  of the light-time corrections in unitTestLightTimeSolution.
 *  \param transmitterState State of transmitter at transmission time.
 *  \param receiverState State of receiver at reception time.
 *  \param transmissionTime Time of signal transmission.
 *  \param receptionTime Time of signal reception.
 *  \return Light-time correction value using made-up relation.
 */
double getTimeDifferenceLightTimeCorrection(
        const Eigen::Vector6d transmitterState,
        const Eigen::Vector6d receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    return 1.0E-4 * ( transmissionTime - receptionTime );
}

//! Test correction function for light-time solution unit test.
/*!
 *  Test correction function for light-time solution unit test.
 *  Return value has no physical meaning, but is used to check the correct implementation
 *  of the light-time corrections in unitTestLightTimeSolution.
 *  \param transmitterState State of transmitter at transmission time.
 *  \param receiverState State of receiver at reception time.
 *  \param transmissionTime Time of signal transmission.
 *  \param receptionTime Time of signal reception.
 *  \return Light-time correction value using made-up relation.
 */
double getVelocityDifferenceLightTimeCorrection(
        const Eigen::Vector6d transmitterState,
        const Eigen::Vector6d receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    return 1.0E-4 * ( transmitterState - receiverState ).segment( 3, 3 ).norm( );
}

//! Test correction function for light-time solution unit test
/*!
 *  Test correction function for light-time solution unit test.
 *  Return value has no physical meaning, but is used to check the correct implementation
 *  of the light-time corrections in unitTestLightTimeSolution.
 *  \param transmitterState State of transmitter at transmission time.
 *  \param receiverState State of receiver at reception time.
 *  \param transmissionTime Time of signal transmission.
 *  \param receptionTime Time of signal reception.
 *  \return Light-time correction value using made-up relation.
 */
double getPositionDifferenceLightTimeCorrection(
        const Eigen::Vector6d transmitterState,
        const Eigen::Vector6d receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    return 1.0E-12 * ( transmitterState - receiverState ).x( );
}

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_LIGHT_TIME_CORRECTIONS_H
