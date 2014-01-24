/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      130226    D. Dirkx          File created.
 *      130522    E.D. Brandon      Minor changes during code check.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_TEST_LIGHT_TIME_CORRECTIONS_H
#define TUDAT_TEST_LIGHT_TIME_CORRECTIONS_H

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

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
        const basic_mathematics::Vector6d transmitterState,
        const basic_mathematics::Vector6d receiverState,
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
        const basic_mathematics::Vector6d transmitterState,
        const basic_mathematics::Vector6d receiverState,
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
        const basic_mathematics::Vector6d transmitterState,
        const basic_mathematics::Vector6d receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    return 1.0E-12 * ( transmitterState - receiverState ).x( );
}

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_LIGHT_TIME_CORRECTIONS_H
