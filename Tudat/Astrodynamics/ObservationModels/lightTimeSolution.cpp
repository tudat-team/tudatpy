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
 *      130226    D. Dirkx          Migrated from personal code.
 *      130522    E.D. Brandon      Minor changes during code check.
 *
 *    References
 *
 *    Notes
 *
 */

#include <iostream>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{
namespace observation_models
{

//! Function to calculate the light time.
double LightTimeCalculator::calculateLightTime( const double ephemerisTime,
                                                const bool isTimeAtReception,
                                                const double tolerance )
{
    // Declare and initialize variables for receiver and transmitter state (returned by reference).
    basic_mathematics::Vector6d receiverState = basic_mathematics::Vector6d::Zero( );
    basic_mathematics::Vector6d transmitterState = basic_mathematics::Vector6d::Zero( );

    // Calculate light time.
    const double lightTime = calculateLightTimeWithLinkEndsStates(
                receiverState, transmitterState, ephemerisTime, isTimeAtReception, tolerance );
    return lightTime;
}

//! Function to calculate the 'measured' vector from transmitter to receiver.
Eigen::Vector3d LightTimeCalculator::calculateRelativeRangeVector( const double ephemerisTime,
                                                                   const bool isTimeAtReception,
                                                                   const double tolerance )
{
    // Declare and initialize variables for receiver and transmitter state (returned by reference).
    basic_mathematics::Vector6d receiverState = basic_mathematics::Vector6d::Zero( );
    basic_mathematics::Vector6d transmitterState = basic_mathematics::Vector6d::Zero( );

    // Calculate link end states and the determine range vector.
    calculateLightTimeWithLinkEndsStates( receiverState, transmitterState,
                                          ephemerisTime, isTimeAtReception, tolerance );
    return ( receiverState - transmitterState ).segment( 0, 3 );
}

//! Function to calculate the light time and link-ends states.
double LightTimeCalculator::calculateLightTimeWithLinkEndsStates(
        basic_mathematics::Vector6d& receiverStateOutput,
        basic_mathematics::Vector6d& transmitterStateOutput,
        const double ephemerisTime,
        const bool isTimeAtReception,
        const double tolerance )
{
    using physical_constants::SPEED_OF_LIGHT;
    using std::fabs;

    // Initialize reception and transmission times and states to initial guess (zero light time)
    double receptionTime = ephemerisTime;
    double transmissionTime = ephemerisTime;
    basic_mathematics::Vector6d receiverState = stateFunctionOfReceivingBody_( receptionTime );
    basic_mathematics::Vector6d transmitterState =
            stateFunctionOfTransmittingBody_( transmissionTime );

    // Set initial light-time correction.
    setTotalLightTimeCorrection(
                transmitterState, receiverState, transmissionTime, receptionTime );

    // Calculate light-time solution assuming infinte speed of signal as initial estimate.
    double previousLightTimeCalculation =
            calculateNewLightTimeEstime( receiverState, transmitterState );

    // Set variables for iteration
    double newLightTimeCalculation = 0.0;
    bool isToleranceReached = false;

    // Recalculate light-time solution until tolerance is reached.
    int counter = 0;

    // Set variable determining whether to update the light time each iteration.
    bool updateLightTimeCorrections = false;
    if( iterateCorrections_ )
    {
        updateLightTimeCorrections = true;
    }

    // Iterate until tolerance reached.
    while( !isToleranceReached )
    {
        // Update light-time corrections, if necessary.
        if( updateLightTimeCorrections )
        {
            setTotalLightTimeCorrection(
                        transmitterState, receiverState, transmissionTime, receptionTime );
        }

        // Update light-time estimate for this iteration.
        if( isTimeAtReception )
        {
            receptionTime = ephemerisTime;
            transmissionTime = ephemerisTime - previousLightTimeCalculation;
            transmitterState = ( stateFunctionOfTransmittingBody_( transmissionTime ) );
        }
        else
        {
            receptionTime = ephemerisTime + previousLightTimeCalculation;
            transmissionTime = ephemerisTime;
            receiverState = ( stateFunctionOfReceivingBody_( receptionTime ) );
        }
        newLightTimeCalculation = calculateNewLightTimeEstime( receiverState, transmitterState );

        // Check for convergence.
        if( fabs( newLightTimeCalculation - previousLightTimeCalculation ) < tolerance )
        {
            // If convergence reached, but light-time corrections not iterated,
            // perform 1 more iteration to check for change in correction.
            if( !updateLightTimeCorrections )
            {
                updateLightTimeCorrections = true;
            }
            else
            {
                isToleranceReached = true;
            }
        }
        else
        {
            // Get out of infinite loop (for instance due to low accuracy state functions,
            // to stringent tolerance or limit case for trop. corrections).
            if( counter == 20 )
            {
                isToleranceReached = true;
                std::cerr << "Warning, light time unconverged at level " <<
                             fabs( newLightTimeCalculation - previousLightTimeCalculation )
                          << std::endl << "Current light-time corrections are: " <<
                           currentCorrection_<< " and input time was " << ephemerisTime
                          << std::endl;
            }

            // Update light time for new iteration.
            previousLightTimeCalculation = newLightTimeCalculation;
        }

        counter++;
    }

    // Set output variables and return the light time.
    receiverStateOutput = receiverState;
    transmitterStateOutput = transmitterState;
    return newLightTimeCalculation;
}

//! Function to reset the currentCorrection_ variable during current iteration.
void LightTimeCalculator::setTotalLightTimeCorrection(
        const basic_mathematics::Vector6d& transmitterState,
        const basic_mathematics::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    double totalLightTimeCorrections = 0.0;
    for( unsigned int i = 0; i < correctionFunctions_.size( ); i++ )
    {
        totalLightTimeCorrections += correctionFunctions_[ i ]
                ( transmitterState, receiverState, transmissionTime, receptionTime );
    }
    currentCorrection_ = totalLightTimeCorrections;
}

//! Function to calculate a new light-time estimate from the link-ends states.
double LightTimeCalculator::calculateNewLightTimeEstime(
        basic_mathematics::Vector6d receiverState,
        basic_mathematics::Vector6d transmitterState ) const
{
    return ( receiverState - transmitterState ).segment( 0, 3 ).norm( ) /
            physical_constants::SPEED_OF_LIGHT + currentCorrection_;
}

} // namespace observation_models
} // namespace tudat
