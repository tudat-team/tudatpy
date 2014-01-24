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

#ifndef TUDAT_LIGHT_TIME_SOLUTIONS_H
#define TUDAT_LIGHT_TIME_SOLUTIONS_H

#include <map>
#include <vector>

#include <boost/function.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace observation_models
{

//! Typedef for function calculating light-time correction in light-time calculation loop.
typedef boost::function< double(
        const basic_mathematics::Vector6d, const basic_mathematics::Vector6d,
        const double, const double ) > LightTimeCorrectionFunction;

//! Class to calculate the light time between two points.
/*!
 *  This class calculates the light time between two points, of which the state functions
 *  have to be provided. Additionally, light-time corrections (such as tropospheric or
 *  relatvistic corrections) can be applied. The motion of the ends of the link during the
 *  light time is taken into account in the calculations.
 */
class LightTimeCalculator
{
public:

    //! Class constructor.
    /*!
     *  This constructor is used to initialize the state functions and light-time correction
     *  functions.
     *  \param positionFunctionOfTransmittingBody State function of transmitter.
     *  \param positionFunctionOfReceivingBody State function of receiver.
     *  \param correctionFunctions List of light-time correction functions.
     *  \param iterateCorrections Boolean determining whether to recalculate the light-time
     *  correction during each iteration.
     */
    LightTimeCalculator(
            const boost::function< basic_mathematics::Vector6d( const double ) >
            positionFunctionOfTransmittingBody,
            const boost::function< basic_mathematics::Vector6d( const double ) >
            positionFunctionOfReceivingBody,
            const std::vector< LightTimeCorrectionFunction > correctionFunctions =
            std::vector< LightTimeCorrectionFunction >( ),
            const bool iterateCorrections = false ):
        stateFunctionOfTransmittingBody_( positionFunctionOfTransmittingBody ),
        stateFunctionOfReceivingBody_( positionFunctionOfReceivingBody ),
        correctionFunctions_( correctionFunctions ),
        iterateCorrections_( iterateCorrections ),
        currentCorrection_( 0.0 ){ }

    //! Function to calculate the light time.
    /*!
     *  This function calculates the light time between the link ends defined in the constructor.
     *  The input time can be either at transmission or at reception (default) time.
     *  \param ephemerisTime Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param tolerance Maximum allowed light-time difference between two subsequent iterations
     *  for which solution is accepted.
     *  \return The value of the light time between the link ends.
     */
    double calculateLightTime( const double ephemerisTime,
                               const bool isTimeAtReception = true,
                               const double tolerance = 1.0E-12 );

    //! Function to calculate the 'measured' vector from transmitter to receiver.
    /*!
     *  Function to calculate the vector from transmitter at transmitter time to receiver at
     *  reception time.
     *  The input time can be either at transmission or reception (default) time.
     *  \param ephemerisTime Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param tolerance Maximum allowed light-time difference between two subsequent iterations
     *  for which solution is accepted.
     *  \return The vector from the transmitter to the reciever.
     */
    Eigen::Vector3d calculateRelativeRangeVector( const double ephemerisTime,
                                                  const bool isTimeAtReception = true,
                                                  const double tolerance = 1.0E-12 );

    //! Function to calculate the light time and link-ends states.
    /*!
     *  Function to calculate the transmitter state at transmission time, the receiver state at
     *  reception time, and the light time.
     *  The input time can be either at transmission or reception (default) time.
     *  \param receiverStateOutput Output by reference of receiver state.
     *  \param transmitterStateOutput Output by reference of transmitter state.
     *  \param ephemerisTime Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param tolerance Maximum allowed light-time difference between two subsequent iterations
     *  for which solution is accepted.
     *  \return The value of the light time between the reciever state and the transmitter state.
     */
    double calculateLightTimeWithLinkEndsStates(
            basic_mathematics::Vector6d& receiverStateOutput,
            basic_mathematics::Vector6d& transmitterStateOutput,
            const double ephemerisTime,
            const bool isTimeAtReception = true,
            const double tolerance = 1.0E-12 );

protected:

    //! Transmitter state function.
    /*!
     *  Transmitter state function.
     */
    boost::function< basic_mathematics::Vector6d( const double ) >
    stateFunctionOfTransmittingBody_;

    //! Receiver state function.
    /*!
     *  Receiver state function.
     */
    boost::function< basic_mathematics::Vector6d( const double ) >
    stateFunctionOfReceivingBody_;

    //! List of light-time correction functions.
    /*!
     *  List of light-time correction functions, i.e. tropospheric, relativistic, etc.
     */
    std::vector< LightTimeCorrectionFunction > correctionFunctions_;

    //! Boolean deciding whether to recalculate the correction during each iteration.
    /*!
     *  Boolean deciding whether to recalculate the correction during each iteration.
     *  If it is set true, the corrections are calculated during each iteration of the
     *  light-time calculations. If it is set to false, it is calculated once at the begining.
     *  Additionally, when convergence is reached, it is recalculated to check
     *  whether the light time with new correction violates the convergence. If so,
     *  another iteration is performed.
     */
    bool iterateCorrections_;

    //! Current light-time correction.
    /*!
     *  Current light-time correction.
     */
    double currentCorrection_;

    //! Function to calculate a new light-time estimate from the link-ends states.
    /*!
     *  Function to calculate a new light-time estimate from the states of the two ends of the
     *  link. This function recalculates the light time each iteration from the assumed
     *  receiver/transmitter state, as well as the currentCorrection_ variable.
     *  \param receiverState Assumed state of receiver.
     *  \param transmitterState Assumed state of transmitter.
     *  \return New value of the light-time estimate.
     */
    double calculateNewLightTimeEstime( basic_mathematics::Vector6d receiverState,
                                        basic_mathematics::Vector6d transmitterState ) const;

    //! Function to reset the currentCorrection_ variable during current iteration.
    /*!
     *  Function to reset the currentCorrection_ variable during current iteration, representing
     *  the sum of all corrections causing the light time to deviate from the Euclidean value.
     *  \param transmitterState State of transmitter.
     *  \param receiverState State of receiver.
     *  \param transmissionTime Time at transmission.
     *  \param receptionTime Time at reception.
     */
    void setTotalLightTimeCorrection( const basic_mathematics::Vector6d& transmitterState,
                                      const basic_mathematics::Vector6d& receiverState,
                                      const double transmissionTime,
                                      const double receptionTime );
private:
};

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_LIGHT_TIME_SOLUTIONS_H
