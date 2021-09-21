/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LIGHT_TIME_SOLUTIONS_H
#define TUDAT_LIGHT_TIME_SOLUTIONS_H

#include <memory>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{
namespace observation_models
{

//! Function to retrieve the default tolerance for the light-time equation solution.
/*!
 *  Function to retrieve the default tolerance for the light-time equation solution. This tolerance denotes the
 *  difference between two subsequent light time solutions (in s) that is deemed acceptable for convergence/
 *  \return Default light-time tolerance for given template arguments.
 */
template< typename ObservationScalarType = double >
ObservationScalarType getDefaultLightTimeTolerance( );


//! Typedef for function calculating light-time correction in light-time calculation loop.
typedef std::function< double(
        const Eigen::Vector6d&, const Eigen::Vector6d&,
        const double, const double ) > LightTimeCorrectionFunction;

//! Class for wrapping a custom light-time correction function
class LightTimeCorrectionFunctionWrapper: public LightTimeCorrection
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrectionFunction Custom light-time correction functions, as a function of transmitter and receiver
     * state and time.
     */
    LightTimeCorrectionFunctionWrapper(
            const LightTimeCorrectionFunction lightTimeCorrectionFunction ):
        LightTimeCorrection( function_wrapper_light_time_correction ),
        lightTimeCorrectionFunction_( lightTimeCorrectionFunction ),
    isWarningProvided_( false ){ }

    //! Function to compute the light-time correction
    /*!
     * Function to compute the custom light-time correction
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \return Light-time correction
     */
    double calculateLightTimeCorrection(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime )
    {
        return lightTimeCorrectionFunction_(
                    transmitterState, receiverState, transmissionTime, receptionTime );
    }

    //! Function to compute the partial derivative of the light-time correction w.r.t. observation time
    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. observation time. NOTE: FUNCTION IS NOT
     * YET IMPLEMENTED, EACH OBJECT PRINTS A WARNING ONCE WHEN THIS FUNCTION IS CALLED.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param fixedLinkEnd Reference link end for observation
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
     * \return Light-time correction w.r.t. observation time
     */
    double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType fixedLinkEnd,
            const LinkEndType linkEndAtWhichPartialIsEvaluated )
    {
        if( !isWarningProvided_ )
        {
            std::cerr << "Warning, light-time partial not yet implemented in LightTimeCorrectionFunctionWrapper." << std::endl;
            isWarningProvided_ = true;
        }

        return 0.0;
    }

    //! Function to compute the partial derivative of the light-time correction w.r.t. link end position
    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. link end position. NOTE: FUNCTION IS NOT
     * YET IMPLEMENTED, EACH OBJECT PRINTS A WARNING ONCE WHEN THIS FUNCTION IS CALLED.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
     * \return Light-time correction w.r.t. link end position
     */
    Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated )
    {
        if( !isWarningProvided_ )
        {
            std::cerr << "Warning, light-time partial not yet implemented in LightTimeCorrectionFunctionWrapper." << std::endl;
            isWarningProvided_ = true;
        }

        return Eigen::Matrix< double, 3, 1 >::Zero( );
    }

private:

    //! Custom light-time correction functions, as a function of transmitter and receiver state and time.
    LightTimeCorrectionFunction lightTimeCorrectionFunction_;

    //! Boolean denoting whether a warning has been provided when calling the partial derivative function(s)
    bool isWarningProvided_;
};

//! Class to calculate the light time between two points.
/*!
 *  This class calculates the light time between two points, of which the state functions
 *  have to be provided. Additionally, light-time corrections (such as tropospheric or
 *  relatvistic corrections) can be applied. The motion of the ends of the link during the
 *  light time is taken into account in the calculations.
 */
template< typename ObservationScalarType = double,
          typename TimeType = double >
class LightTimeCalculator
{
public:


    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    //! Class constructor.
    /*!
     *  This constructor is used to initialize the state functions and light-time correction
     *  objects.
     *  \param positionFunctionOfTransmittingBody State function of transmitter.
     *  \param positionFunctionOfReceivingBody State function of receiver.
     *  \param correctionFunctions List of light-time correction objects.
     *  \param iterateCorrections Boolean determining whether to recalculate the light-time
     *  correction during each iteration.
     */
    LightTimeCalculator(
            const std::function< StateType( const TimeType ) > positionFunctionOfTransmittingBody,
            const std::function< StateType( const TimeType ) > positionFunctionOfReceivingBody,
            const std::vector< std::shared_ptr< LightTimeCorrection > > correctionFunctions =
            std::vector< std::shared_ptr< LightTimeCorrection > >( ),
            const bool iterateCorrections = false ):
        stateFunctionOfTransmittingBody_( positionFunctionOfTransmittingBody ),
        stateFunctionOfReceivingBody_( positionFunctionOfReceivingBody ),
        correctionFunctions_( correctionFunctions ),
        iterateCorrections_( iterateCorrections ),
        currentCorrection_( 0.0 ){ }

    //! Class constructor.
    /*!
     *  This constructor is used to initialize the state functions and light-time functions
     *  \param positionFunctionOfTransmittingBody State function of transmitter.
     *  \param positionFunctionOfReceivingBody State function of receiver.
     *  \param correctionFunctions List of light-time correction functions.
     *  \param iterateCorrections Boolean determining whether to recalculate the light-time
     *  correction during each iteration.
     */
    LightTimeCalculator(
            const std::function< StateType( const TimeType ) > positionFunctionOfTransmittingBody,
            const std::function< StateType( const TimeType ) > positionFunctionOfReceivingBody,
            const std::vector< LightTimeCorrectionFunction > correctionFunctions,
            const bool iterateCorrections = false ):
        stateFunctionOfTransmittingBody_( positionFunctionOfTransmittingBody ),
        stateFunctionOfReceivingBody_( positionFunctionOfReceivingBody ),
        iterateCorrections_( iterateCorrections ),
        currentCorrection_( 0.0 )
    {
        for( unsigned int i = 0; i < correctionFunctions.size( ); i++ )
        {
            correctionFunctions_.push_back(
                        std::make_shared< LightTimeCorrectionFunctionWrapper >(
                                                correctionFunctions.at( i ) ) );
        }
    }

    //! Function to calculate the light time.
    /*!
     *  This function calculates the light time between the link ends defined in the constructor.
     *  The input time can be either at transmission or at reception (default) time.
     *  \param time Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param tolerance Maximum allowed light-time difference between two subsequent iterations
     *  for which solution is accepted.
     *  \return The value of the light time between the link ends.
     */
    ObservationScalarType calculateLightTime( const TimeType time,
                               const bool isTimeAtReception = true,
                               const ObservationScalarType tolerance =
            getDefaultLightTimeTolerance< ObservationScalarType >( ) )
    {
        // Declare and initialize variables for receiver and transmitter state (returned by reference).
        StateType receiverState;
        StateType transmitterState;

        // Calculate light time.
        ObservationScalarType lightTime = calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, time, isTimeAtReception, tolerance );
        return lightTime;
    }

    //! Function to calculate the 'measured' vector from transmitter to receiver.
    /*!
     *  Function to calculate the vector from transmitter at transmitter time to receiver at
     *  reception time.
     *  The input time can be either at transmission or reception (default) time.
     *  \param time Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param tolerance Maximum allowed light-time difference between two subsequent iterations
     *  for which solution is accepted.
     *  \return The vector from the transmitter to the reciever.
     */
    PositionType calculateRelativeRangeVector( const TimeType time,
                                               const bool isTimeAtReception = 1 ,
                                               const ObservationScalarType tolerance =
            getDefaultLightTimeTolerance< ObservationScalarType >( ) )
    {
        // Declare and initialize variables for receiver and transmitter state (returned by reference).
        StateType receiverState;
        StateType transmitterState;

        // Calculate link end states and the determine range vector.
        calculateLightTimeWithLinkEndsStates( receiverState, transmitterState,
                                              time, isTimeAtReception, tolerance );
        return ( receiverState - transmitterState ).segment( 0, 3 );
    }

    //! Function to calculate the light time and link-ends states.
    /*!
     *  Function to calculate the transmitter state at transmission time, the receiver state at
     *  reception time, and the light time.
     *  The input time can be either at transmission or reception (default) time.
     *  \param receiverStateOutput Output by reference of receiver state.
     *  \param transmitterStateOutput Output by reference of transmitter state.
     *  \param time Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param tolerance Maximum allowed light-time difference between two subsequent iterations
     *  for which solution is accepted.
     *  \return The value of the light time between the reciever state and the transmitter state.
     */
    ObservationScalarType calculateLightTimeWithLinkEndsStates(
            StateType& receiverStateOutput,
            StateType& transmitterStateOutput,
            const TimeType time,
            const bool isTimeAtReception = 1,
            const ObservationScalarType tolerance =
            ( getDefaultLightTimeTolerance< ObservationScalarType >( ) ) )
    {
        // Initialize reception and transmission times and states to initial guess (zero light time)
        TimeType receptionTime = time;
        TimeType transmissionTime = time;
        StateType receiverState = stateFunctionOfReceivingBody_( receptionTime );
        StateType transmitterState =
                stateFunctionOfTransmittingBody_( transmissionTime );

        // Set initial light-time correction.
        setTotalLightTimeCorrection(
                    transmitterState, receiverState, transmissionTime, receptionTime );

        // Calculate light-time solution assuming infinte speed of signal as initial estimate.
        ObservationScalarType previousLightTimeCalculation =
                calculateNewLightTimeEstime( receiverState, transmitterState );

        // Set variables for iteration
        ObservationScalarType newLightTimeCalculation = 0.0;
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
                receptionTime = time;
                transmissionTime = time - previousLightTimeCalculation;
                transmitterState = ( stateFunctionOfTransmittingBody_( transmissionTime ) );
            }
            else
            {
                receptionTime = time + previousLightTimeCalculation;
                transmissionTime = time;
                receiverState = ( stateFunctionOfReceivingBody_( receptionTime ) );
            }
            newLightTimeCalculation = calculateNewLightTimeEstime( receiverState, transmitterState );

            // Check for convergence.
            if( std::fabs( newLightTimeCalculation - previousLightTimeCalculation ) < tolerance )
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
                if( counter == 50 )
                {
                    isToleranceReached = true;
//                    std::string errorMessage  =
//                            "Warning, light time unconverged at level " +
//                            boost::lexical_cast< std::string >(
//                                std::fabs( newLightTimeCalculation - previousLightTimeCalculation ) ) +
//                            "; current light-time corrections are: "  +
//                            std::to_string( currentCorrection_ ) + " and input time was " +
//                            std::to_string( static_cast< double >( time ) );
//                   std::cerr << errorMessage << std::endl;
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

    //! Function to get the part wrt linkend position
    /*!
     *  Function to get the part wrt linkend position
     *  \param transmitterState State of transmitter.
     *  \param receiverState State of receiver.
     *  \param transmitterTime Time at transmission.
     *  \param receiverTime Time at reiver.
     *  \param isPartialWrtReceiver If partial is to be calculated w.r.t. receiver or transmitter.
     */
    Eigen::Matrix< ObservationScalarType, 1, 3 > getPartialOfLightTimeWrtLinkEndPosition(
            const StateType& transmitterState,
            const StateType& receiverState,
            const TimeType transmitterTime,
            const TimeType receiverTime,
            const bool isPartialWrtReceiver )
    {
        setTotalLightTimeCorrection( transmitterState, receiverState, transmitterTime, receiverTime );

        Eigen::Matrix< ObservationScalarType, 3, 1 > relativePosition =
                receiverState.segment( 0, 3 ) - transmitterState.segment( 0, 3 );
        return ( relativePosition.normalized( ) ).transpose( ) *
                ( mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) +
                  currentCorrection_ / relativePosition.norm( ) ) *
                ( isPartialWrtReceiver ? mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) :
                                         mathematical_constants::getFloatingInteger< ObservationScalarType >( -1 ) );
    }

    //! Function to get list of light-time correction functions
    /*!
     * Function to get list of light-time correction functions
     * \return List of light-time correction functions
     */
    std::vector< std::shared_ptr< LightTimeCorrection > > getLightTimeCorrection( )
    {
        return correctionFunctions_;
    }

protected:

    //! Transmitter state function.
    /*!
     *  Transmitter state function.
     */
    std::function< StateType( const double ) >
    stateFunctionOfTransmittingBody_;

    //! Receiver state function.
    /*!
     *  Receiver state function.
     */
    std::function< StateType( const double ) >
    stateFunctionOfReceivingBody_;

    //! List of light-time correction functions.
    /*!
     *  List of light-time correction functions, i.e. tropospheric, relativistic, etc.
     */
    std::vector< std::shared_ptr< LightTimeCorrection > > correctionFunctions_;

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
    ObservationScalarType calculateNewLightTimeEstime(
            const StateType& receiverState,
            const StateType& transmitterState ) const
    {
        return ( ( ( receiverState - transmitterState ).segment( 0, 3 ) ).
                 template cast< ObservationScalarType >( ) ).norm( ) /
                physical_constants::getSpeedOfLight< ObservationScalarType >( ) + currentCorrection_;
    }

    //! Function to reset the currentCorrection_ variable during current iteration.
    /*!
     *  Function to reset the currentCorrection_ variable during current iteration, representing
     *  the sum of all corrections causing the light time to deviate from the Euclidean value.
     *  \param transmitterState State of transmitter.
     *  \param receiverState State of receiver.
     *  \param transmissionTime Time at transmission.
     *  \param receptionTime Time at reception.
     */
    void setTotalLightTimeCorrection( const StateType& transmitterState,
                                      const StateType& receiverState,
                                      const TimeType transmissionTime ,
                                      const TimeType receptionTime )
    {
        ObservationScalarType totalLightTimeCorrections = mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
        for( unsigned int i = 0; i < correctionFunctions_.size( ); i++ )
        {
            totalLightTimeCorrections += static_cast< ObservationScalarType >(
                        correctionFunctions_[ i ]->calculateLightTimeCorrection(
                            transmitterState.template cast< double >( ), receiverState.template cast< double >( ),
                            static_cast< double >( transmissionTime ), static_cast< double >( receptionTime ) ) );
        }
        currentCorrection_ = totalLightTimeCorrections;
    }
private:
};

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_LIGHT_TIME_SOLUTIONS_H
