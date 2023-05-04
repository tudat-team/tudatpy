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

#include <functional>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/astro/observation_models/observationModel.h"

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
        const double, const double ) > LightTimeCorrectionFunctionSingleLeg;

typedef std::function< double(
        const std::vector< Eigen::Vector6d >&, const std::vector< double >&,
        const unsigned int ) > LightTimeCorrectionFunctionMultiLeg;

enum LightTimeFailureHandling
{
    accept_without_warning,
    print_warning_and_accept,
    throw_exception
};

class LightTimeConvergenceCriteria
{
public:
    LightTimeConvergenceCriteria(
            const bool iterateCorrections = false,
            const int maximumNumberOfIterations = 50,
            const double fractionOfLightTimeTolerance = TUDAT_NAN,
            const LightTimeFailureHandling failureHandling = accept_without_warning):
            iterateCorrections_( iterateCorrections ),
            maximumNumberOfIterations_( maximumNumberOfIterations ),
            failureHandling_( failureHandling ),
            fractionOfLightTimeTolerance_( fractionOfLightTimeTolerance ) { }

    virtual ~LightTimeConvergenceCriteria( ){ }

    template< typename ScalarType = double >
    double getFractionOfLightTimeTolerance( )
    {
        if( fractionOfLightTimeTolerance_ == fractionOfLightTimeTolerance_ )
        {
            return fractionOfLightTimeTolerance_;
        }
        else
        {
            return getDefaultLightTimeTolerance< ScalarType >( );
        }

    }
    bool iterateCorrections_;

    int maximumNumberOfIterations_;

    LightTimeFailureHandling failureHandling_;

protected:
    double fractionOfLightTimeTolerance_;

};

class MultiLegLightTimeConvergenceCriteria: public LightTimeConvergenceCriteria
{
public:
    MultiLegLightTimeConvergenceCriteria(
            const bool iterateCorrections = false,
            const bool iterateMultiLegLightTime = true,
            const int maximumNumberOfIterations = 50,
            const double absoluteTolerance = TUDAT_NAN,
            const LightTimeFailureHandling failureHandling = accept_without_warning):
        LightTimeConvergenceCriteria( iterateCorrections, maximumNumberOfIterations, absoluteTolerance, failureHandling ),
        iterateMultiLegLightTime_( iterateMultiLegLightTime )
        { }

    MultiLegLightTimeConvergenceCriteria(
            const LightTimeConvergenceCriteria singleLegConvergenceCriteria,
            const bool iterateMultiLegLightTime = true ):
        LightTimeConvergenceCriteria( singleLegConvergenceCriteria ),
        iterateMultiLegLightTime_( iterateMultiLegLightTime )
        { }

    bool iterateMultiLegLightTime_;

protected:

};

template< typename ObservationScalarType, typename TimeType >
bool isSingleLegLightTimeSolutionConverged(
        const std::shared_ptr< LightTimeConvergenceCriteria > convergenceCriteria,
        const ObservationScalarType previousLightTimeCalculation,
        const ObservationScalarType newLightTimeCalculation,
        const int numberOfIterations,
        const double currentCorrection,
        const TimeType& currentTime,
        bool& updateLightTimeCorrections )
{
    bool isToleranceReached = false;
    // Check for convergence.
    if( std::fabs( newLightTimeCalculation - previousLightTimeCalculation ) <
            convergenceCriteria->getFractionOfLightTimeTolerance< ObservationScalarType >( ) * newLightTimeCalculation )
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
        if( numberOfIterations == convergenceCriteria->maximumNumberOfIterations_ )
        {
            std::string errorMessage  =
                    "light time unconverged at level " +
                    std::to_string( std::fabs( newLightTimeCalculation - previousLightTimeCalculation ) ) +
                    "; current light-time corrections are: "  +
                    std::to_string( currentCorrection ) + " and current time was " +
                    std::to_string( static_cast< double >( currentTime ) );
            switch( convergenceCriteria->failureHandling_ )
            {
            case accept_without_warning:
                isToleranceReached = true;
                break;
            case print_warning_and_accept:
                std::cerr<<"Warning, "<<errorMessage<<std::endl;
                break;
            case throw_exception:
                throw std::runtime_error( "Error, " + errorMessage );
                break;
            default:
                throw std::runtime_error( "Error, did not recognize light ime failure handling; " + errorMessage );
            }
        }
    }
    return isToleranceReached;
}

template< typename ObservationScalarType, typename TimeType >
bool isMultiLegLightTimeSolutionConverged(
        const std::shared_ptr< LightTimeConvergenceCriteria > convergenceCriteria,
        const ObservationScalarType previousLightTimeCalculation,
        const ObservationScalarType newLightTimeCalculation,
        const int numberOfIterations,
        const TimeType& currentTime )
{
    bool isToleranceReached = false;
    // Check for convergence.
    if( std::fabs( newLightTimeCalculation - previousLightTimeCalculation ) <
            convergenceCriteria->getFractionOfLightTimeTolerance< ObservationScalarType >( ) * newLightTimeCalculation )
    {
        isToleranceReached = true;
    }
    else
    {
        // Get out of infinite loop (for instance due to low accuracy state functions,
        // to stringent tolerance or limit case for trop. corrections).
        if( numberOfIterations == convergenceCriteria->maximumNumberOfIterations_ )
        {
            std::string errorMessage  =
                    "multi-leg light time unconverged at level " +
                    std::to_string( std::fabs( newLightTimeCalculation - previousLightTimeCalculation ) ) +
                    " seconds. Current reference time was " +
                    std::to_string( static_cast< double >( currentTime ) );
            switch( convergenceCriteria->failureHandling_ )
            {
            case accept_without_warning:
                isToleranceReached = true;
                break;
            case print_warning_and_accept:
                std::cerr<< "Warning, " << errorMessage << std::endl;
                break;
            case throw_exception:
                throw std::runtime_error( "Error, " + errorMessage );
                break;
            default:
                throw std::runtime_error( "Error, did not recognize light time failure handling; " + errorMessage );
            }
        }
    }
    return isToleranceReached;
}

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
            const LightTimeCorrectionFunctionSingleLeg lightTimeCorrectionFunction ):
        LightTimeCorrection( function_wrapper_light_time_correction ),
    isWarningProvided_( false )
    {
        lightTimeCorrectionFunction_ = [=] (
                const std::vector< Eigen::Vector6d >& linkEndsStates,
                const std::vector< double >& linkEndsTimes,
                const unsigned int currentMultiLegTransmitterIndex ) -> double
        {
            return lightTimeCorrectionFunction( linkEndsStates.at( currentMultiLegTransmitterIndex ),
                                                linkEndsStates.at( currentMultiLegTransmitterIndex + 1 ),
                                                linkEndsTimes.at( currentMultiLegTransmitterIndex ),
                                                linkEndsTimes.at( currentMultiLegTransmitterIndex + 1) );
        };
    }

    LightTimeCorrectionFunctionWrapper(
            const LightTimeCorrectionFunctionMultiLeg lightTimeCorrectionFunction ):
        LightTimeCorrection( function_wrapper_light_time_correction ),
        lightTimeCorrectionFunction_( lightTimeCorrectionFunction ),
    isWarningProvided_( false )
    { }

    //! Function to compute the light-time correction
    /*!
     * Function to compute the custom light-time correction
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \return Light-time correction
     */
     double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex ) override
    {
         return lightTimeCorrectionFunction_(
                 linkEndsStates, linkEndsTimes, currentMultiLegTransmitterIndex );
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
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
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
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) override
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
    LightTimeCorrectionFunctionMultiLeg lightTimeCorrectionFunction_;

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
     */
    LightTimeCalculator(
            const std::function< StateType( const TimeType ) > positionFunctionOfTransmittingBody,
            const std::function< StateType( const TimeType ) > positionFunctionOfReceivingBody,
            const std::vector< std::shared_ptr< LightTimeCorrection > > correctionFunctions =
                std::vector< std::shared_ptr< LightTimeCorrection > >( ),
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        stateFunctionOfTransmittingBody_( positionFunctionOfTransmittingBody ),
        stateFunctionOfReceivingBody_( positionFunctionOfReceivingBody ),
        correctionFunctions_( correctionFunctions ),
        lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ),
        currentCorrection_( 0.0 )
        { }

    //! Class constructor.
    /*!
     *  This constructor is used to initialize the state functions and light-time functions
     *  \param positionFunctionOfTransmittingBody State function of transmitter.
     *  \param positionFunctionOfReceivingBody State function of receiver.
     *  \param correctionFunctions List of light-time correction functions.
     */
    LightTimeCalculator(
            const std::function< StateType( const TimeType ) > positionFunctionOfTransmittingBody,
            const std::function< StateType( const TimeType ) > positionFunctionOfReceivingBody,
            const std::vector< LightTimeCorrectionFunctionMultiLeg > correctionFunctions,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        stateFunctionOfTransmittingBody_( positionFunctionOfTransmittingBody ),
        stateFunctionOfReceivingBody_( positionFunctionOfReceivingBody ),
        lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ),
        currentCorrection_( 0.0 )
    {
        for( unsigned int i = 0; i < correctionFunctions.size( ); i++ )
        {
            correctionFunctions_.push_back(
                        std::make_shared< LightTimeCorrectionFunctionWrapper >(
                            correctionFunctions.at( i ) ) );
        }
    }

    LightTimeCalculator(
            const std::function< StateType( const TimeType ) > positionFunctionOfTransmittingBody,
            const std::function< StateType( const TimeType ) > positionFunctionOfReceivingBody,
            const std::vector< LightTimeCorrectionFunctionSingleLeg > correctionFunctions,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        stateFunctionOfTransmittingBody_( positionFunctionOfTransmittingBody ),
        stateFunctionOfReceivingBody_( positionFunctionOfReceivingBody ),
        lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ),
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
                               const bool isTimeAtReception = true )
    {
        // Declare and initialize variables for receiver and transmitter state (returned by reference).
        StateType receiverState;
        StateType transmitterState;

        // Calculate light time.
        ObservationScalarType lightTime = calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, time, isTimeAtReception );
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
                                               const bool isTimeAtReception = 1 )
    {
        // Declare and initialize variables for receiver and transmitter state (returned by reference).
        StateType receiverState;
        StateType transmitterState;

        // Calculate link end states and the determine range vector.
        calculateLightTimeWithLinkEndsStates( receiverState, transmitterState,
                                              time, isTimeAtReception );
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
            const bool isTimeAtReception = 1 )
    {
        std::vector< StateType > linkEndsStates( 2, StateType::Constant( TUDAT_NAN ) );
        std::vector< TimeType > linkEndsTimes( 2, TUDAT_NAN );

        ObservationScalarType lightTime = calculateLightTimeWithMultiLegLinkEndsStates(
                linkEndsStates, linkEndsTimes, time, isTimeAtReception, 0 );

        transmitterStateOutput = linkEndsStates.at( 0 );
        receiverStateOutput = linkEndsStates.at( 1 );

        return lightTime;
    }

    //! Function to calculate the light time and link-ends states, given an initial guess for all legs.
    /*!
     *  Function to calculate the transmitter state at transmission time, the receiver state at
     *  reception time, and the light time. Calculation uses the states and times at each link end of the model (which
     *  are provided as argument) as initial guess.
     *  The input time can be either at transmission or reception (default) time.
     *  \param linkEndsStates Input/output link end states over all legs of model.
     *  \param linkEndsTimes Input/output link end times over all legs of model.
     *  \param time Time at reception or transmission.
     *  \param isTimeAtReception True if input time is at reception, false if at transmission.
     *  \param currentMultiLegTransmitterIndex Index of current transmitter in multi-leg model
     *  \return The value of the light time between the reciever state and the transmitter state.
     */
    ObservationScalarType calculateLightTimeWithMultiLegLinkEndsStates(
        std::vector< StateType >& linkEndsStates,
        std::vector< TimeType >& linkEndsTimes,
        const TimeType time,
        const bool isTimeAtReception = true,
        const unsigned int currentMultiLegTransmitterIndex = 0 )
    {
        const unsigned int currentMultiLegReceiverIndex = currentMultiLegTransmitterIndex + 1;

        if ( linkEndsStates.size( ) != linkEndsTimes.size( ) || currentMultiLegReceiverIndex > linkEndsTimes.size( ) - 1 )
        {
            throw std::runtime_error( "Error when calculating light time with multi-leg information: size of provided"
                                      "state and time vectors is inconsistent." );
        }

        // Initialize reception and transmission times and states to initial guess: zero light time
        TimeType receptionTime = time, transmissionTime = time;
        StateType receiverState, transmitterState;
        // If initial guess is provided as input, use that instead
        if ( !isTimeAtReception && !std::isnan( linkEndsTimes.at( currentMultiLegReceiverIndex ) ) )
        {
            receptionTime = linkEndsTimes.at( currentMultiLegReceiverIndex );
        }
        else if ( isTimeAtReception && !std::isnan( linkEndsTimes.at( currentMultiLegTransmitterIndex ) ) )
        {
            transmissionTime = linkEndsTimes.at( currentMultiLegTransmitterIndex );
        }
        receiverState = stateFunctionOfReceivingBody_( receptionTime );
        transmitterState = stateFunctionOfTransmittingBody_( transmissionTime );

        // Set initial light-time correction.
        updateCurrentLinkEndStatesAndTimes(
                linkEndsTimes, linkEndsStates, currentMultiLegTransmitterIndex, receptionTime, transmissionTime,
                receiverState, transmitterState );
        setTotalLightTimeCorrection( linkEndsStates, linkEndsTimes, currentMultiLegTransmitterIndex );

        // Calculate light-time solution
        ObservationScalarType previousLightTimeCalculation = calculateNewLightTimeEstimate(
                receiverState, transmitterState );

        // Set variables for iteration
        ObservationScalarType newLightTimeCalculation = 0.0;
        bool isToleranceReached = false;

        // Recalculate light-time solution until tolerance is reached.
        int counter = 0;

        // Set variable determining whether to update the light time each iteration.
        bool updateLightTimeCorrections = false;
        if( lightTimeConvergenceCriteria_->iterateCorrections_ )
        {
            updateLightTimeCorrections = true;
        }

        // Iterate until tolerance reached.
        while( !isToleranceReached )
        {
            // Update light-time corrections, if necessary.
            if( updateLightTimeCorrections )
            {
                updateCurrentLinkEndStatesAndTimes(
                    linkEndsTimes, linkEndsStates, currentMultiLegTransmitterIndex, receptionTime, transmissionTime,
                    receiverState, transmitterState );
                setTotalLightTimeCorrection( linkEndsStates, linkEndsTimes, currentMultiLegTransmitterIndex );
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
            newLightTimeCalculation = calculateNewLightTimeEstimate( receiverState, transmitterState );
            isToleranceReached = isSingleLegLightTimeSolutionConverged(
                    lightTimeConvergenceCriteria_, previousLightTimeCalculation, newLightTimeCalculation, counter,
                    currentCorrection_, time, updateLightTimeCorrections );

            // Update light time for new iteration.
            previousLightTimeCalculation = newLightTimeCalculation;
            counter++;
        }

        // Set output variables and return the light time.
        updateCurrentLinkEndStatesAndTimes(
                linkEndsTimes, linkEndsStates, currentMultiLegTransmitterIndex, receptionTime, transmissionTime,
                receiverState, transmitterState );

        std::cerr << std::setprecision(18) << "Light time: " << newLightTimeCalculation <<
            " (" << counter + 1 << " iter) " << std::endl;

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

    //! Function to get the current ideal light time (distance divided by the speed of light)
    ObservationScalarType getCurrentIdealLightTime( )
    {
        return currentIdealLightTime_;
    }

    //! Function to get the value of the current light time corrections
    ObservationScalarType getCurrentLightTimeCorrection( )
    {
        return currentCorrection_;
    }

protected:

    //! Transmitter state function.
    /*!
     *  Transmitter state function.
     */
    std::function< StateType( const TimeType ) > stateFunctionOfTransmittingBody_;

    //! Receiver state function.
    /*!
     *  Receiver state function.
     */
    std::function< StateType( const TimeType ) > stateFunctionOfReceivingBody_;

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

    std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria_;

    //! Current ideal light-time (i.e. without corrections)
    ObservationScalarType currentIdealLightTime_;

    //! Current light-time correction.
    ObservationScalarType currentCorrection_;

    //! Function to calculate a new light-time estimate from the link-ends states.
    /*!
     *  Function to calculate a new light-time estimate from the states of the two ends of the
     *  link. This function recalculates the light time each iteration from the assumed
     *  receiver/transmitter state, as well as the currentCorrection_ variable.
     *  \param receiverState Assumed state of receiver.
     *  \param transmitterState Assumed state of transmitter.
     *  \return New value of the light-time estimate.
     */
    ObservationScalarType calculateNewLightTimeEstimate(
            const StateType& receiverState,
            const StateType& transmitterState )
    {
        currentIdealLightTime_ = ( receiverState - transmitterState ).segment( 0, 3 ).norm( ) /
                physical_constants::getSpeedOfLight< ObservationScalarType >( );

        return currentIdealLightTime_ + currentCorrection_;
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
        std::vector< StateType > linkEndStates = { transmitterState, receiverState };
        std::vector< TimeType > linkEndTimes = { transmissionTime, receptionTime };
        unsigned int currentMultiLegTransmitterIndex = 0;

        setTotalLightTimeCorrection( linkEndStates, linkEndTimes, currentMultiLegTransmitterIndex );
    }

    void setTotalLightTimeCorrection(
            std::vector< StateType >& linkEndStates,
            std::vector< TimeType >& linkEndTimes,
            const unsigned int currentMultiLegTransmitterIndex )
    {
        ObservationScalarType totalLightTimeCorrections = mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );

        std::vector< double > linkEndTimesDouble( linkEndTimes.size( ) );
        std::vector< Eigen::Vector6d > linkEndStatesDouble( linkEndStates.size( ) );
        for ( unsigned int i = 0; i < linkEndTimesDouble.size( ); ++i )
        {
            linkEndTimesDouble.at( i ) = static_cast< double >( linkEndTimes.at( i ) );
            linkEndStatesDouble.at( i ) = linkEndStates.at( i ).template cast< double >( );
        }

        for( unsigned int i = 0; i < correctionFunctions_.size( ); i++ )
        {
            totalLightTimeCorrections += static_cast< ObservationScalarType >(
                        correctionFunctions_[ i ]->calculateLightTimeCorrectionWithMultiLegLinkEndStates(
                            linkEndStatesDouble, linkEndTimesDouble, currentMultiLegTransmitterIndex ) );
        }
        currentCorrection_ = totalLightTimeCorrections;
    }

private:

    void updateCurrentLinkEndStatesAndTimes( std::vector< TimeType >& linkEndsTimes,
                                             std::vector< StateType >& linkEndsStates,
                                             const unsigned int currentMultiLegTransmitterIndex,
                                             const TimeType receptionTime,
                                             const TimeType transmissionTime,
                                             const StateType receiverState,
                                             const StateType transmitterState )
    {
        linkEndsTimes.at( currentMultiLegTransmitterIndex + 1 ) = receptionTime;
        linkEndsTimes.at( currentMultiLegTransmitterIndex ) = transmissionTime;
        linkEndsStates.at( currentMultiLegTransmitterIndex + 1 ) = receiverState;
        linkEndsStates.at( currentMultiLegTransmitterIndex ) = transmitterState;
    }
};

template< typename ObservationScalarType = double, typename TimeType = double >
class MultiLegLightTimeCalculator
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    MultiLegLightTimeCalculator(
            const std::vector< std::shared_ptr< LightTimeCalculator < ObservationScalarType, TimeType > > >
                    lightTimeCalculators,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< MultiLegLightTimeConvergenceCriteria >( ) ):
        lightTimeCalculators_( lightTimeCalculators ),
        numberOfLinks_( lightTimeCalculators.size( ) ),
        numberOfLinkEnds_( lightTimeCalculators.size( ) + 1 )
    {
        // Check if convergence criteria is of type MultiLegLightTimeConvergenceCriteria. If not, create it. This will
        // effectively set the value of whether to iterate the multi-leg corrections to the
        // MultiLegLightTimeConvergenceCriteria default
        lightTimeConvergenceCriteria_ = std::dynamic_pointer_cast< MultiLegLightTimeConvergenceCriteria >( lightTimeConvergenceCriteria );
        if ( lightTimeConvergenceCriteria_ == nullptr )
        {
            lightTimeConvergenceCriteria_ = std::make_shared< MultiLegLightTimeConvergenceCriteria >(
                    *lightTimeConvergenceCriteria );
        }
    }

    // Constructor for a single leg
    MultiLegLightTimeCalculator(
            const std::function< StateType( const TimeType ) > positionFunctionOfTransmittingBody,
            const std::function< StateType( const TimeType ) > positionFunctionOfReceivingBody,
            const std::vector< std::shared_ptr< LightTimeCorrection > > correctionFunctions =
                std::vector< std::shared_ptr< LightTimeCorrection > >( ),
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< MultiLegLightTimeConvergenceCriteria >( ) ):
        numberOfLinks_( 1 ),
        numberOfLinkEnds_( 2 )
    {
        lightTimeCalculators_.clear( );
        lightTimeCalculators_.push_back( std::make_shared< LightTimeCalculator< ObservationScalarType, TimeType > >(
                positionFunctionOfTransmittingBody, positionFunctionOfReceivingBody, correctionFunctions,
                lightTimeConvergenceCriteria ) );

        // Check if convergence criteria is of type MultiLegLightTimeConvergenceCriteria. If not, create it. This will
        // effectively set the value of whether to iterate the multi-leg corrections to the
        // MultiLegLightTimeConvergenceCriteria default
        lightTimeConvergenceCriteria_ = std::dynamic_pointer_cast< MultiLegLightTimeConvergenceCriteria >( lightTimeConvergenceCriteria );
        if ( lightTimeConvergenceCriteria_ == nullptr )
        {
            lightTimeConvergenceCriteria_ = std::make_shared< MultiLegLightTimeConvergenceCriteria >(
                    *lightTimeConvergenceCriteria );
        }
    }

    ObservationScalarType calculateLightTimeWithLinkEndsStates(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndsTimesOutput,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndsStatesOutput,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings = nullptr )
    {

        if( ancillarySettings != nullptr )
        {
            transmissionReceptionDelays_ = ancillarySettings->getAncilliaryDoubleVectorData( transmission_reception_delays );

            // Delays vector already including delays at receiving and transmitting stations
            if ( transmissionReceptionDelays_.size( ) == numberOfLinkEnds_ )
            { }
            // Delays vector not including delays at receiving and transmitting stations: set them to 0
            else if ( transmissionReceptionDelays_.size( ) == numberOfLinkEnds_ - 2 )
            {
                transmissionReceptionDelays_.insert( transmissionReceptionDelays_.begin( ), 0.0 );
                transmissionReceptionDelays_.push_back( 0.0 );
            }
            else
            {
                throw std::runtime_error(
                        "Error when computing multi-leg light time: size of retransmission delays (" +
                        std::to_string( transmissionReceptionDelays_.size() ) + ") is invalid, should be " +
                        std::to_string( numberOfLinkEnds_ ) + " or " + std::to_string( numberOfLinkEnds_ - 2) + "." );
            }
        }
        else
        {
            for( unsigned int i = 0; i < numberOfLinkEnds_; i++ )
            {
                transmissionReceptionDelays_.push_back( 0.0 );
            }
        }

        // Retrieve index of link end where to start.
        unsigned int startLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndAssociatedWithTime, numberOfLinkEnds_ );

        // If start is not at transmitter or receiver, compute and add retransmission delay.
        if( ( startLinkEndIndex != 0 ) && ( startLinkEndIndex != numberOfLinkEnds_ - 1 ) )
        {
            if ( transmissionReceptionDelays_.at( startLinkEndIndex ) != 0.0 )
            {
                throw std::runtime_error(
                        "Error when computing light time with reference link end that is not receiver or transmitter: "
                        "dealing with non-zero retransmission delays at the reference link end is not implemented. It "
                        "would require distinguishing between reception and transmission delays." );
            }
        }

        // Initialize vectors with states and times
        std::vector< TimeType > linkEndTimes( 2 * numberOfLinks_, TUDAT_NAN );
        std::vector< StateType > linkEndStates( 2 * numberOfLinks_, StateType::Constant( TUDAT_NAN ) );

        ObservationScalarType previousLightTimeEstimate = calculateMultiLegLightTimeEstimate(
                time, startLinkEndIndex, linkEndTimes, linkEndStates );
        ObservationScalarType newLightTimeEstimate;

        // Iterate light times only if necessary. Don't iterate if the model is constituted by a single leg
        iterationCounter_ = 1;
        if ( numberOfLinks_ == 1 || !lightTimeConvergenceCriteria_->iterateMultiLegLightTime_ )
        {
            newLightTimeEstimate = previousLightTimeEstimate;
        }
        else
        {
            bool isToleranceReached = false;

            while( !isToleranceReached )
            {
                newLightTimeEstimate = calculateMultiLegLightTimeEstimate(
                        time, startLinkEndIndex, linkEndTimes, linkEndStates );

                isToleranceReached = isMultiLegLightTimeSolutionConverged(
                    lightTimeConvergenceCriteria_, previousLightTimeEstimate, newLightTimeEstimate, iterationCounter_,
                    time );

                // Update light time for next iteration
                previousLightTimeEstimate = newLightTimeEstimate;

                ++iterationCounter_;
            }
        }

        std::cerr << "Iterations: " << iterationCounter_ << std::endl;


        // Save output
        linkEndsTimesOutput.clear( );
        linkEndsStatesOutput.clear( );
        linkEndsTimesOutput.resize( 2 * numberOfLinks_ );
        linkEndsStatesOutput.resize( 2 * numberOfLinks_ );
        for ( unsigned int i = 0; i < linkEndTimes.size( ); ++i )
        {
            linkEndsStatesOutput.at( i ) = linkEndStates.at( i ).template cast< double >( );
            linkEndsTimesOutput.at( i ) = linkEndTimes.at( i );
        }

        return newLightTimeEstimate;
    }

    ObservationScalarType getTotalIdealLightTime(  )
    {
        ObservationScalarType lightTime = 0;
        for ( unsigned int i = 0; i < lightTimeCalculators_.size( ); ++i )
        {
            lightTime += lightTimeCalculators_.at( i )->getCurrentIdealLightTime( );
        }
        return lightTime;
    }

    ObservationScalarType getTotalLightTimeCorrections(  )
    {
        ObservationScalarType lightTimeCorrections = 0;
        for ( unsigned int i = 0; i < lightTimeCalculators_.size( ); ++i )
        {
            lightTimeCorrections += lightTimeCalculators_.at( i )->getCurrentLightTimeCorrection( );
        }
        return lightTimeCorrections;
    }

    std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > getLightTimeCalculators( )
    {
        return lightTimeCalculators_;
    }

    unsigned int getNumberOfMultiLegIterations( )
    {
        return iterationCounter_;
    }

private:

    ObservationScalarType calculateMultiLegLightTimeEstimate(
            const TimeType time,
            const unsigned int startLinkEndIndex,
            std::vector< TimeType >& linkEndTimes,
            std::vector< StateType >& linkEndStates )
    {
        // Define objects to keep light times
        ObservationScalarType totalLightTime = 0.0;
        ObservationScalarType currentLightTime;

        // Initialize light time with initial delay
        totalLightTime += transmissionReceptionDelays_.at( startLinkEndIndex );

        // Define 'current reception time': time at the receiving antenna
        TimeType currentLinkEndReceptionTime = time - transmissionReceptionDelays_.at( startLinkEndIndex );

        // Move 'backwards' from reference link end to transmitter.
        for( unsigned int currentDownIndex = startLinkEndIndex; currentDownIndex > 0; --currentDownIndex )
        {
            unsigned int transmitterIndex = 2 * ( currentDownIndex - 1 );
            currentLightTime = lightTimeCalculators_.at( currentDownIndex - 1 )->calculateLightTimeWithMultiLegLinkEndsStates(
                        linkEndStates, linkEndTimes, currentLinkEndReceptionTime, true, transmitterIndex );

            // If an additional leg is required, retrieve retransmission delay and update current time
            currentLightTime += transmissionReceptionDelays_.at( currentDownIndex - 1 );
            currentLinkEndReceptionTime -= currentLightTime;

            // Add computed light-time to total time and move to next leg
            totalLightTime += currentLightTime;
        }

        // Define 'current transmission time': time at the transmitting antenna
        TimeType currentLinkEndTransmissionTime = time + transmissionReceptionDelays_.at( startLinkEndIndex );

        // Move 'forwards' from reference link end to receiver.
        for( unsigned int currentUpIndex = startLinkEndIndex; currentUpIndex < numberOfLinkEnds_ - 1; ++currentUpIndex )
        {
            unsigned int transmitterIndex = 2 * currentUpIndex;
            currentLightTime = lightTimeCalculators_.at( currentUpIndex )->calculateLightTimeWithMultiLegLinkEndsStates(
                        linkEndStates, linkEndTimes, currentLinkEndTransmissionTime, false, transmitterIndex );

            // If an additional leg is required, retrieve retransmission delay and update current time
            currentLightTime += transmissionReceptionDelays_.at( currentUpIndex + 1 );
            currentLinkEndTransmissionTime += currentLightTime;

            // Add computed light-time to total time and move to next leg
            totalLightTime += currentLightTime;
        }

        return totalLightTime;
    }

    std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > lightTimeCalculators_;

    std::shared_ptr< MultiLegLightTimeConvergenceCriteria > lightTimeConvergenceCriteria_;

    //! Number of links in multi-leg light-time model
    const unsigned int numberOfLinks_;

    //! Number of link ends in multi-leg light-time model (=number of links + 1)
    const unsigned int numberOfLinkEnds_;

    std::vector< double > transmissionReceptionDelays_;

    unsigned int iterationCounter_;
};

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_LIGHT_TIME_SOLUTIONS_H
