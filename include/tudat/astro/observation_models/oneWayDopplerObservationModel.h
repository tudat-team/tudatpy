/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYDOPPLEROBSERVATIONMODEL_H
#define TUDAT_ONEWAYDOPPLEROBSERVATIONMODEL_H


#include <map>

#include <functional>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include "tudat/astro/relativity/metric.h"
namespace tudat
{

namespace observation_models
{

//! Function to compute the component of a (velocity) vector projected along a unit vector, divided by speed of light.
/*!
 *  Function to compute the component of a (velocity) vector projected along a unit vector, divided by speed of light.
 *  \param lineOfSightUnitVector Unit vector along which velocityVector is to be projected
 *  \param velocityVector Vector for which component along lineOfSightUnitVector is to be computed
 *  \return The component of a (velocity) vector projected along a unit vector, divided by speed of light.
 */
template< typename ObservationScalarType = double >
ObservationScalarType calculateLineOfSightVelocityAsCFraction(
        const Eigen::Matrix< ObservationScalarType, 3, 1 >& lineOfSightUnitVector,
        const Eigen::Matrix< ObservationScalarType, 3, 1 >& velocityVector )
{
    return lineOfSightUnitVector.dot( velocityVector ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );
}

//! Function to compute component of transmitter velocity projected along line-of-sight vector, divided by speed of light.
/*!
 *  Function to compute component of transmitter velocity projected along line-of-sight vector, divided by speed of light,
 *  from the receiver position and transmitter state function. The unit vector is computed in the direction from the
 *  transmitter to the receiver
 *  \param receiverPosition Cartesian position of receiver
 *  \param transmitterStateFunction Function returning the Cartesian state of the transmitter.
 *  \param currentTime Time at which transmitterStateFunction is to be evaluated.
 *  \return Component of transmitter velocity projected along line-of-sight vector to receiver, divided by speed of light.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
ObservationScalarType calculateLineOfSightVelocityAsCFractionFromTransmitterStateFunction(
        const Eigen::Matrix< ObservationScalarType, 3, 1 >& receiverPosition,
        const std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const double ) >& transmitterStateFunction,
        const TimeType currentTime )
{
    Eigen::Matrix< ObservationScalarType, 6, 1 > currentState = transmitterStateFunction( currentTime );
    return calculateLineOfSightVelocityAsCFraction< ObservationScalarType >(
                ( receiverPosition - currentState.segment( 0, 3 ) ).normalized( ), currentState.segment( 3, 3 ) );
}

//! Function to compute component of transmitter velocity projected along line-of-sight vector, divided by speed of light.
/*!
 *  Function to compute component of transmitter velocity projected along line-of-sight vector, divided by speed of light,
 *  from the receiver state function and transmitter position. The unit vector is computed in the direction from the
 *  transmitter to the receiver
 *  \param receiverStateFunction Function returning the Cartesian state of the receiver.
 *  \param transmitterPosition Cartesian position of transmitter
 *  \param currentTime Time at which transmitterStateFunction is to be evaluated.
 *  \return Component of transmitter velocity projected along line-of-sight vector to receiver, divided by speed of light.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
ObservationScalarType calculateLineOfSightVelocityAsCFractionFromReceiverStateFunction(
        const std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const double ) >& receiverStateFunction,
        const Eigen::Matrix< ObservationScalarType, 3, 1 >& transmitterPosition,
        const TimeType currentTime )
{
    Eigen::Matrix< ObservationScalarType, 6, 1 > currentState = receiverStateFunction( currentTime );
    return calculateLineOfSightVelocityAsCFraction< ObservationScalarType >(
                ( currentState.segment( 0, 3 ) - transmitterPosition ).normalized( ), currentState.segment( 3, 3 ) );
}

//! Function to compute first-order (radial) Doppler term from a Taylor series expansion
/*!
 *  Function to compute first-order (radial) Doppler term from a Taylor series expansion. The function computes the
 *  (dt1/dt2 -1) term, with t2 the coordinate reception time and t1 the coordinate transmission time of the signal. Light
 *  time corrections are not included in this function. The Taylor series of the denominator of dt1/dt2 is used in the
 *  calculation, to an order that is provided as input.
 *  \param transmitterState Cartesian state of the transmitter at t1
 *  \param receiverState Cartesian state of the receiver at t2
 *  \param lightTimeWrtTransmitterPositionPartial Light time partial w.r.t. transmitter position
 *  \param lightTimeWrtReceiverPositionPartial Light time partial w.r.t. receiver position
 *  \param taylorSeriesOrder Order to which Taylor series is to be expanded
 *  \return First-order Doppler effect for electromagnetic signal transmission from transmitter to receiver as: (dt1/dt2 -1)
 *   with t2 the coordinate reception time and t1 the coordinate transmission time of the signal.
 */
template< typename ObservationScalarType = double >
ObservationScalarType computeOneWayFirstOrderDopplerTaylorSeriesExpansion(
        Eigen::Matrix< ObservationScalarType, 6, 1 >& transmitterState,
        Eigen::Matrix< ObservationScalarType, 6, 1 >& receiverState,
        Eigen::Matrix< ObservationScalarType, 1, 3 >& lightTimeWrtTransmitterPositionPartial,
        Eigen::Matrix< ObservationScalarType, 1, 3 >& lightTimeWrtReceiverPositionPartial,
        const int taylorSeriesOrder )
{
    // Compute projected velocity components
    ObservationScalarType transmitterTerm  =
            ( -lightTimeWrtTransmitterPositionPartial * ( transmitterState.segment( 3, 3 ) ) /
              physical_constants::getSpeedOfLight< ObservationScalarType >( )  )( 0 );
    ObservationScalarType receiverTerm =
            ( lightTimeWrtReceiverPositionPartial * ( receiverState.segment( 3, 3 ) ) /
              physical_constants::getSpeedOfLight< ObservationScalarType >( ) )( 0 );

    // Compute Taylor series of 1/(1-r21*v2) up to required order
    ObservationScalarType currentTaylorSeriesTerm =
            mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 );
    ObservationScalarType currentTaylorSeries = mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
    for( int i = 0; i < taylorSeriesOrder; i++ )
    {
        currentTaylorSeriesTerm *= transmitterTerm;
        currentTaylorSeries += currentTaylorSeriesTerm;
    }

    // Compute Doppler term
    return -receiverTerm + currentTaylorSeries *
            ( mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) - receiverTerm );
}

//! Function to compute proper time contribution one-way Doppler term from a Taylor series expansion
/*!
 *  Function to compute proper time contribution one-way Doppler term from a Taylor series expansion.
 *  \param transmitterProperTimeRateDifference Derivative of deviation between proper and coordinate time (Delta - t) w.r.t.
 *  coordinate time t at transmitter
 *  \param receiverProperTimeRateDifference Derivative of deviation between proper and coordinate time (Delta - t) w.r.t.
 *  coordinate time t at receiver
 *  \param taylorSeriesOrder Order to which Taylor series is to be expanded
 */
template< typename ObservationScalarType = double >
ObservationScalarType computeDopplerProperTimeInfluenceTaylorSeriesExpansion(
        const ObservationScalarType transmitterProperTimeRateDifference,
        const ObservationScalarType receiverProperTimeRateDifference,
        const int taylorSeriesOrder )
{
    ObservationScalarType currentTaylorSeriesTerm = mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 );
    ObservationScalarType currentTaylorSeries = mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
    for( int i = 0; i < taylorSeriesOrder; i++ )
    {
        currentTaylorSeriesTerm *= -receiverProperTimeRateDifference;
        currentTaylorSeries += currentTaylorSeriesTerm;
    }

    return ( currentTaylorSeries + transmitterProperTimeRateDifference  ) +
            transmitterProperTimeRateDifference * currentTaylorSeries;
}

//! Base class for interface class that is used to compute proper-time rate at a transmitter/receiver for one-way Doppler model
/*!
 *  Base class for interface class that is used to compute proper-time rate at a transmitter/receiver for one-way Doppler model.
 *  Every different calculation method for proper time rate must be implemented in a dedicated derived class.
 */
class DopplerProperTimeRateInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param computationPointLinkEndType Variable that denotes for which link end this object computes the proper time rate
     */
    DopplerProperTimeRateInterface(
            const LinkEndType computationPointLinkEndType ):
        computationPointLinkEndType_( computationPointLinkEndType ){ }

    //! Destructor
    virtual ~DopplerProperTimeRateInterface( ){ }

    //! Function (pure virtual) to compute the proper time rate
    /*!
     *  Function (pure virtual) to compute the proper time rate as (dtau/dt-1), with tau and t proper and coordinate time,
     *  respectively.
     *  \param linkEndTimes Link end coordinate times for one-way Doppler observale for which proper time rate is to be computed.
     *  \param linkEndStates Link end Cartesian states for one-way Doppler observale for which proper time rate is to be computed.
     *  \return Proper time rate difference as (dtau/dt-1).
     */
    virtual double getOberverProperTimeDeviation(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates ) = 0;

    //! Function to retrieve variable that denotes for which link end this object computes the proper time rate
    /*!
     * Function to retrieve variable that denotes for which link end this object computes the proper time rate
     * \return Variable that denotes for which link end this object computes the proper time rate
     */
    LinkEndType getComputationPointLinkEndType( )
    {
        return computationPointLinkEndType_;
    }

protected:

    //! Variable that denotes for which link end this object computes the proper time rate
    LinkEndType computationPointLinkEndType_;
};

//! Class to compute proper time rate for one-way Doppler observation using user-provided custom function
class CustomDopplerProperTimeRateInterface: public DopplerProperTimeRateInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param computationPointLinkEndType Variable that denotes for which link end this object computes the proper time rate
     * \param properTimeRateFunction Function that is used to compute proper time rate
     */
    CustomDopplerProperTimeRateInterface(
            const LinkEndType computationPointLinkEndType,
            const std::function< double( const double ) > properTimeRateFunction ):
        DopplerProperTimeRateInterface( computationPointLinkEndType ),
        properTimeRateFunction_( properTimeRateFunction )
    { }

    //! Destructor
    ~CustomDopplerProperTimeRateInterface( ){ }

    //! Function to coompute the proper time rate
    /*!
     *  Function to coompute the proper time rate as (dtau/dt-1), with tau and t proper and coordinate time,
     *  respectively.
     *  \param linkEndTimes Link end coordinate times for one-way Doppler observale for which proper time rate is to be computed.
     *  \param linkEndStates Link end Cartesian states for one-way Doppler observale for which proper time rate is to be computed.
     *  \return Proper time rate difference as (dtau/dt-1).
     */
    virtual double getOberverProperTimeDeviation(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        // Check input consistency
        if( linkEndTimes.size( ) != 2 || linkEndStates.size( ) != 2 )
        {
            throw std::runtime_error( "Error when getting custom proper time rate for Doppler data, inconsistent input" );
        }

        return properTimeRateFunction_(
                    ( this->computationPointLinkEndType_ == transmitter ) ? linkEndTimes.at( 0 ) : linkEndTimes.at( 1 ) );
    }

private:

    //!  Function that is used to compute proper time rate
    std::function< double( const double ) > properTimeRateFunction_;


};

//! Function to compute first-order approximation of proper-time rate for one-way Doppler observable
/*!
 *  Function to compute first-order approximation of proper-time rate for one-way Doppler observable. Assumes a single static mass
 *  generating the gravity field (Schwarzschild metric) and a c^-2 expansion of proper time rate function.
 */
class DirectFirstOrderDopplerProperTimeRateInterface:
        public DopplerProperTimeRateInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param computationPointLinkEndType Variable that denotes for which link end this object computes the proper time rate
     * \param gravitationalParameterFunction Function that returns the gravitational parameter of the central body.
     * \param referenceBody Name of body generating the gravity field.
     * \param referencePointLinkEndType Link end type of central body (default unidentified_link_end, meaning that the
     * central body is not one of the link ends)
     * \param referencePointStateFunction Function that returns the state of the central body as a function of time,
     * default empty, but must be provided if referencePointLinkEndType equals unidentified_link_end.
     */
    DirectFirstOrderDopplerProperTimeRateInterface(
            const LinkEndType computationPointLinkEndType,
            const std::function< double( ) > gravitationalParameterFunction,
            const std::string& referenceBody,
            const LinkEndType referencePointLinkEndType = unidentified_link_end,
            const std::function< Eigen::Vector6d( const double ) > referencePointStateFunction =
            std::function< Eigen::Vector6d( const double ) >( ) ):
        DopplerProperTimeRateInterface( computationPointLinkEndType ),
        gravitationalParameterFunction_( gravitationalParameterFunction ),
        referenceBody_( referenceBody ),
        referencePointLinkEndType_( referencePointLinkEndType ),
        referencePointStateFunction_( referencePointStateFunction )
    {
        // Check input consistency
        if( this->computationPointLinkEndType_ == referencePointLinkEndType )
        {
            throw std::runtime_error( "Error when creating DirectFirstOrderDopplerProperTimeRateInterface, input link end types must be different" );
        }
        else if( ( this->computationPointLinkEndType_ != receiver ) &&
                 ( this->computationPointLinkEndType_ != transmitter ) )
        {
            throw std::runtime_error( "Error when creating DirectFirstOrderDopplerProperTimeRateInterface, computation point must be receiver or transmitter" );
        }
        else if( ( this->computationPointLinkEndType_ != receiver ) &&
                 ( this->computationPointLinkEndType_ != transmitter ) &&
                 ( this->computationPointLinkEndType_ != unidentified_link_end ) )
        {
            throw std::runtime_error( "Error when creating DirectFirstOrderDopplerProperTimeRateInterface, reference point must be receiver, transmitter or unidentified" );
        }
        else if( ( referencePointLinkEndType == unidentified_link_end ) && ( referencePointStateFunction == nullptr ) )
        {
            throw std::runtime_error( "Error when creating DirectFirstOrderDopplerProperTimeRateInterface, reference point must have state information" );
        }
        else if( ( referencePointLinkEndType != unidentified_link_end ) && !( referencePointStateFunction == nullptr ) )
        {
            throw std::runtime_error( "Error when creating DirectFirstOrderDopplerProperTimeRateInterface, reference point must have unambiguous state information" );
        }
    }

    //! Destructor
    ~DirectFirstOrderDopplerProperTimeRateInterface( ){ }

    //! Function to coompute the proper time rate
    /*!
     *  Function to coompute the proper time rate as (dtau/dt-1), with tau and t proper and coordinate time,
     *  respectively.
     *  \param linkEndTimes Link end coordinate times for one-way Doppler observale for which proper time rate is to be computed.
     *  \param linkEndStates Link end Cartesian states for one-way Doppler observale for which proper time rate is to be computed.
     *  \return Proper time rate difference as (dtau/dt-1).
     */
    double getOberverProperTimeDeviation(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        if( linkEndTimes.size( ) != 2 || linkEndStates.size( ) != 2 )
        {
            throw std::runtime_error( "Error when getting first order direct proper time rate for Doppler data, inconsistent input" );
        }

        // Compute central body state w.r.t. computation point
        Eigen::Vector6d computationPointRelativeState =
                getComputationPointRelativeState( linkEndTimes, linkEndStates );

        // Compute proper time rate
        return relativity::calculateFirstCentralBodyProperTimeRateDifference(
                    computationPointRelativeState, gravitationalParameterFunction_( ),
                    relativity::equivalencePrincipleLpiViolationParameter );
    }

    //! Function to compute the state of the computation point w.r.t. the central body
    /*!
     * Function to compute the state of the computation point w.r.t. the central body
     * \param linkEndTimes Link end coordinate times for one-way Doppler observale for which proper time rate is to be computed.
     * \param linkEndStates Link end Cartesian states for one-way Doppler observale for which proper time rate is to be computed.
     * \return The state of the computation point w.r.t. the central body
     */
    Eigen::Vector6d  getComputationPointRelativeState(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        return ( ( this->computationPointLinkEndType_ == transmitter ) ? linkEndStates.at( 0 ) : linkEndStates.at( 1 ) ) -
                getReferencePointState( linkEndTimes, linkEndStates );
    }

    //! Function to compute the state of the central body (e.g. origin for positions and velocities in proper time calculations)
    /*!
     * Function to compute the state of the central body (e.g. origin for positions and velocities in proper time calculations)
     * \param linkEndTimes Link end coordinate times for one-way Doppler observale for which proper time rate is to be computed.
     * \param linkEndStates Link end Cartesian states for one-way Doppler observale for which proper time rate is to be computed.
     * \return The state of the central body
     */
    Eigen::Vector6d getReferencePointState(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        if( referencePointLinkEndType_ == unidentified_link_end )
        {
            return referencePointStateFunction_(
                        ( ( this->computationPointLinkEndType_ == transmitter ) ? linkEndTimes.at( 0 ) : linkEndTimes.at( 1 ) ) );
        }
        else
        {
            return ( ( referencePointLinkEndType_ == transmitter ) ? linkEndStates.at( 0 ) : linkEndStates.at( 1 ) );
        }
    }

    //! Function to retrieve central body gravitational parameter
    /*!
     * Function to retrieve central body gravitational parameter
     * \return Central body gravitational parameter
     */
    double getGravitationalParameter( )
    {
        return gravitationalParameterFunction_( );
    }

    //! Function to return the name of body generating the gravity field.
    /*!
     * Function to return the name of body generating the gravity field
     * \return Name of body generating the gravity field
     */
    std::string getCentralBody( )
    {
        return referenceBody_;
    }

private:

    //! Function that returns the gravitational parameter of the central body.
    std::function< double( ) > gravitationalParameterFunction_;

    //! Name of body generating the gravity field.
    std::string referenceBody_;

    //! Link end type of central body (unidentified_link_end if central body is not one of the link ends)
    LinkEndType referencePointLinkEndType_;

    //! Function that returns the state of the central body as a function of time.
    /*!
     *  Function that returns the state of the central body as a function of time,  must be provided if
     *  referencePointLinkEndType equals unidentified_link_end.
     */
    std::function< Eigen::Vector6d( const double ) > referencePointStateFunction_;


};

//! Computes observable the (simplified) one-way Doppler observation between two link ends, omitting proper time rates and
//! light time corrections.
/*!
 *  Computes observable the (simplified) one-way Doppler observation between two link ends, omitting proper time rates and
 *  light time corrections. The observable is defined as d f_{B}/d_f_{A} - 1, with A the transmitter and B the receiver, and
 *  f the frequency of the signal.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayDopplerObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     *  \param transmitterProperTimeRateFunction Function to compute derivative of deviation between proper and coordinate time
     *  at transmitter, w.r.t. coordinate time.
     *  \param receiverProperTimeRateFunction Function to compute derivative of deviation between proper and coordinate time
     *  at receiver w.r.t. coordinate time.
     */
    OneWayDopplerObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
            lightTimeCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::function< ObservationScalarType( const TimeType ) > transmitterProperTimeRateFunction
            = std::function< ObservationScalarType( const TimeType ) >( ),
            const std::function< ObservationScalarType( const TimeType ) > receiverProperTimeRateFunction
            = std::function< ObservationScalarType( const TimeType ) >( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >( one_way_doppler, linkEnds, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator ),
        transmitterProperTimeRateCalculator_(
            ( transmitterProperTimeRateFunction == nullptr ) ?
                nullptr : std::make_shared< CustomDopplerProperTimeRateInterface >(
                transmitter, transmitterProperTimeRateFunction ) ),
        receiverProperTimeRateCalculator_(
            ( receiverProperTimeRateFunction == nullptr ) ?
                nullptr : std::make_shared< CustomDopplerProperTimeRateInterface >(
                receiver, receiverProperTimeRateFunction ) )
    {
        one_ = mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 );
        taylorSeriesExpansionOrder_ = 3;
    }

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param transmitterProperTimeRateCalculator Object to compute derivative of deviation between proper and coordinate time
     *  at transmitter, w.r.t. coordinate time.
     *  \param receiverProperTimeRateFunction Object to compute derivative of deviation between proper and coordinate time
     *  at receiver w.r.t. coordinate time.
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    OneWayDopplerObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
            lightTimeCalculator,
            const std::shared_ptr< DopplerProperTimeRateInterface > transmitterProperTimeRateCalculator,
            const std::shared_ptr< DopplerProperTimeRateInterface > receiverProperTimeRateFunction,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( one_way_doppler, linkEnds, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator ),
        transmitterProperTimeRateCalculator_( transmitterProperTimeRateCalculator ),
        receiverProperTimeRateCalculator_( receiverProperTimeRateFunction )
    {
        if( ( transmitterProperTimeRateCalculator == nullptr ) || (
                    receiverProperTimeRateFunction == nullptr ) )
//        {
//            throw std::runtime_error( "Error when making one-way Doppler model, input proper time rates are zero" );
//        }
        one_ = mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 );
        taylorSeriesExpansionOrder_ = 3;
    }

    //! Destructor
    ~OneWayDopplerObservationModel( ){ }


    //! Function to compute one-way Doppler observable without any corrections.
    /*!
     *  Function to compute one-way Doppler  observable without any corrections, i.e. the true physical Doppler as computed
     *  from the defined link ends. It does not include system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in double precision. These states and times are returned by
     *  reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way Doppler observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        ObservationScalarType lightTime = TUDAT_NAN;
        TimeType transmissionTime = TUDAT_NAN, receptionTime = TUDAT_NAN;

        // Compute light time
        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState_, transmitterState_, time, true );
            transmissionTime = time - lightTime;
            receptionTime = time;
            break;

        case transmitter:
            lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState_, transmitterState_, time, false );
            transmissionTime = time;
            receptionTime = time + lightTime;
            break;
        default:
            throw std::runtime_error(
                        "Error when calculating one way Doppler observation, link end is not transmitter or receiver" );
        }

        linkEndTimes.clear( );
        linkEndStates.clear( );

        // Save link end times and states
        linkEndTimes.push_back( transmissionTime );
        linkEndTimes.push_back( receptionTime );

        linkEndStates.push_back( transmitterState_.template cast< double >( ) );
        linkEndStates.push_back( receiverState_.template cast< double >( ) );

        // Compute transmitter and receiver proper time rate
        ObservationScalarType transmitterProperTimeDifference =
                mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
        if( transmitterProperTimeRateCalculator_ != nullptr )
        {
            transmitterProperTimeDifference = static_cast< ObservationScalarType >(
                        transmitterProperTimeRateCalculator_->getOberverProperTimeDeviation(
                        linkEndTimes, linkEndStates ) );
        }
        ObservationScalarType receiverProperTimeDifference =
                mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
        if( receiverProperTimeRateCalculator_ != nullptr )
        {
            receiverProperTimeDifference =  static_cast< ObservationScalarType >(
                        receiverProperTimeRateCalculator_->getOberverProperTimeDeviation(
                        linkEndTimes, linkEndStates ) );
        }

        // Compute proper time correction term
        ObservationScalarType properTimeCorrectionTerm =
                computeDopplerProperTimeInfluenceTaylorSeriesExpansion(
                    transmitterProperTimeDifference, receiverProperTimeDifference, taylorSeriesExpansionOrder_ );

        // Compute first-order (geometrical) one-way Doppler contribution
        lightTimePartialWrtReceiverPosition_ =
                lightTimeCalculator_->getPartialOfLightTimeWrtLinkEndPosition(
                    transmitterState_, receiverState_, transmissionTime, receptionTime, true );
        lightTimePartialWrtTransmitterPosition_ =
                lightTimeCalculator_->getPartialOfLightTimeWrtLinkEndPosition(
                    transmitterState_, receiverState_, transmissionTime, receptionTime, false );
        ObservationScalarType firstOrderDopplerObservable =
                computeOneWayFirstOrderDopplerTaylorSeriesExpansion<
                ObservationScalarType >(
                    transmitterState_, receiverState_,
                    lightTimePartialWrtTransmitterPosition_, lightTimePartialWrtReceiverPosition_,
                    taylorSeriesExpansionOrder_ );

        // Compute full Doppler observable and return
        ObservationScalarType totalDopplerObservable = firstOrderDopplerObservable *
                ( mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) + properTimeCorrectionTerm ) +
                properTimeCorrectionTerm;
        return ( Eigen::Matrix<  ObservationScalarType, 1, 1  >( ) << totalDopplerObservable ).finished( );
    }

    //! Function to return the object to calculate light time.
    /*!
     * Function to return the object to calculate light time.
     * \return Object to calculate light time.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

    //! Function to retrieve object to compute derivative of deviation between proper and coordinate time at transmitter
    /*!
     *  Function to retrieve object to compute derivative of deviation between proper and coordinate time at transmitter
     * \return Object to compute derivative of deviation between proper and coordinate time at transmitter
     */
    std::shared_ptr< DopplerProperTimeRateInterface > getTransmitterProperTimeRateCalculator( )
    {
        return transmitterProperTimeRateCalculator_;
    }

    //! Function to retrieve object to compute derivative of deviation between proper and coordinate time at receiver
    /*!
     *  Function to retrieve object to compute derivative of deviation between proper and coordinate time at receiver
     * \return Object to compute derivative of deviation between proper and coordinate time at receiver
     */
    std::shared_ptr< DopplerProperTimeRateInterface > getReceiverProperTimeRateCalculator( )
    {
        return receiverProperTimeRateCalculator_;
    }



private:

    //! Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator_;

    //! Templated precision value of 1.0
    ObservationScalarType one_;

    //! Order to which Doppler effect Taylor series is to be expanded.
    int taylorSeriesExpansionOrder_;

    //! Pre-declared receiver state, to prevent many (de-)allocations
    StateType receiverState_;

    //! Pre-declared transmitter state, to prevent many (de-)allocations
    StateType transmitterState_;

    //! Object to compute derivative of deviation between proper and coordinate time at transmitter, w.r.t. coordinate time.
    std::shared_ptr< DopplerProperTimeRateInterface > transmitterProperTimeRateCalculator_;

    //! Object to compute derivative of deviation between proper and coordinate time at receiver, w.r.t. coordinate time.
    std::shared_ptr< DopplerProperTimeRateInterface > receiverProperTimeRateCalculator_;

    //! Pre-declared light-time partial w.r.t. receiver sensitivity (used fopr first-order Doppler)
    Eigen::Matrix< ObservationScalarType, 1, 3 > lightTimePartialWrtReceiverPosition_;

    //! Pre-declared light-time partial w.r.t. transmitter sensitivity (used fopr first-order Doppler)
    Eigen::Matrix< ObservationScalarType, 1, 3 > lightTimePartialWrtTransmitterPosition_;

};

}

}

#endif // TUDAT_ONEWAYDOPPLEROBSERVATIONMODEL_H
