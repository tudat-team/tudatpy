/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INITIALTRANSLATIONALSTATE_H
#define TUDAT_INITIALTRANSLATIONALSTATE_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of an initial translational state.
template< typename InitialStateParameterType = double >
class InitialTranslationalStateParameter: public EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:

    //! Constructor, sets initial value of translational state.
    /*!
     * Constructor, sets initial value of translational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTranslationalState Current value of initial state (w.r.t. centralBody)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateParameter(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialTranslationalState,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >(
            initial_body_state, associatedBody ),
        initialTranslationalState_( initialTranslationalState ), centralBody_( centralBody ),
        frameOrientation_( frameOrientation )
    { }

    //! Function to get the current value of initial state w.r.t. centralBody.
    /*!
     * Function to get the current value of initial state w.r.t. centralBody.
     * \return The current value of initial state w.r.t. centralBody.
     */
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialTranslationalState_;
    }

    //! Function to reset the current value of initial state w.r.t. centralBody.
    /*!
     * Function to reset the current value of initial state w.r.t. centralBody.
     * \param parameterValue The new value of initial state w.r.t. centralBody.
     */
    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  parameterValue )
    {
        initialTranslationalState_ = parameterValue;
    }

    //! Function to retrieve the size of the parameter (always 6).
    /*!
     *  Function to retrieve the size of the parameter (always 6).
     *  \return Size of parameter value (always 6).
     */
    int getParameterSize( )
    {
        return 6;
    }

    //! Function to get the name of the body w.r.t. which the initial state is to be estimated.
    /*!
     * Function to get the name of the body w.r.t. which the initial state is to be estimated.
     * \return Name of the body w.r.t. which the initial state is to be estimated.
     */
    std::string getCentralBody( )
    {
        return centralBody_;
    }

private:

    //! Current value of initial state (w.r.t. centralBody)
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalState_;

    //!  Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

};

//! Interface class for the estimation of an arcwise initial translational state.
template< typename InitialStateParameterType = double >
class ArcWiseInitialTranslationalStateParameter: public EstimatableParameter< Eigen::Matrix<
        InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Body for which initial state is to be estimated.
     * \param arcStartTimes Start times for separate arcs
     * \param initialTranslationalState List of initial arc states
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ArcWiseInitialTranslationalStateParameter(
            const std::string& associatedBody,
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >& initialTranslationalState,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >(
            arc_wise_initial_body_state, associatedBody ),
        arcStartTimes_( arcStartTimes ), centralBody_( centralBody ), frameOrientation_( frameOrientation )
    {
        if( arcStartTimes_.size( ) != initialTranslationalState.size( ) )
        {
            throw std::runtime_error(
                        "Error when creating arc-wise initial translational state parameters, incompatible sizes are found" );
        }
        else
        {
            for( unsigned int i = 0; i < initialTranslationalState.size( ); i++ )
            {
                initialTranslationalState_.segment( i * 6, 6 ) = initialTranslationalState.at( i );
            }
        }
    }

    //! Constructor
    /*!
     * Constructor
     * \param associatedBody Body for which initial state is to be estimated.
     * \param arcStartTimes Start times for separate arcs
     * \param initialTranslationalStates Current value of initial arc states (concatenated in same order as arcs)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ArcWiseInitialTranslationalStateParameter(
            const std::string& associatedBody,
            const std::vector< double >& arcStartTimes,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalStates,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >(
            arc_wise_initial_body_state, associatedBody ),
        initialTranslationalState_( initialTranslationalStates ),
        arcStartTimes_( arcStartTimes ), centralBody_( centralBody ), frameOrientation_( frameOrientation )
    {
        if( 6 * static_cast< int >( arcStartTimes_.size( ) ) != initialTranslationalStates.rows( ) )
        {
            throw std::runtime_error(
                        "Error when creatiung arc-wise initial translational state parameters, incompatible sizes " );
        }
    }

    //! Function to get the current values of initial states w.r.t. centralBody for each arc, concatenated in arc order.
    /*!
     *Function to get the current values of initial states w.r.t. centralBody for each arc, concatenated in arc order
     * \return The current values of initial states w.r.t. centralBody for each arc, concatenated in arc order
     */
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialTranslationalState_;
    }

    //! Function to reset the current values of initial states w.r.t. centralBody for each arc, concatenated in arc order.
    /*!
     *Function to reset the current values of initial states w.r.t. centralBody for each arc, concatenated in arc order
     * \param parameterValue New values of initial states w.r.t. centralBody for each arc, concatenated in arc order
     */
    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  parameterValue )
    {
        initialTranslationalState_ = parameterValue;
    }

    //! Function to retrieve the size of the parameter (always 6 times number of arcs).
    /*!
     *  Function to retrieve the size of the parameter (always 6 times number of arcs).
     *  \return Size of parameter value (always 6 times number of arcs).
     */
    int getParameterSize( )
    {
        return 6 * arcStartTimes_.size( );
    }

    //! Function to retrieve the number of state arcs
    /*!
     *  Function to retrieve the number of state arcs
     *  \return Number of state arcs
     */
    int getNumberOfStateArcs( )
    {
        return arcStartTimes_.size( );
    }

    //! Function to get the name of the body w.r.t. which the initial states are to be estimated.
    /*!
     * Function to get the name of the body w.r.t. which the initial states are to be estimated.
     * \return Name of the body w.r.t. which the initial states are to be estimated.
     */
    std::string getCentralBody( )
    {
        return centralBody_;
    }

    //! Function to get the start times for separate arcs
    /*!
     * Function to get the start times for separate arcs
     * \return Start times for separate arcs
     */
    std::vector< double > getArcStartTimes( )
    {
        std::vector< double > arcStartTimes;
        for( unsigned int i = 0; i < arcStartTimes_.size( ); i++ )
        {
            arcStartTimes.push_back( arcStartTimes_.at( i ) );
        }

        return arcStartTimes;
    }

private:

    //! The current values of initial states w.r.t. centralBody for each arc, concatenated in arc order
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalState_;

    //! Start times for separate arcs
    std::vector< double > arcStartTimes_;

    //! Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

};

//! Function to retrieve the size of the estimatable parameter set.
/*!
 *  Function to retrieve the size of the estimatable parameter set.
 *  \param estimatableParameterSet Set of estimatable parameters.
 *  \return Size of parameter set.
 */
template< typename InitialStateParameterType = double >
int getSingleArcParameterSetSize(
        boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameterSet )
{
    int totalParameterSetSize = estimatableParameterSet->getEstimatedParameterSetSize( );
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialStateParameters = estimatableParameterSet->getEstimatedInitialStateParameters( );

    for( unsigned int i = 0; i < initialStateParameters.size( ); i++ )
    {
        if( initialStateParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state )
        {
            totalParameterSetSize -=
                    ( boost::dynamic_pointer_cast< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > >(
                        initialStateParameters.at( i ) )->getNumberOfStateArcs( ) - 1 ) * 6;
        }
        else if( ( initialStateParameters.at( i )->getParameterName( ).first != initial_body_state ) )
        {
            throw std::runtime_error(
                        "Error when getting single arc paramater vector, did not recognize initial state parameter " +
                        std::to_string( initialStateParameters.at( i )->getParameterName( ).first ) );
        }
    }
    return totalParameterSetSize;
}

//! Function to retrieve the size of the dynamical state.
/*!
 *  Function to retrieve the size of the dynamical state.
 *  \param estimatableParameterSet Set of estimatable parameters.
 *  \return Size of the initial dynamical state.
 */
template< typename InitialStateParameterType = double >
int getSingleArcInitialDynamicalStateParameterSetSize(
        boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameterSet )
{
    return getSingleArcParameterSetSize( estimatableParameterSet ) -
            ( estimatableParameterSet->getEstimatedParameterSetSize( ) -
              estimatableParameterSet->getInitialDynamicalStateParameterSize( ) );
}

//! Function to get arc start times from list of estimated parameters
/*!
 *  Function to get arc start times from list of estimated parameters. Function throws an error if multiple arc-wise
 *  estimations are found, but arc times are not compatible
 *  \param estimatableParameters List of estimated parameters
 *  \return Start times for estimation arcs
 */
template< typename InitialStateParameterType >
std::vector< double > getMultiArcStateEstimationArcStartTimes(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< boost::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    std::vector< double > arcStartTimes;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state )
        {
            boost::shared_ptr< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > > arcWiseStateParameter =
            boost::dynamic_pointer_cast< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > >(
                        initialDynamicalParameters.at( i ) );
            if( arcWiseStateParameter == NULL )
            {
                throw std::runtime_error( "Error when getting arc times from estimated parameters, parameter is inconsistent" );
            }

            if( arcStartTimes.size( ) == 0 )
            {
                arcStartTimes = arcWiseStateParameter->getArcStartTimes( );
            }
            else
            {
                // Check arc consistency
                std::vector< double > currentArcStartTimes = arcWiseStateParameter->getArcStartTimes( );
                if( currentArcStartTimes.size( ) != arcStartTimes.size( ) )
                {
                    throw std::runtime_error(
                                "Error when getting arc start times for parameters, number of arcs is inconsistent" );
                }
                else
                {
                    for( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
                    {
                        if( arcStartTimes.at( j ) != currentArcStartTimes.at( j ) )
                        {
                            throw std::runtime_error(
                                        "Error when getting arc start times for parameters, arc times are inconsistent" );
                        }
                    }
                }
            }

        }
        else
        {
            throw std::runtime_error( "Error when getting arc times from estimated parameters, soingle arc dynamics found" );
        }
    }

    return arcStartTimes;
}


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_INITIALTRANSLATIONALSTATE_H
