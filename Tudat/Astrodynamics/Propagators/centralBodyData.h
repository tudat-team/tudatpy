/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CENTRALBODYDATA_H
#define TUDAT_CENTRALBODYDATA_H

#include <vector>
#include <map>

#include <boost/function.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace propagators
{

//! Types of integration origins that can be used in the simulations.
enum OriginType
{
    global_frame_origin,// origin is inertial.
    from_ephemeris, // origin is moving with origin of body which is not propagated.
    from_integration // origin is moving with origin of body which is propagated.
};

//! This class acts as a data container for properties of central bodies in a numerical integration
/*!
 *  This class acts as a data container for properties of central bodies in a numerical integration.
 *  It makes a distinction between central bodies that are integrated and those for which the state
 *  is taken from ephemeris The state of the central bodies in an inertial frame can be retrieved
 *  from it.
 */
template< typename StateScalarType = double, typename TimeType = double >
class CentralBodyData
{
public:
    //! Constructor of the central bodies and the integrated bodies
    /*!
     *  Constructor, takes the names of the central bodies and the integrated bodies, as well as a
     *  list of body objects.  Checks which central bodies are integrated bodies and sets the update
     *  order of the bodies' states accordingly.
     *  \param centralBodies Names of central bodies, belonging to the entries in the
     *         bodiesToIntegrate vector of same index.
     *  \param bodiesToIntegrate Names of bodies that are to be integrated numerically.
     *  \param bodyStateFunctions List of functions for the origins of selected bodies.
     *  \param globalFrameOriginBarycentricStateFunction State function of global frame origin w.r.t. barycenter
     *  \param globalFrameOrigin Origin of global frame.
     */
    CentralBodyData( const std::vector< std::string >& centralBodies,
                     const std::vector< std::string >& bodiesToIntegrate,
                     const std::map< std::string,
                     boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >&
                     bodyStateFunctions,
                     const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >
                     globalFrameOriginBarycentricStateFunction,
                     const std::string globalFrameOrigin ):
        centralBodies_( centralBodies ), globalFrameOriginBarycentricStateFunction_( globalFrameOriginBarycentricStateFunction )
    {
        // Check consistency of input.
        if( centralBodies.size( ) != bodiesToIntegrate.size( ) )
        {
            throw std::runtime_error(
                "Error in CentralBodyData, number of central bodies not equal to number of bodies to integrate " );
        }


        bodyOriginType_.resize( bodiesToIntegrate.size( ) );

        // Iterate over all integrated bodies to set body origin, its type and associated data.
        for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
        {
            // Check if central body is inertial.
            if( centralBodies.at( i ) == globalFrameOrigin )
            {
                bodyOriginType_.at( i ) = global_frame_origin;
            }
            else
            {
                // Check if central body of current integrated body is also integrated.
                int centralBodyIndex = -1;
                for( unsigned int j = 0; j < bodiesToIntegrate.size( ); j++ )
                {
                    // If there is a match between central body and an integrated body, store body
                    // index.
                    if( centralBodies.at( i ) == bodiesToIntegrate[ j ] )
                    {
                        centralBodyIndex = j;
                        if( i == j )
                        {
                            throw std::runtime_error( "Error, body " + bodiesToIntegrate[ j ] +
                                                      " cannot be its own central body" );
                        }
                    }
                }

                // If no integrated central body foundl set body origin as being from an ephemeris.
                if( centralBodyIndex == -1 )
                {
                    bodyOriginType_.at( i ) = from_ephemeris;
                    centralBodiesFromEphemerides_[ i ]
                        = bodyStateFunctions.at( centralBodies.at( i ) );
                }
                // Else, set body origin as being from another integrated body, set indices of both
                // bodies in centralBodiesFromIntegration_
                else
                {
                    bodyOriginType_.at( i ) = from_integration;
                    centralBodiesFromIntegration_[ i ] = centralBodyIndex;
                }
            }
        }

        std::vector< int > numericalBodies_;

        // Set initial update order of bodies.
        updateOrder_.resize( bodiesToIntegrate.size( ) );
        int currentUpdateIndex = 0;
        for( unsigned int i = 0; i < bodyOriginType_.size( ); i++ )
        {
            // If body origin is not from another integrated body, order in update is irrelevant,
            // set at current index of vector.
            if( bodyOriginType_.at( i ) == global_frame_origin || bodyOriginType_.at( i ) == from_ephemeris  )
            {
                updateOrder_[ currentUpdateIndex ] = i;
                currentUpdateIndex++;
            }
            // Else, store index of body having another integrated body as central body
            else if( bodyOriginType_.at( i ) == from_integration )
            {
                numericalBodies_.push_back( i );
            }
        }

        // Switch body order if there are dependencies between bodies,
        int bodyToMove;
        for( unsigned int i = 0; i < numericalBodies_.size( ); i++ )
        {
            for( unsigned int j = 0; j < i; j++ )
            {
                if( centralBodies[ numericalBodies_[ j ] ]
                    == bodiesToIntegrate[ numericalBodies_.at( i ) ] )
                {
                    // Move central body to index before integrated body.
                    bodyToMove = numericalBodies_.at( i );
                    numericalBodies_.erase( numericalBodies_.begin( ) + i );
                    numericalBodies_.insert( numericalBodies_.begin( ) + j, bodyToMove );
                    break;
                }
            }
        }

        // Add propagated bodies to list.
        for( unsigned int i = 0; i < numericalBodies_.size( ); i++ )
        {
            updateOrder_[ currentUpdateIndex ] = numericalBodies_.at( i );
            currentUpdateIndex++;
        }

        localInternalState_.resize( 6 * bodiesToIntegrate.size( ), 1 );
    }


    //! Function to return the state of the central bodies in an inertial frame.
    /*!
     *  Function to return the state of the central bodies in an inertial frame. The states of the
     *  integrated bodies in either their local frame or the inertial frame, and the current time
     *  have to be passed as arguments to this function.
     *  \param internalState States of bodies that are numerically integrated, size should be 6 *
     *  size of bodiesToIntegrate, with entries in the order of the bodies in the bodiesToIntegrate
     *  vector.
     *  \param time Current time (used for retrieving states from ephemerides)
     *  \param referenceFrameOriginStates Vector of states of the reference frame origins for each body
     *  (returned by reference).
     *  \param areInputStateLocal True if the internalState vector is given in the local frames of the integrated
     *   bodies, or the global frame.
     */
    void getReferenceFrameOriginInertialStates(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalState, const TimeType time,
            std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& referenceFrameOriginStates,
            const bool areInputStateLocal = true )
    {
        localInternalState_ =  internalState;
        if( referenceFrameOriginStates.size( ) != updateOrder_.size( ) )
        {
            referenceFrameOriginStates.resize( updateOrder_.size( ) );
        }

        // Update state in correct order.
        for( unsigned int i = 0; i < updateOrder_.size( ); i++ )
        {
            getSingleReferenceFrameOriginInertialState(
                        localInternalState_, time, updateOrder_.at( i ),
                        referenceFrameOriginStates.at( updateOrder_.at( i ) ));

            // Modify current input state to global frame if input is local (in propagation frame).
            if( areInputStateLocal )
            {
                localInternalState_.segment( 6 * updateOrder_.at( i ), 6 ) +=
                        referenceFrameOriginStates.at( updateOrder_.at( i ) );
            }
        }
    }


    //! Function to get the order in which the body states are to be called
    /*!
     * Function to get the order in which the body states are to be called when getting global
     * states.
     * \return Order in which the body states are to be called when getting global states.
     */
    std::vector< int > getUpdateOrder( ){ return updateOrder_; }

    //! Function to get the type of reference frame origin for each of the propagated bodies.
    /*!
     * Function to get the type of reference frame origin for each of the propagated bodies.
     * \return Type of reference frame origin for each of the propagated bodies.
     */
    std::vector< OriginType > getBodyOriginType( ){ return bodyOriginType_; }

    //! Function to get the names of central bodies
    /*!
     * Function to get the names of central bodies, belonging to the entries in the
     * bodiesToIntegrate vector of same index.
     * \return Names of central bodies, belonging to the entries in the bodiesToIntegrate vector of
     * same index.
     */
    std::vector< std::string > getCentralBodies( ){ return centralBodies_; }


private:

    //! Names of central bodies, belonging to the entries in the bodiesToIntegrate vector of same
    //! index.
    std::vector< std::string > centralBodies_;

    boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > globalFrameOriginBarycentricStateFunction_;

    //! Order in which the body states are to be called when getting global states (taking into
    //! account frame origin dependencies).
    std::vector< int > updateOrder_;

    //! Type of reference frame origin for each of the propagated bodies.
    std::vector< OriginType > bodyOriginType_;

    //!  List of functions for the origins of selected bodies.
    std::map< int, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >
        centralBodiesFromEphemerides_;

    //! Map defining frame origin body index, for bodies having one of the other propagated bodies
    //! as propagation origin
    std::map< int, int > centralBodiesFromIntegration_;

    //! State of propagated bodies, with the origins translated to the global origin
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > localInternalState_;

    //! Function to get the global origin of the propagation center of a single body.
    /*!
     *  Function to get the global origin of the propagation center of a single body, where the
     *  origin may be inertial, from an ephemeris, or from one of the other propagated bodies.
     *  \param internalSolution Full current state of the bodies being propagated, modified so that
     *  the origin of the current body is already translated to the global origin (if from
     *  integration)
     *  \sa getReferenceFrameOriginInertialStates
     *  \sa updateOrder
     *  \param time Current time.
     *  \param bodyIndex Index of the body for which the global origin state is to be retrieved
     *  \param originState Global origin state of the requested body (returned by reference).
     */
    void getSingleReferenceFrameOriginInertialState(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution,
            const TimeType time,
            const int bodyIndex,
            Eigen::Matrix< StateScalarType, 6, 1 >& originState )
    {
        // Check origin type.
        switch( bodyOriginType_[ bodyIndex ] )
        {
        case global_frame_origin:
            originState.setZero( );
            break;
        case from_ephemeris:
            originState = centralBodiesFromEphemerides_.at( bodyIndex )( static_cast< double >( time ) );
            break;
        case from_integration:
            originState = internalSolution.segment( centralBodiesFromIntegration_.at( bodyIndex ) * 6, 6 );
            break;
        default:
            throw std::runtime_error( "Error, do not recognize boy origin type " +
                              std::to_string( bodyOriginType_[ bodyIndex ] ) );
            break;
        }
    }
};


} // namespace propagators

} // namespace tudat
#endif // TUDAT_CENTRALBODYDATA_H
