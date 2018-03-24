/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_FRAMEMANAGER_H
#define TUDAT_FRAMEMANAGER_H

#include <map>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Function to determine if a given frame is an inertial frame.
/*!
 *  Function to determine if a given frame is an inertial frame. Currently a frame identified as
 *  "SSB" (solar system barycenter), "inertial" or "" (empty) is recognized as inertial.
 *  \param frame Name of frame for which it is to be determined whether it is inertial.
 *  \return True if inertial, false if not.
 */
bool isFrameInertial( const std::string& frame );

//! Function to return base frame
/*!
 * Function to return base frame, i.e. in which the states of the bodies are defined during the integration.
 * \return Base frame for simulations.
 */
std::string getBaseFrameName( );


//! Class to retrieve translation functions between different frames
/*!
 * Class to retrieve translation functions between different frames, as calculated from a list of
 * ephemeris objects.  Using this class, the various Ephemeris objects may be 'pasted' together to
 * obtain the state of one body w.r.t. any other body.
 */
class ReferenceFrameManager
{
public:

    //! Constructor from named list of ephemerides.
    /*!
     *  Constructor from named list of ephemerides.
     *  \param ephemerisMap List of ephemerides per body.
     */
    ReferenceFrameManager( const std::map< std::string, boost::shared_ptr< Ephemeris > >& ephemerisMap );

    //! Function to retrieve the ephemeris of a body with a requested frame origin.
    /*!
     *  Function to retrieve the ephemeris of a body with a requested frame origin. Both the body
     *  and the origin must be loaded into the frame manager.

     *  \param origin Origin of ephemeris
     *  \param body Body for which ephemeris is requested.
     *  \return Ephemeris of requested body qith requested frame origin
     */
    template< typename StateScalarType = double, typename TimeType = double >
    boost::shared_ptr< Ephemeris > getEphemeris(
            const std::string& origin, const std::string& body )
    {
        typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;
        boost::shared_ptr< Ephemeris > ephemerisBetweenFrames;

        // If requested 'body' is global base frame, return constant zero ephemeris.
        if( body == origin )
        {
            //NOTE: Should generalize to long double state type.
            ephemerisBetweenFrames = boost::make_shared< ConstantEphemeris >(
                        boost::lambda::constant( Eigen::Vector6d::Zero( ) ), origin, "ECLIPJ2000" );
        }
        else
        {
            // Set frames between which ephemeris must be created.
            std::vector< std::string > framesToCheck;
            framesToCheck.push_back( origin );
            framesToCheck.push_back( body );

            // Find nearest common frame between frames.
            std::pair< std::string, int > nearestCommonFrame = getNearestCommonFrame( framesToCheck );

            // Initialize list of ephemeris functions for composite ephemeris creation
            std::map< int, std::pair< boost::function< StateType( const TimeType& ) >, bool > >
                 totalEphemerisList;

            // If body is nearest common frame, get set of ephemeris and set to subtract them when
            // making composite ephemeris.
            std::vector< boost::shared_ptr< Ephemeris > > ephemerisList;
            if( nearestCommonFrame.first == body )
            {
                ephemerisList = getDirectEphemerisFromLowerToUpperFrame( body, origin );
                for( unsigned int i = 0; i < ephemerisList.size( ); i++ )
                {
                    totalEphemerisList[ i ] = std::make_pair(
                                boost::bind( &Ephemeris::getTemplatedStateFromEphemeris
                                             < StateScalarType, TimeType >,
                                             ephemerisList[ i ], _1 ), false );
                }
            }
            // If origin is nearest common frame, get set of ephemeris and set to add them when
            // making composite ephemeris.
            else if( nearestCommonFrame.first == origin )
            {
                ephemerisList = getDirectEphemerisFromLowerToUpperFrame( origin, body );
                for( unsigned int i = 0; i < ephemerisList.size( ); i++ )
                {
                    totalEphemerisList[ i ] = std::make_pair(
                                boost::bind( &Ephemeris::getTemplatedStateFromEphemeris
                                             < StateScalarType, TimeType >,
                                             ephemerisList[ i ], _1 ), true );
                }
            }
            // If nearest common frame is neither input, create link from both to nearest common frame.
            else
            {
                // Get set of ephemeris from nearest common frame to body and set to add them when
                // making composite ephemeris.
                ephemerisList = getDirectEphemerisFromLowerToUpperFrame( nearestCommonFrame.first, body );
                for( unsigned int i = 0; i < ephemerisList.size( ); i++ )
                {
                    totalEphemerisList[ i ] = std::make_pair(
                                boost::bind( &Ephemeris::getTemplatedStateFromEphemeris
                                             < StateScalarType, TimeType >,
                                             ephemerisList[ i ], _1 ), true );
                }
                int firstListSize = ephemerisList.size( );

                // Get set of ephemeris from nearest common frame to origin and set to subtract them
                // when making composite ephemeris.
                ephemerisList = getDirectEphemerisFromLowerToUpperFrame( nearestCommonFrame.first, origin );
                for( unsigned int i = 0; i < ephemerisList.size( ); i++ )
                {
                    totalEphemerisList[ i + firstListSize ] = std::make_pair(
                                boost::bind( &Ephemeris::getTemplatedStateFromEphemeris
                                             < StateScalarType, TimeType >,
                                             ephemerisList[ i ], _1 ), false );
                }
            }

            // Create composite ephemeris
            ephemerisBetweenFrames = boost::make_shared< CompositeEphemeris< TimeType, StateScalarType > >(
                        totalEphemerisList,
                        std::map< int, boost::function< StateType( const TimeType, const StateType& ) > >( ),
                        origin );
        }

        return ephemerisBetweenFrames;
    }

    //! Return the level at which the requested ephemeris is in the hierarchy.
    /*!
     *  Return the level at which the requested ephemeris is in the hierarchy.
     *  \param frame Frame for which the frame level is requested.
     *  \return Pair of frame level and boolean. Boolean is true if frame exists,
     *          false if not (and frame level NAN).
     */
    std::pair< int, bool > getFrameLevel( const std::string& frame );

    //! Returns the nearest common frame between frames.
    /*!
     *  Returns the nearest common frame between frames.
     *  \param frameList Vector of frame names for which the nearest common frame is to be found.
     *  \return Nearest common frame, with its frame level.
     */
    std::pair< std::string, int > getNearestCommonFrame( std::vector< std::string > frameList );

    //! Get name of the base frame for a given body.
    /*!
     * Get name of the base frame for a given body.
     * \param bodyName Name of body for which the base frame is to be returned.
     * \return Name od base frame of bodyName.
     */
    std::string getBaseFrameNameOfBody( const std::string& bodyName );

    //! Get ephemeris origins (base frame names) for a list of bodies.
    /*!
     * Get ephemeris origins (base frame names) for a list of bodies, calls getBaseFrameNameOfBody
     * for each entry of bodyList.
     * \param bodyList List of bodies for which ephemeris origins are to be calculated.
     * \return Ephemeris origins of bodyList
     */
    std::vector< std::string > getEphemerisOrigins( const std::vector< std::string >& bodyList );

private:

    //! Vector of frames with associated base frames, ordered by frame level.
    /*!
     *  Vector of frames with associated base frames, ordered by frame level. Index of vector
     *  indicates at which frame level the frames in the maps at that index are. The map contains
     *  frame names (key) with their associated base frames (value). Each base frame must be a frame
     *  of one level lower; base frame of level 0 frames must be global base frame.
     */
    std::vector< std::map< std::string, std::string > > baseFrameList_;

    //! Map of ephemerides (values) of bodies (keys)
    /*!
     *  Map of ephemerides (values) of bodies (keys). The frame origins and orientations of these
     *  ephemerides can be retrievd from the objects themselves.
     */
    std::map< std::string, boost::shared_ptr< Ephemeris > > availableEphemerides_;

    //! Map giving the frame level for each frame name.
    /*!
     *  Map giving the frame level for each frame name.
     */
    std::map< std::string, int > frameIndexList_;

    //! Returns an ephemeris along a single line of the hierarchy tree.
    /*!
     *  Returns an ephemeris along a single line of the hierarchy tree, i.e. returned ephemeris
     *  constituent frame levels must be continuously increasing
     */
    std::vector< boost::shared_ptr< Ephemeris > > getDirectEphemerisFromLowerToUpperFrame(
            const std::string& lowerFrame, const std::string& upperFrame );

    //! Function to determine frame levels and base frames of all frames.
    /*!
     *  Function to determine frame levels and base frames of all frames; called by constructor.
     */
    void setEphemerides( const std::map< std::string,
                         boost::shared_ptr< Ephemeris > >& additionalEphemerides );

};

template< typename StateScalarType = double, typename TimeType = double >
//! Function to get a list of translation functions from integration frames to ephemeris frames.
/*!
 * Function to get a list of translation functions from integration frames to ephemeris frames. The
 * output functions provide the state of the ephemeris origin in the integration origin.
 * \param centralBodies List of integration origins.
 * \param bodiesToIntegrate List of bodies for which the origins are considered.
 * \param frameManager Object to retrieve translations between origins
 * \return List of translation functions from integration frames to ephemeris frames.
 */
std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >
getTranslationFunctionsFromIntegrationFrameToEphemerisFrame(
        const std::vector< std::string >& centralBodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager )
{
    std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >
            translationFunctionMap;

    // Check consistency of input.
    if( centralBodies.size( ) != bodiesToIntegrate.size( ) )
    {
        throw std::runtime_error(
         "Error when making translation functrions from integration to ephemeris frame, input vector sizes inconsistent" );
    }
    else
    {
        for( unsigned int i = 0; i < centralBodies.size( ); i++ )
        {
            // If ephemeris and integration origin are not equal, provide translation function.
            if( centralBodies.at( i ) != frameManager->getBaseFrameNameOfBody( bodiesToIntegrate.at( i ) ) )
            {
                translationFunctionMap[ bodiesToIntegrate.at( i ) ] = boost::bind(
                            &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< StateScalarType, TimeType >,
                            frameManager->getEphemeris< StateScalarType, TimeType >(
                                centralBodies.at( i ),
                                frameManager->getBaseFrameNameOfBody( bodiesToIntegrate.at( i ) )  ), _1 );
            }
        }
    }

    return translationFunctionMap;
}

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_FRAMEMANAGER_H
