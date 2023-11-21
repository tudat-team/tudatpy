/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "tudat/astro/ephemerides/frameManager.h"


namespace tudat
{

namespace ephemerides
{

std::string getBaseFrameName( )
{
    return "SSB";
}

//! Function to determine if a given frame is an inertial frame.
bool isFrameInertial( const std::string& frame )
{
    bool isFrameInertial_;
    if( frame == "SSB" || frame == "" || frame == "Inertial" )
    {
        isFrameInertial_ = true;
    }
    else
    {
        isFrameInertial_ = false;
    }
    return isFrameInertial_;
}

//! Constructor from named list of ephemerides.
ReferenceFrameManager::ReferenceFrameManager(
        const std::map< std::string, std::shared_ptr< Ephemeris > >& ephemerisMap )
{
    // Set name of global base frame.
    frameIndexList_[ getBaseFrameName( ) ] = -1;

    // Initialize frame book-keeping variables.
    setEphemerides( ephemerisMap );
}

//! Function to determine frame levels and base frames of all frames.
void ReferenceFrameManager::setEphemerides(
        const std::map< std::string, std::shared_ptr< Ephemeris > >& additionalEphemerides )
{
    // Set list of all frames for which frame level and base frame have not yet been determined.
    std::map< std::string, std::shared_ptr< Ephemeris > > unhandledFrames_ = additionalEphemerides;

    // Set list of available ephemerides and check whether it already exists (not possible incurrent
    // implementation)
    for( std::map< std::string, std::shared_ptr< Ephemeris > >::const_iterator ephemerisIterator =
         additionalEphemerides.begin( ); ephemerisIterator != additionalEphemerides.end( );
         ephemerisIterator++ )
    {
        if( availableEphemerides_.count( ephemerisIterator->first ) == 0 )
        {
            availableEphemerides_[ ephemerisIterator->first ] = ephemerisIterator->second;
        }
        else
        {
            throw std::runtime_error( "Error when adding ephemerides to frame manager, frame already exists" );
        }
    }

    // Start at frame level 0
    int currentLevel = 0;
    std::map< std::string, std::string > singleLevelList;

    // While there are unhandled frames, continue at next frame level.
    while( unhandledFrames_.size( ) != 0 )
    {
        // Clear frame list for new level.
        singleLevelList.clear( );

        // If frame level is 0, base frame equals baseFrameName.
        if( currentLevel == 0 )
        {
            for( std::map< std::string, std::shared_ptr< Ephemeris > >::iterator
                         frameIterator = unhandledFrames_.begin( );
                 frameIterator != unhandledFrames_.end( ); frameIterator++ )
            {
                // If current frame is at level 0, add to list of levels.
                if( frameIterator->second->getReferenceFrameOrigin( ) == getBaseFrameName( ) )
                {
                    singleLevelList[ frameIterator->first ] = getBaseFrameName( );
                }
            }
        }
        // Else check if frame is at current level.
        else
        {
            std::map< std::string, std::string >::iterator previousLevelIterator;
            for( std::map< std::string, std::shared_ptr< Ephemeris > >::iterator
                         frameIterator = unhandledFrames_.begin( );
                 frameIterator != unhandledFrames_.end( ); frameIterator++ )
            {
                // If base frame of current ephemeris is on previous level, add to list of current
                // level.
                previousLevelIterator = baseFrameList_[ currentLevel - 1 ].find(
                            frameIterator->second->getReferenceFrameOrigin( ) );
                if( previousLevelIterator != baseFrameList_[ currentLevel - 1 ].end( ) )
                {
                    singleLevelList[ frameIterator->first ]
                            = frameIterator->second->getReferenceFrameOrigin( );
                }
            }
        }

        // Go through all ephemerides set at current level, add to frame list and remove from
        // unhandled frame list if present.
        std::map< std::string, std::shared_ptr< Ephemeris > >::iterator frameIterator;
        for( std::map< std::string, std::string >::iterator singleListIterator = singleLevelList.begin( );
             singleListIterator != singleLevelList.end( ); singleListIterator++ )
        {
            // Find current frame in unhandles frame list.
            frameIterator = unhandledFrames_.find( singleListIterator->first );
            if( frameIterator != unhandledFrames_.end( ) )
            {
                unhandledFrames_.erase( frameIterator );
                frameIndexList_[ singleListIterator->first ] = currentLevel;
            }
            else
            {
                throw std::runtime_error( "Error when making frame manager, could not find frame " +
                                          singleListIterator->first + " when deleting" );
            }
        }

        // If no new frame matches found and unhandled frame list is not empty: stuck in infinite loop.
        if( singleLevelList.size( ) == 0 )
        {
            throw std::runtime_error( "Error, found no new frame matches at current level, but list is not empty" );
        }

        // Add to base frame list.
        baseFrameList_.push_back( singleLevelList );
        currentLevel++;
    }

    // Check if all frames have same orientation.
    std::string firstFrameOrientation = availableEphemerides_.begin( )->second->getReferenceFrameOrientation( );
    for( std::map< std::string, std::shared_ptr< Ephemeris > >::iterator ephemerisIterator =
         availableEphemerides_.begin( ); ephemerisIterator != availableEphemerides_.end( ); ephemerisIterator++ )
    {
        if( ephemerisIterator->second->getReferenceFrameOrientation( ) != firstFrameOrientation )
        {
            throw std::runtime_error(
                        "Error, multiple reference frame orientations of ephemerides currently not supported" +
                        firstFrameOrientation + ", " + ephemerisIterator->second->getReferenceFrameOrientation( ) );
        }
    }
}


//! Returns an ephemeris along a single line of the hierarchy tree.
std::vector< std::shared_ptr< Ephemeris > > ReferenceFrameManager::getDirectEphemerisFromLowerToUpperFrame(
        const std::string& lowerFrame, const std::string& upperFrame )
{
    // Get indices of frames.
    int upperIndex = frameIndexList_.at( upperFrame );
    int lowerIndex = frameIndexList_.at( lowerFrame );

    std::vector< std::shared_ptr< Ephemeris > > ephemerisList;

    // Check validity of input (i.e. upper > lower)
    if( upperIndex < lowerIndex )
    {
        throw std::runtime_error(
            "Error when making direct ephemeris link in frame manager, upper index is smaller than lower index" );
    }
    // If frames are not equal, make list of ephemeris
    else if( upperIndex != lowerIndex )
    {
        // Start list creation at upper frame.
        int currentIndex = upperIndex;
        std::string currentFrame = upperFrame;

        // Continue while lower frame is reached
        while( currentIndex > lowerIndex )
        {
            // Check whether current frame is consistent with current frame index.
            if( frameIndexList_[ currentFrame ] != currentIndex )
            {
                throw std::runtime_error(
                            "Error when making direct constituent ephemeris, frame index inconsistent." );
            }

            // Add base frame of current frame.
            std::string currentBase = baseFrameList_[ currentIndex ][ currentFrame ];
            ephemerisList.push_back( availableEphemerides_[ currentFrame ] );

            // Decrement frame level and move to one frame level lower.
            currentIndex--;
            currentFrame = currentBase;
        }
    }
    return ephemerisList;
}

//! Return the level at which the requested ephemeris is in the hierarchy.
std::pair< int, bool > ReferenceFrameManager::getFrameLevel( const std::string& frame )
{
    std::pair< int, bool > returnValue;

    // Find level of frame.
    std::map< std::string, int >::iterator frameIndexIterator;
    frameIndexIterator = frameIndexList_.find( frame );

    // If not found, return false and NaN frame index.
    if( frameIndexIterator  == frameIndexList_.end( ) )
    {
        returnValue = std::make_pair( TUDAT_NAN, 0 );
    }
    // If found, return found value.
    else
    {
        returnValue = std::make_pair( frameIndexIterator->second, 1 );
    }
    return returnValue;
}

//! Returns the nearest common frame between frames.
std::pair< std::string, int > ReferenceFrameManager::getNearestCommonFrame( std::vector< std::string > frameList )
{
    // Make lists of current frames and frame levels that are currently being investigated for each
    // frame in frameList.
    std::vector< int > currentFrameLevels;
    currentFrameLevels.resize( frameList.size( ) );
    std::vector< std::string > currentFrames;
    currentFrames.resize( frameList.size( ) );

    // Initialize current maximum frame level to nonsense value.
    int currentMaximumFrameLevel = -1E8;

    // Check frame level of all input frames and find maximum level.
    for( unsigned int i = 0; i < frameList.size( ); i++ )
    {
        // Get current frame level.
        currentFrameLevels[ i ] = frameIndexList_[ frameList[ i ] ];
        currentFrames[ i ] = frameList[ i ];

        if( currentFrameLevels[ i ] > currentMaximumFrameLevel )
        {
            currentMaximumFrameLevel = currentFrameLevels[ i ];
        }
    }

    bool isConverged = false;
    std::string firstFrame;
    bool areAllFramesEqual = true;

    // Check if all frames are equal, i.e. if already converged.
    firstFrame = currentFrames[ 0 ];
    for( unsigned int i = 1; i < currentFrames.size( ); i++ )
    {
        if( currentFrames[ i ] != firstFrame )
        {
            areAllFramesEqual = false;
        }
    }
    if( areAllFramesEqual )
    {
        isConverged = true;
    }

    // Decrement maximum frame levels and associated frames until all frames are equal.
    while( isConverged == false && currentMaximumFrameLevel >= 0 )
    {
        // Iterate over all frames.
        for( unsigned int i = 0; i < frameList.size( ); i++ )
        {
            // If frame level of current frame is at maximum, reduce it by one level.
            if( currentFrameLevels[ i ] == currentMaximumFrameLevel )
            {
                // Find base frame of current level.
                currentFrames[ i ] = baseFrameList_[ currentFrameLevels[ i ] ][ currentFrames[ i ] ];

                // Reduce level by 1
                currentFrameLevels[ i ]--;
            }
        }
        currentMaximumFrameLevel--;

        // Check if all frames are equal.
        areAllFramesEqual = true;
        firstFrame = currentFrames[ 0 ];
        for( unsigned int i = 1; i < currentFrames.size( ); i++ )
        {
            if( currentFrames[ i ] != firstFrame )
            {
                areAllFramesEqual = false;
            }
        }

        if( areAllFramesEqual )
        {
            isConverged = true;
        }
    }

    // If not converged, set global base frame as common frame.
    std::pair< std::string, int > output;
    if( isConverged == false )
    {
        output = std::make_pair( getBaseFrameName( ), -1 );
    }
    else
    {
        output = std::make_pair( currentFrames[ 0 ], currentMaximumFrameLevel );
    }

    return output;
}

//! Get name of the base frame for a given body.
std::string ReferenceFrameManager::getBaseFrameNameOfBody( const std::string& bodyName )
{
    std::string frameName;

    if( availableEphemerides_.count( bodyName ) == 0 )
    {
        throw std::runtime_error(
                    "Error when getting base frame name of body, body " + bodyName + " not found to have an ephemeris" );
    }
    else
    {
        frameName = availableEphemerides_.at( bodyName )-> getReferenceFrameOrigin( );
    }
    return frameName;
}

//! Get ephemeris origins (base frame names) for a list of bodies.
std::vector< std::string > ReferenceFrameManager::getEphemerisOrigins(
        const std::vector< std::string >& bodyList )
{
    std::vector< std::string > ephemerisOrigins;
    for( unsigned int i = 0; i < bodyList.size( ); i++ )
    {
        ephemerisOrigins.push_back( getBaseFrameNameOfBody( bodyList.at( i ) ) );
    }
    return ephemerisOrigins;
}


} // namespace ephemerides

} // namespace tudat
