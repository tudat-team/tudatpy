/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_COMPOSITEEPHEMERIS_H
#define TUDAT_COMPOSITEEPHEMERIS_H

#include <vector>


#include <memory>


#include <Eigen/Core>

#include "tudat/astro/ephemerides/ephemeris.h"

#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"



namespace tudat
{

namespace ephemerides
{


//! Class that combines a series of translational and rotational ephemeris functions
/*!
 *  Class that combines a series of translational and rotational ephemeris functions to yield a
 *  single translational ephemeris.  By using this class, a single object can be used for calling a
 *  series of translations and rotations.  The option is provided for both adding or subtracting a
 *  given translational state.
 */
template< typename TimeType = double, typename StateScalarType = double >
class CompositeEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianLongState;
    using Ephemeris::getCartesianState;

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    //! Constructor from series of translational and rotational ephemeris functions.
    /*!
     *  Constructor from series of translational and rotational ephemeris functions. Input is
     *  provided as maps, with key being the index providing the position in the order where the
     *  given function should be applied. Each index must be present only once (i.e. either in
     *  translation or rotation), must start with 0 (which must be a translational ephemeris), and
     *  must be increasing by 1.
     *  \param translationalEphemerides List of translational ephemerides.
     *  \param rotationalEphemerides List of rotational ephemerides.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    CompositeEphemeris(
            const std::map< int, std::function< StateType( const TimeType& ) > >
            translationalEphemerides,
            const std::map< int, std::function< StateType( const TimeType, const StateType& ) > >
            rotationalEphemerides,
            const std::string referenceFrameOrigin = "SSB",
            const std::string referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
    {
        // Create iterators over ephemeris functions.
        typename std::map< int, std::function< StateType( const TimeType& ) > >::const_iterator
                translationIterator = translationalEphemerides.begin( );
        typename std::map< int, std::function< StateType( const TimeType, const StateType& ) > >::const_iterator
                rotationIterator = rotationalEphemerides.begin( );

        // Check whether chain starts with translation.
        if( translationIterator->first != 0 )
        {
            std::string errorMessage = "Error, composite ephemeris must start with translation";
            throw std::runtime_error( errorMessage );
        }

        // Run over all indices and set order.
        int currentIndex = 0;
        while( currentIndex < static_cast< int >( translationalEphemerides.size( )
                                                  + rotationalEphemerides.size( ) ) )
        {
            // If current ephemeris is translational, add to translations list (set as addition) and
            // set translation flag at current index to true.
            if( translationIterator != translationalEphemerides.end( )
                    && translationIterator->first == currentIndex )
            {
                translationalEphemerides_.push_back( std::make_pair( translationIterator->second, 1 ) );
                isCurrentEphemerisTranslational_.push_back( 1 );
                translationIterator++;
            }
            // If current ephemeris is rotational, add to rotations list and set translation flag at
            // current index to false.
            else if( rotationIterator != rotationalEphemerides.end( )
                     && rotationIterator->first == currentIndex )
            {
                rotationalEphemerides_.push_back( rotationIterator->second );
                isCurrentEphemerisTranslational_.push_back( 0 );
                rotationIterator++;
            }
            // If index is not found, display error message.
            else
            {
                std::string errorMessage = "Error when  making composite ephemeris, input indices inconsistentn";
                throw std::runtime_error( errorMessage );
            }
            currentIndex++;
        }

    }

    //! Constructor from series of translational and rotational ephemeris functions
    /*!
     *  Constructor from series of translational and rotational ephemeris functions, allowing either
     *  addition or subtraction of translation ep Input is provided as maps, with key being the
     *  index providing the position in the order where the given function should be applied. Each
     *  index must be present only once (i.e. either in translation or rotation), must start with 0
     *  (which must be a translational ephemeris), and must be increasing by 1.  The second element
     *  of the pair that is the translational map value denotes whether to add (1) or subtract(0)
     *  it.
     *  \param translationalEphemerides List of translational ephemerides, with subtraction/addition indicator.
     *  \param rotationalEphemerides List of rotational ephemerides.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    CompositeEphemeris(
            const std::map< int, std::pair< std::function< StateType( const TimeType& ) >, bool > >
            translationalEphemerides,
            const std::map< int, std::function< StateType( const TimeType, const StateType& ) > >
            rotationalEphemerides,
            const std::string referenceFrameOrigin = "SSB",
            const std::string referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
    {
        // Create iterators over ephemeris functions.
        typename std::map< int, std::pair< std::function< StateType( const TimeType& ) >, bool > >
                ::const_iterator translationIterator = translationalEphemerides.begin( );
        typename std::map< int, std::function< StateType( const TimeType, const StateType& ) > >
                ::const_iterator rotationIterator = rotationalEphemerides.begin( );

        // Check whether chain starts with translation.
        if( translationIterator->first != 0 )
        {
            std::string errorMessage = "Error, composite ephemeris must start with translation";
            throw std::runtime_error( errorMessage );
        }

        // Run over all indices and set order.
        int currentIndex = 0;
        while( currentIndex < static_cast< int >( translationalEphemerides.size( )
                                                  + rotationalEphemerides.size( ) ) )
        {
            // If current ephemeris is translational, add to translations list and set translation
            // flag at current index to true.
            if( translationIterator != translationalEphemerides.end( )
                    && translationIterator->first == currentIndex )
            {
                int addCurrentEphemeris = ( translationIterator->second.second == true ) ? 1 : -1;

                translationalEphemerides_.push_back( std::make_pair( translationIterator->second.first,
                                                                     addCurrentEphemeris ) );
                isCurrentEphemerisTranslational_.push_back( 1 );
                translationIterator++;
            }
            // If current ephemeris is rotational, add to rotations list and set translation flag at
            // current index to false.
            else if( rotationIterator != rotationalEphemerides.end( )
                     && rotationIterator->first == currentIndex )
            {
                rotationalEphemerides_.push_back( rotationIterator->second );
                isCurrentEphemerisTranslational_.push_back( 0 );
                rotationIterator++;
            }
            // If index is not found, display error message.
            else
            {
                std::string errorMessage = "Error when  making composite ephemeris, input indices inconsistent";
                throw std::runtime_error( errorMessage );
            }
            currentIndex++;
        }
    }

    //! Destructor
    ~CompositeEphemeris( ){ }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Constant state given by combined rotations and translations.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch )
    {
        return getTemplatedCartesianStateFromCompositeEphemeris< double, double >( secondsSinceEpoch );
    }

    //! Get state from ephemeris (with long double as state scalar).
    /*!
     * Returns state from ephemeris with long double as state scalar at given time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Constant state with long double as state scalar given by combined rotations and translations.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongState(
            const double secondsSinceEpoch )
    {
        return getTemplatedCartesianStateFromCompositeEphemeris< double, long double >( secondsSinceEpoch );
    }

    //! Get state from ephemeris (with double as state scalar and Time as time type).
    /*!
     * Returns state from ephemeris with double as state scalar at given time (as custom Time type).
     * \param currentTime Time at which state is to be evaluated
     * \return State from ephemeris with long double as state scalar
     */
    Eigen::Matrix< double, 6, 1 > getCartesianStateFromExtendedTime(
            const Time& currentTime )
    {
        return getTemplatedCartesianStateFromCompositeEphemeris< Time, double >( currentTime );
    }

    //! Get state from ephemeris (with long double as state scalar and Time as time type).
    /*!
     * Returns state from ephemeris with long double as state scalar at given time (as custom Time type).
     * \param currentTime Time at which state is to be evaluated
     * \return State from ephemeris with long double as state scalar
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime(
            const Time& currentTime )
    {
        return getTemplatedCartesianStateFromCompositeEphemeris< Time, long double >( currentTime );
    }

    //! Templated function to get the state from tabulated ephemeris.
    /*!
     *  Templated function to get the state from tabulated ephemeris. This function is called with the appropriate
     *  template arguments by each of the specific state functions. Its function is to have only a single function
     *  implemented, while ensuring that the numerical precision is at least that of the CompositeEphemeris and requested
     *  time/state scalar types.
     *  \param currentTime Seconds since epoch at which ephemeris is to be evaluated.
     *  \return State given by combined rotations and translations, at requested precision.
     */
    template< typename OutputTimeType, typename OutputStateScalarType >
    Eigen::Matrix< OutputStateScalarType, 6, 1 > getTemplatedCartesianStateFromCompositeEphemeris(
            const OutputTimeType& currentTime )
    {

        // Initialize state to zero;
        Eigen::Matrix< StateScalarType, 6, 1 > state = Eigen::Matrix< StateScalarType, 6, 1 >::Zero( );

        // Initialize current indices of translation and rotation to zero.
        int currentTranslationIndex = 0;
        int currentRotationIndex = 0;

        // Loop over all ephemerides
        for( unsigned int i = 0; i < isCurrentEphemerisTranslational_.size( ); i++ )
        {
            // If current ephemeris is translational, add it and increment currentTranslationIndex
            if( isCurrentEphemerisTranslational_[ i ] == true )
            {
                state += translationalEphemerides_[ currentTranslationIndex ].first( static_cast< TimeType >( currentTime ) )*
                        static_cast< double >( translationalEphemerides_[ currentTranslationIndex ].second );
                currentTranslationIndex++;
            }
            // If current ephemeris is rotational, multiply position and state by rotation.
            else
            {
                state = rotationalEphemerides_[ currentRotationIndex ]( static_cast< TimeType >( currentTime ), state );
                currentRotationIndex++;
            }
        }

        return state.template cast< OutputStateScalarType >( );
    }

    //! Add an additional translational ephemeris at the end of the chain.
    /*!
     *  Function to add an additional translational ephemeris at the end of the chain.
     *  \param stateFunction Translational ephemeris function to add.
     *  \param add Identifier setting whether to add (1) or subtract (-1) translation.
     */
    void addTranslationalEphemeris( const std::function< StateType( const TimeType& ) > stateFunction,
                                    const int add = 1 )
    {
        //Check validity of input.
        if( add != 1 && add != -1 )
        {
            std::string errorMessage = "Error when adding to composite ephemeris";
            throw std::runtime_error( errorMessage );
        }

        // Add translational ephemeris to end of chain
        isCurrentEphemerisTranslational_.push_back( 1 );
        translationalEphemerides_.push_back( std::make_pair( stateFunction, add ) );
    }

private:

    //! Vector of translational ephemeris functions.
    /*!
     *  Vector of translational ephemeris functions and addition (1) or subtraction (-1) indicator.
     */
    std::vector< std::pair< std::function< StateType( const TimeType& ) >, int > > translationalEphemerides_;

    //! Vector of rotational ephemeris functions.
    /*!
     *  Vector of rotational ephemeris functions.
     */
    std::vector< std::function< StateType( const TimeType, const StateType& ) > > rotationalEphemerides_;

    //! Vector indicating order of translational and rotational ephemeris.
    /*!
     *  Vector indicating order of translational and rotational ephemeris (0 is first, highest is last)
     */
    std::vector< bool > isCurrentEphemerisTranslational_;
};

//! Interface function used to change the state scalar type of the output of a state function
/*!
 *  Interface function used to change the state scalar type of the output of a state function
 *  \param originalStateFunction State function with original scalar type
 *  \param currentTime Time at which state function is to be evaluated.
 *  \return State from originalStateFunction at currentTime, cast to NewStateScalarType.
 */
template< typename OldStateScalarType, typename NewStateScalarType, typename TimeType, int StateSize >
Eigen::Matrix< NewStateScalarType, StateSize, 1 > convertStateFunctionStateScalarOutput(
        const std::function< Eigen::Matrix< OldStateScalarType, StateSize, 1 >( const double& ) >
        originalStateFunction,
        const TimeType currentTime )
{
    return originalStateFunction( currentTime ).template cast< NewStateScalarType >( );
}

//! Function to create the state function of a reference point.
/*!
 *  Function to create the state function of a reference point, for instance the state of a ground station in a
 *  barycentric frame.
 *  \param bodyEphemeris Global ephemeris of body on which reference point is located
 *  \param bodyRotationModel Rotation model between global and body-fixed frame
 *  \param referencePointRelativeStateFunction Function returning the body-fixed state of the reference point as
 *  a function of time
 *  \return Ephemeris describbing the state of the local reference point in the global frame (expressed global frame's
 *  coordinates, i.e. same as bodyEphemeris).
 */
template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< Ephemeris > createReferencePointCompositeEphemeris(
        std::shared_ptr< Ephemeris > bodyEphemeris,
        std::shared_ptr< RotationalEphemeris > bodyRotationModel,
        std::function< Eigen::Vector6d( const double& ) > referencePointRelativeStateFunction )
{
    if( bodyEphemeris == nullptr )
    {
        throw std::runtime_error( "Error when creating reference point composite ephemeris, no body ephemeris is provided" );
    }

    if( bodyRotationModel == nullptr )
    {
        throw std::runtime_error( "Error when creating reference point composite ephemeris, no body rotation model is provided" );
    }

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    // Cast state fucntion of body (global) and reference point (local) into correct form.
    std::map< int, std::function< StateType( const TimeType& ) > > referencePointEphemerisVector;
    referencePointEphemerisVector[ 2 ] = std::bind(
                &Ephemeris::getTemplatedStateFromEphemeris< StateScalarType, TimeType >, bodyEphemeris, std::placeholders::_1 );
    referencePointEphemerisVector[ 0 ] = std::bind(
                &convertStateFunctionStateScalarOutput< double, StateScalarType, TimeType, 6 >,
                referencePointRelativeStateFunction, std::placeholders::_1 );


    // Crate rotation functions from local to global frame.
    std::function< Eigen::Quaterniond( const TimeType ) > rotationToFrameFunction =
            std::bind( &RotationalEphemeris::getRotationToBaseFrameTemplated< TimeType >, bodyRotationModel, std::placeholders::_1 );
    std::function< Eigen::Matrix3d( const TimeType ) > rotationMatrixToFrameDerivativeFunction =
            std::bind( &RotationalEphemeris::getDerivativeOfRotationToBaseFrameTemplated< TimeType >, bodyRotationModel, std::placeholders::_1 );

    // Create ephemeris
    std::map< int, std::function< StateType( const TimeType, const StateType& ) > > referencePointRotationVector;
    referencePointRotationVector[ 1 ] = std::bind(
                transformStateToFrameFromRotationTimeFunctions< StateScalarType, TimeType >,
                std::placeholders::_2, std::placeholders::_1, rotationToFrameFunction, rotationMatrixToFrameDerivativeFunction );

    return std::make_shared< CompositeEphemeris< TimeType, StateScalarType > >(
                referencePointEphemerisVector, referencePointRotationVector, "SSB", "ECLIPJ2000" );
}

//extern template class CompositeEphemeris< double, double >;

} // namespace ephemerides

} // namespace tudat
#endif // TUDAT_COMPOSITEEPHEMERIS_H
