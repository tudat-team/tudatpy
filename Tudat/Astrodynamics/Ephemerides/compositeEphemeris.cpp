#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Get state from ephemeris.
template< >
basic_mathematics::Vector6d CompositeEphemeris< double, double >::getCartesianStateFromEphemeris(
        const double ephemerisTime, const double julianDayAtEpoch )
{
    // Initialize state to zero;
    basic_mathematics::Vector6d state = basic_mathematics::Vector6d::Zero( );

    // Initialize current indices of translation and rotation to zero.
    int currentTranslationIndex = 0;
    int currentRotationIndex = 0;

    // Loop over all ephemerides
    for( unsigned int i = 0; i < isCurrentEphemerisTranslational_.size( ); i++ )
    {
        // If current ephemeris is translational, add it and increment currentTranslationIndex
        if( isCurrentEphemerisTranslational_[ i ] == true )
        {
            state += translationalEphemerides_[ currentTranslationIndex ].first( ephemerisTime ) *
                    static_cast< double >( translationalEphemerides_[ currentTranslationIndex ].second );
            currentTranslationIndex++;
        }
        // If current ephemeris is rotational, multiply position and state by rotation.
        else
        {
            state = rotationalEphemerides_[ currentRotationIndex ]( ephemerisTime, state );
            currentRotationIndex++;
        }
    }

    return state;
}


template< >
Eigen::Matrix< long double, 6, 1 > CompositeEphemeris< double, double >::getCartesianLongStateFromEphemeris(
        const double ephemerisTime, const double julianDayAtEpoch )
{
    return getCartesianStateFromEphemeris( ephemerisTime, julianDayAtEpoch ).cast< long double >( );
}


template< >
Eigen::Matrix< long double, 6, 1 > CompositeEphemeris< double, long double >::getCartesianLongStateFromEphemeris(
        const double ephemerisTime, const double julianDayAtEpoch )

{
    // Initialize state to zero;
    Eigen::Matrix< long double, 6, 1 > state = Eigen::Matrix< long double, 6, 1 >::Zero( );

    // Initialize current indices of translation and rotation to zero.
    int currentTranslationIndex = 0;
    int currentRotationIndex = 0;

    // Loop over all ephemerides
    for( unsigned int i = 0; i < isCurrentEphemerisTranslational_.size( ); i++ )
    {
        // If current ephemeris is translational, add it and increment currentTranslationIndex
        if( isCurrentEphemerisTranslational_[ i ] == true )
        {
            state += translationalEphemerides_[ currentTranslationIndex ].first( ephemerisTime ) *
                    static_cast< long double >( translationalEphemerides_[ currentTranslationIndex ].second );
            currentTranslationIndex++;
        }
        // If current ephemeris is rotational, multiply position and state by rotation.
        else
        {
            state = rotationalEphemerides_[ currentRotationIndex ]( ephemerisTime, state );
            currentRotationIndex++;
        }
    }
    return state;
}


template< >
basic_mathematics::Vector6d CompositeEphemeris< double, long double >::getCartesianStateFromEphemeris(
        const double ephemerisTime, const double julianDayAtEpoch )

{
    return getCartesianLongStateFromEphemeris( ephemerisTime, julianDayAtEpoch ).cast< double >( );
}

}

}
