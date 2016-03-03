#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

bool isBodyAccelerationDependentOnBody( const basic_astrodynamics::SingleBodyAccelerationMap& bodyAccelerations, const std::string bodyToCheck )
{
    bool isAccelerationDependentOnRequestedBody = 0;
    for( basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationIterator = bodyAccelerations.begin( );
         accelerationIterator != bodyAccelerations.end( ); accelerationIterator++ )
    {
        if( accelerationIterator->first == bodyToCheck )
        {
            isAccelerationDependentOnRequestedBody = 1;
        }
    }
    return isAccelerationDependentOnRequestedBody;
}

bool isBodyAccelerationDependentOnBody( const basic_astrodynamics::SingleBodyAccelerationMap& bodyAccelerations,
                                        const std::vector< std::string > bodiesToCheck )
{
    bool isAccelerationDependentOnRequestedBody = 0;
    for( unsigned int i = 0; i < bodiesToCheck.size( ); i++ )
    {
        if( isBodyAccelerationDependentOnBody( bodyAccelerations, bodiesToCheck.at( i ) ) )
        {
            isAccelerationDependentOnRequestedBody = 1;
        }
    }
    return isAccelerationDependentOnRequestedBody;
}

bool isBodyListAccelerationDependentOnBody( const basic_astrodynamics::AccelerationMap& bodyAccelerations,
                                            const std::vector< std::string > bodiesToCheck )
{
    bool isAccelerationDependentOnRequestedBody = 0;
    for( basic_astrodynamics::AccelerationMap::const_iterator accelerationIterator = bodyAccelerations.begin( );
         accelerationIterator != bodyAccelerations.end( ); accelerationIterator++ )
    {
        if( isBodyAccelerationDependentOnBody( accelerationIterator->second, bodiesToCheck ) )
        {
            isAccelerationDependentOnRequestedBody = 1;
        }
    }
    return isAccelerationDependentOnRequestedBody;
}

}

}
