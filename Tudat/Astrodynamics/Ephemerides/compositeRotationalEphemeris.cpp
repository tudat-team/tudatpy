#include <iostream>

#include "Astrodynamics/Ephemerides/compositeRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

CompositeRotationalEphemeris::CompositeRotationalEphemeris( std::vector< boost::shared_ptr< RotationalEphemeris > > subRotations ):
    subRotations_( subRotations )
{
    for( unsigned int i = 0; i < subRotations_.size( ) - 1; i++ )
    {
        if( subRotations_[ i ]->getTargetFrameOrientation( ) != subRotations_[ i + 1 ]->getBaseFrameOrientation( ) )
        {
            std::cerr<<"Error when making composite rotational ephemeris, target and original frames of subsequent rotations"<<
                       " are incompatible, namely "<<subRotations_[ i ]->getTargetFrameOrientation( )<<" and "<<
                       subRotations_[ i ]->getBaseFrameOrientation( )<<std::endl;
        }
    }

    baseFrameOrientation_ = subRotations_[ 0 ]->getBaseFrameOrientation( );
    targetFrameOrientation_ = subRotations_[ subRotations_.size( ) - 1 ]->getTargetFrameOrientation( );
}

Eigen::Quaterniond CompositeRotationalEphemeris::getRotationToBaseFrame( const double ephemerisTime )
{
    Eigen::Quaterniond rotation = subRotations_[ 0 ]->getRotationToBaseFrame( ephemerisTime );

    for( unsigned int i = 1; i < subRotations_.size( ); i++ )
    {
        rotation *= subRotations_[ i ]->getRotationToBaseFrame( ephemerisTime );
    }

    return rotation;
}

Eigen::Matrix3d CompositeRotationalEphemeris::getDerivativeOfRotationFromFrame( const double ephemerisTime )
{
    Eigen::Matrix3d totalRotationDerivative = Eigen::Matrix3d::Zero( );

    std::vector< Eigen::Matrix3d > subRotationMatrices;
    for( int i = 0; i < static_cast< int >( subRotations_.size( ) ); i++ )
    {
        subRotationMatrices.push_back( ( subRotations_[ i ]->getRotationToBaseFrame( ephemerisTime ) ).toRotationMatrix( ) );
    }

    Eigen::Matrix3d singleLevelAddition;
    for( int i = 0; i < static_cast< int >( subRotations_.size( ) ); i++ )
    {
        singleLevelAddition.setIdentity( );

        for( int j = 0; j <= i - 1; j++ )
        {
            singleLevelAddition *= subRotationMatrices.at( j );
        }

        singleLevelAddition *= subRotations_[ i ]->getDerivativeOfRotationFromFrame( ephemerisTime );

        for( int j = i + 1; j < static_cast< int >( subRotations_.size( ) ); j++ )
        {
            singleLevelAddition *=  subRotationMatrices.at( j );
        }

        totalRotationDerivative += singleLevelAddition;
    }
    return totalRotationDerivative;
}

}

}
