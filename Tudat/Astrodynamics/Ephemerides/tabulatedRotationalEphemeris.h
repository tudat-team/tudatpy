#ifndef TABULATEDROTATIONALEPHEMERIS_H
#define TABULATEDROTATIONALEPHEMERIS_H


#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Mathematics/Interpolators/slerpInterpolator.h"

namespace tudat
{

namespace ephemerides
{

////! Class for calculating rotational state of a body from a time-tabulated set of quaternions.
///*!
// *  Class for calculating rotational state of a body from a time-tabulated set of quaternions. The class performs interpolation of the
// *  tabulated quaternions by using spherical linear interpolation (SLERP), which interpolates over the unit sphere linearly along a great
// *  circle arc, ensuring the norm condition of the rotation quaternion.
// */
//class SlerpTabulatedRotationalEphemeris : public RotationalEphemeris
//{
//public:

//    //! Constructor taking quaternion map
//    /*!
//     *  Constructor taking quaternion map, with time values as key. The map is split into two vectors by the constructor to facilitate the
//     *  interpolation.
//     *  \param quaternionMap Map of quaternions, denoting the rotation from the local(target) to the inertial(original) frame, with keys
//     *  denoting the times of teh rotations.
//     *  \param originalFrameOrientation Base frame identifier.
//     *  \param targetFrameOrientation Target frame identifier.
//     */
//    SlerpTabulatedRotationalEphemeris( std::map< double, Eigen::Quaterniond >& quaternionMap,
//                                       const std::string& originalFrameOrientation = "ECLIPJ2000",
//                                       const std::string& targetFrameOrientation = "" );

//    SlerpTabulatedRotationalEphemeris( std::map< double, Eigen::Quaterniond >& quaternionMap,
//                                       std::map< double, Eigen::Matrix3d >& rotationDerivativeMap,
//                                       const std::string& originalFrameOrientation = "ECLIPJ2000",
//                                       const std::string& targetFrameOrientation = "" );

//    ~SlerpTabulatedRotationalEphemeris( ){ }

//    //! Function to reset the tabulated quaternion values.
//    /*!
//     *  Function to reset the tabulated quaternion value map, with time values as key. The map is split into two vectors by the constructor
//     *  to facilitate the interpolation.
//     *  \param quaternionMap Map of quaternions, denoting the rotation from the local(target) to the inertial(original) frame, with keys
//     *  denoting the times of teh rotations.
//     */
//    void reset( const std::map< double, Eigen::Quaterniond >& quaternionMap );

//    //! Function to calculate the rotation quaternion from target frame to original frame.
//    /*!
//     *  Function to calculate the rotation quaternion from target frame to original frame at specified time. Note that the input quaternions to
//     *  the constructor or the reset function tabulate the rotation calculated by this function, not the inverse getRotationToTargetFrame function.
//     *  \param ephemerisTime Time at which rotation is to be calculated.
//     *  \return Rotation from target (typically local) to original (typically global) frame at specified time.
//     */
//    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime );

//    Eigen::Quaterniond getRotationToTargetFrame(
//            const double secondsSinceEpoch )
//    {
//        return getRotationToBaseFrame( secondsSinceEpoch ).inverse( );
//    }

//    //! Function to calculate the time derivative of the rotation matrix from target frame to original frame at specified time.
//    /*!
//     *  Function to calculate the time derivative of the rotation matrix from target frame to original frame at specified time.
//     */
//    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime )
//    {
//        if( rotationMatrixDerivativeInterpolator_ == NULL )
//        {
//            std::cerr<<"Error, derivative of tabulated rotation not yet implemented"<<std::endl;
//            return Eigen::Matrix3d::Zero( );
//        }
//        else
//        {
//            return rotationMatrixDerivativeInterpolator_->interpolate( ephemerisTime );
//        }
//    }

//    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double ephemerisTime )
//    {
//        return getDerivativeOfRotationToBaseFrame( ephemerisTime ).transpose( );
//    }

//    //! Get fraction in tabulated time interval for last call to class
//    /*!
//     *  Get fraction in tabulated time interval for last call to class
//     *  \return Fraction in tabulated time interval for last call to class
//     */
//    double getDifference( )
//    {
//        return quaternionInterpolator_->getPreviousDifference( );
//    }

//private:

//    boost::shared_ptr< interpolators::SlerpInterpolator< double > > quaternionInterpolator_;

//    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix3d > > rotationMatrixDerivativeInterpolator_;

//};

template< typename StateScalarType, typename TimeType >
class TabulatedRotationalEphemeris : public RotationalEphemeris
{
public:

    typedef Eigen::Matrix< StateScalarType, 7, 1 > StateType;
    typedef Eigen::Matrix< StateScalarType, 4, 1 > OrientationType;

    //! Constructor, sets data interpolator and frame data.
    /*!
     *  \param interpolator Interpolator that returns the interpolated state as a function of time.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    TabulatedRotationalEphemeris( const boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > >
                                  interpolatorToBaseFrame,
                                  const std::string& originalFrameOrientation = "ECLIPJ2000",
                                  const std::string& targetFrameOrientation = "" ):
        RotationalEphemeris( originalFrameOrientation, targetFrameOrientation ), interpolator_( interpolatorToBaseFrame )
    {  }

    ~TabulatedRotationalEphemeris( ){ }

    //! Function to reset the state interpolator.
    /*!
     *  Function to reset the state interpolator, for instance following an update of the states of the body
     *  after a new numerical integration.
     *  \param interpolator New interpolator that returns the interpolated state as a function of time.
     */
    void reset( const boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > > interpolator )
    {
        interpolator_ = interpolator;
    }


    boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > > getInterpolator( )
    {
        return interpolator_;
    }

    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime )
    {
        updateInterpolator( ephemerisTime );

        return currentRotationToBaseFrame_.template cast< double >( );
    }

    Eigen::Quaterniond getRotationToTargetFrame( const double ephemerisTime )
    {
        return ( getRotationToBaseFrame( ephemerisTime ) ).inverse( );
    }

    Eigen::Vector3d getRotationalVelocityVectorInBaseFrame( const double ephemerisTime )
    {
        updateInterpolator( ephemerisTime );

        return ( currentRotationToBaseFrame_ * currentRotationalVelocityVectorInTargetFrame_ ).template cast< double >( );
    }

    Eigen::Vector3d getRotationalVelocityVectorInTargetFrame( const double ephemerisTime )
    {
        updateInterpolator( ephemerisTime );

        return currentRotationalVelocityVectorInTargetFrame_.template cast< double >( );
    }

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double ephemerisTime )
    {
        updateInterpolator( ephemerisTime );

        return getDerivativeOfRotationMatrixToFrame(
                    ( currentRotationToBaseFrame_.inverse( ) ).toRotationMatrix( ).template cast< double >( ),
                     ( currentRotationToBaseFrame_ * currentRotationalVelocityVectorInTargetFrame_ ).template cast< double >( ) );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime )
    {
        updateInterpolator( ephemerisTime );

        return getDerivativeOfRotationToTargetFrame( ephemerisTime ).transpose( );
    }


private:

    void updateInterpolator( const double time )
    {
        if( ! ( time == currentTime_ ) )
        {
            currentRotationalState_ = interpolator_->interpolate( time );

            double quaternionNorm = ( currentRotationalState_.block( 0, 0, 4, 1 ) ).norm( );
            currentRotationalState_.block( 0, 0, 4, 1 ) = currentRotationalState_.block( 0, 0, 4, 1 ) / quaternionNorm;
            currentRotationToBaseFrame_ = Eigen::Quaternion< StateScalarType >(
                        currentRotationalState_( 0 ), currentRotationalState_( 1 ),
                        currentRotationalState_( 2 ), currentRotationalState_( 3 ) );
            currentRotationalVelocityVectorInTargetFrame_ = currentRotationalState_.block( 4, 0, 3, 1 );
            currentTime_ = time;
        }
    }

    //! Interpolator that returns body state as a function of time.
    /*!
     *  Interpolator that returns body state as a function of time by calling the interpolate function (i.e. time as
     *  independent variable and states as dependent variables ).
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > > interpolator_;

    double currentTime_;

    StateType currentRotationalState_;

    Eigen::Matrix< StateScalarType, 3, 1 > currentRotationalVelocityVectorInTargetFrame_;

    Eigen::Quaternion< StateScalarType > currentRotationToBaseFrame_;

};

}

}

#endif // TABULATEDROTATIONALEPHEMERIS_H

