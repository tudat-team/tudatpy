#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Ephemerides/tabulatedRotationalEphemeris.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace ephemerides
{

////! Constructor taking quaternion map
//SlerpTabulatedRotationalEphemeris::SlerpTabulatedRotationalEphemeris( std::map< double, Eigen::Quaterniond >& quaternionMap,
//                                                            const std::string& originalFrameOrientation,
//                                                            const std::string& targetFrameOrientation ):
//    RotationalEphemeris( originalFrameOrientation, targetFrameOrientation )
//{
//    // Set quaternion values and times.
//    reset( quaternionMap );
//}

//SlerpTabulatedRotationalEphemeris::SlerpTabulatedRotationalEphemeris( std::map< double, Eigen::Quaterniond >& quaternionMap,
//                                                            std::map< double, Eigen::Matrix3d >& rotationDerivativeMap,
//                                                            const std::string& originalFrameOrientation,
//                                                            const std::string& targetFrameOrientation ):
//    RotationalEphemeris( originalFrameOrientation, targetFrameOrientation )
//{
//    reset( quaternionMap );
//    rotationMatrixDerivativeInterpolator_ = boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix3d > >(
//                rotationDerivativeMap );

//}

////! Function to reset the tabulated quaternion values.
//void SlerpTabulatedRotationalEphemeris::reset( const std::map< double, Eigen::Quaterniond >& quaternionMap )
//{
//    quaternionInterpolator_ = boost::make_shared< interpolators::SlerpInterpolator< double > >(
//                quaternionMap );
//}

////! Function to calculate the rotation quaternion from target frame to original frame.
//Eigen::Quaterniond SlerpTabulatedRotationalEphemeris::getRotationToBaseFrame( const double ephemerisTime )
//{
//    return quaternionInterpolator_->interpolate( ephemerisTime );
//}

}

}
