/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace reference_frames
{

Eigen::Vector3d computeBodyFixedAeroAngles(
        const Eigen::Matrix3d& inertialToBodyFixedFrame,
        const Eigen::Matrix3d& trajectoryToInertialFrame )
{
    // Retrieve rotation matrix that is to be converted to orientation angles.
    Eigen::Matrix3d currentRotationFromBodyToTrajectoryFrame_ =
            ( inertialToBodyFixedFrame * trajectoryToInertialFrame ).transpose( );

    // Compute associated Euler angles and set as orientation angles.
    Eigen::Vector3d eulerAngles = basic_mathematics::get132EulerAnglesFromRotationMatrix(
                currentRotationFromBodyToTrajectoryFrame_ );
    return ( Eigen::Vector3d( ) << eulerAngles( 0 ), eulerAngles( 1 ), -eulerAngles( 2 ) ).finished( );
}

}

namespace ephemerides
{


////! Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
//void setAerodynamicDependentOrientationCalculatorClosure(
//        const std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame,
//        std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator )
//{
//    using namespace reference_frames;
//    std::shared_ptr< AerodynamicAnglesClosure > aerodynamicAnglesClosure =
//            std::make_shared< AerodynamicAnglesClosure >(
//                imposedRotationFromInertialToBodyFixedFrame, aerodynamicAngleCalculator );
////    aerodynamicAngleCalculator->setOrientationAngleFunctions(
////                std::bind( &AerodynamicAnglesClosure::getCurrentAngleOfAttack, aerodynamicAnglesClosure ),
////                std::bind( &AerodynamicAnglesClosure::getCurrentAngleOfSideslip, aerodynamicAnglesClosure ),
////                std::bind( &AerodynamicAnglesClosure::getCurrentBankAngle, aerodynamicAnglesClosure ),
////                std::bind( &AerodynamicAnglesClosure::updateAngles, aerodynamicAnglesClosure, std::placeholders::_1 ) );
//}

////! Function to make aerodynamic angle computation consistent with existing DependentOrientationCalculator
//void setAerodynamicDependentOrientationCalculatorClosure(
//        std::shared_ptr< DependentOrientationCalculator > dependentOrientationCalculator,
//        std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
//{
//    std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame =
//            std::bind( &DependentOrientationCalculator::computeAndGetRotationToLocalFrame, dependentOrientationCalculator, std::placeholders::_1 );
//    setAerodynamicDependentOrientationCalculatorClosure(
//                imposedRotationFromInertialToBodyFixedFrame,
//                aerodynamicAngleCalculator );
//}

////! Function to make aerodynamic angle computation consistent with existing rotational ephemeris
//void setAerodynamicDependentOrientationCalculatorClosure(
//        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
//        std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator )
//{
//    using namespace reference_frames;
//    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > angleBasedRotationModel =
//            std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >( rotationalEphemeris );

//    if( angleBasedRotationModel == nullptr )
//    {
//        setAerodynamicDependentOrientationCalculatorClosure(
//                    std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
//                               rotationalEphemeris, std::placeholders::_1 ),
//                    aerodynamicAngleCalculator );
//    }
//    else if( aerodynamicAngleCalculator != angleBasedRotationModel->getAerodynamicAngleCalculator( ) )
//    {
//        std::cout<<aerodynamicAngleCalculator<<" "<<angleBasedRotationModel->getAerodynamicAngleCalculator( )<<std::endl;
//        std::cout<<aerodynamicAngleCalculator.get( )<<" "<<angleBasedRotationModel->getAerodynamicAngleCalculator( ).get( )<<std::endl;

//        throw std::runtime_error( "Error, body has AerodynamicAngleRotationalEphemeris and FlightConditions, but angle calculators are not compatible" );
//    }
//}

void verifyAerodynamicDependentOrientationCalculatorClosure(
        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{

}

} // namespace tudat
} // namespace ephemerides
