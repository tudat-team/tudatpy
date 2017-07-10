/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H
#define TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H

#include <iomanip>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{

namespace gravitation
{

Eigen::Vector3d calculateSecondDegreeGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodySubjectToTorque,
        const double gravitationalParameterOfAttractingBody,
        const Eigen::Matrix3d& inertiaTensorOfRotatingBody );

class SecondDegreeGravitationalTorqueModel: public basic_astrodynamics::TorqueModel
{
public:
    SecondDegreeGravitationalTorqueModel(
            const boost::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction,
            const boost::function< double( ) > gravitationalParameterOfAttractingBodyFunction,
            const boost::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction,
            const boost::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction =
            boost::lambda::constant( Eigen::Vector3d::Zero( ) ),
            const boost::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ) ):
        positionOfBodySubjectToTorqueFunction_( positionOfBodySubjectToTorqueFunction ),
        gravitationalParameterOfAttractingBodyFunction_( gravitationalParameterOfAttractingBodyFunction ),
        inertiaTensorOfRotatingBodyFunction_( inertiaTensorOfRotatingBodyFunction ),
        positionOfBodyExertingTorqueFunction_( positionOfBodyExertingTorqueFunction ),
        rotationToBodyFixedFrameFunction_( rotationToBodyFixedFrameFunction )
    {

    }

    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    void updateMembers( const double currentTime )
    {
        currentRelativePositionOfBodySubjectToTorque_ = positionOfBodyExertingTorqueFunction_( ) - positionOfBodySubjectToTorqueFunction_( );
        currentRotationToBodyFixedFrameFunction_ = rotationToBodyFixedFrameFunction_( );

        currentGravitationalParameterOfAttractingBody_ = gravitationalParameterOfAttractingBodyFunction_( );
        currentInertiaTensorOfRotatingBody_ = inertiaTensorOfRotatingBodyFunction_( );

        currentTorque_ = calculateSecondDegreeGravitationalTorque(
                    currentRotationToBodyFixedFrameFunction_ * currentRelativePositionOfBodySubjectToTorque_,
                    currentGravitationalParameterOfAttractingBody_,
                    currentInertiaTensorOfRotatingBody_ );
    }

protected:

private:
    boost::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction_;

    boost::function< double( ) > gravitationalParameterOfAttractingBodyFunction_;

    boost::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction_;

    boost::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction_;

    const boost::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction_;


    Eigen::Vector3d currentRelativePositionOfBodySubjectToTorque_;

    double currentGravitationalParameterOfAttractingBody_;

    Eigen::Matrix3d currentInertiaTensorOfRotatingBody_;

    Eigen::Quaterniond currentRotationToBodyFixedFrameFunction_;

    Eigen::Vector3d currentTorque_;
};

}

}

#endif // TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H
