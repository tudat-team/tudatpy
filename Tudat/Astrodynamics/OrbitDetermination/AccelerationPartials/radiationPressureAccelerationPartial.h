/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSUREACCELERATIONPARTIAL_H
#define TUDAT_RADIATIONPRESSUREACCELERATIONPARTIAL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

Eigen::Vector3d computePartialOfCannonBallRadiationPressureAccelerationWrtRadiationPressureCoefficient(
        const double radiationPressure,
        const double area,
        const double bodyMass,
        const Eigen::Vector3d& vectorToSource );

class CannonBallRadiationPressurePartial: public AccelerationPartial
{
public:
    CannonBallRadiationPressurePartial(
            const boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface,
            const boost::function< double( ) > massFunction, const std::string& acceleratedBody, const std::string& acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody,
                             basic_astrodynamics::cannon_ball_radiation_pressure ),
        sourceBodyState_( radiationPressureInterface->getSourcePositionFunction( ) ),
        acceleratedBodyState_( radiationPressureInterface->getTargetPositionFunction( ) ),
        areaFunction_( boost::bind( &electro_magnetism::RadiationPressureInterface::getArea, radiationPressureInterface ) ),
        radiationPressureCoefficientFunction_( boost::bind( &electro_magnetism::RadiationPressureInterface::getRadiationPressureCoefficient, radiationPressureInterface ) ),
        radiationPressureFunction_( boost::bind( &electro_magnetism::RadiationPressureInterface::getCurrentRadiationPressure,
                                                 radiationPressureInterface ) ),
        acceleratedBodyMassFunction_( massFunction )
    { }

    ~CannonBallRadiationPressurePartial( ){ }

    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }


    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }

    std::pair< boost::function< Eigen::MatrixXd( ) >, int >
    getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );


    std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    Eigen::MatrixXd wrtRadiationPressureCoefficient( )
    {
        return computePartialOfCannonBallRadiationPressureAccelerationWrtRadiationPressureCoefficient(
                    radiationPressureFunction_( ), areaFunction_( ), acceleratedBodyMassFunction_( ),
                    ( sourceBodyState_( ) - acceleratedBodyState_( ) ).normalized( ) );
    }


    void update( const double currentTime = 0.0 )
    {
        if( !( currentTime_ == currentTime ) )
        {
            Eigen::Vector3d rangeVector = ( acceleratedBodyState_( ) - sourceBodyState_( ) );
            double range = rangeVector.norm( );
            double rangeInverse = 1.0 / ( range );
            currentPartialWrtPosition_ =  ( radiationPressureCoefficientFunction_( ) * areaFunction_( ) * radiationPressureFunction_( ) /
                                            acceleratedBodyMassFunction_( ) ) * ( Eigen::Matrix3d::Identity( ) * rangeInverse - 3.0 *
                                                                                  rangeVector * rangeVector.transpose( ) * rangeInverse / (
                                                                                      range * range ) );
            currentTime_ = currentTime;
        }
    }

private:

    boost::function< Eigen::Vector3d( ) > sourceBodyState_;

    boost::function< Eigen::Vector3d( )> acceleratedBodyState_;

    boost::function< double( ) > areaFunction_;

    boost::function< double( ) > radiationPressureCoefficientFunction_;

    boost::function< double( ) > radiationPressureFunction_;

    boost::function< double( ) > acceleratedBodyMassFunction_;

    Eigen::Matrix3d currentPartialWrtPosition_;
};

}

}

}

#endif // TUDAT_RADIATIONPRESSUREACCELERATIONPARTIAL_H
