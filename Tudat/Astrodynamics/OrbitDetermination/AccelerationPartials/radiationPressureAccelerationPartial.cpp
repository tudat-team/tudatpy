#include <iostream>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/radiationPressureAccelerationPartial.h"

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
        const Eigen::Vector3d& vectorToSource )
{
    return -radiationPressure * area / bodyMass * vectorToSource;
}


Eigen::Matrix3d CannonBallRadiationPressurePartial::wrtPositionOfAcceleratedBody( )
{
    return currentPartialWrtPosition_;
}

Eigen::Matrix3d CannonBallRadiationPressurePartial::wrtPositionOfAcceleratingBody( )
{
    return -currentPartialWrtPosition_;
}

std::pair< boost::function< Eigen::MatrixXd( ) >, int > CannonBallRadiationPressurePartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    boost::function< Eigen::MatrixXd( ) > partialFunction;
    int numberOfRows = 0;
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::radiation_pressure_coefficient:

            partialFunction = boost::bind( &CannonBallRadiationPressurePartial::wrtRadiationPressureCoefficient,
                                           this );
            numberOfRows = 1;

            break;
        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );
}


std::pair< boost::function< Eigen::MatrixXd( ) >, int > CannonBallRadiationPressurePartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    boost::function< Eigen::MatrixXd( ) > partialFunction;
    int numberOfRows = 0;
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {

        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );

}


}

}

}


