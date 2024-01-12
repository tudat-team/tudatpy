/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/acceleration_partials/radiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{


//! Calculates partial derivative of cannon ball radiation pressure acceleration wrt radiation pressure coefficient.
Eigen::Vector3d computePartialOfCannonBallRadiationPressureAccelerationWrtRadiationPressureCoefficient(
    const double radiationPressure,
    const double area,
    const double bodyMass,
    const Eigen::Vector3d& vectorToSource )
{
    return -radiationPressure * area / bodyMass * vectorToSource;
}

void CannonBallRadiationPressurePartial::wrtRadiationPressureCoefficient( Eigen::MatrixXd& partial )
{
    partial = computePartialOfCannonBallRadiationPressureAccelerationWrtRadiationPressureCoefficient(
        radiationPressureFunction_( ), areaFunction_( ), acceleratedBodyMassFunction_( ),
        ( sourceBodyState_( ) - acceleratedBodyState_( ) ).normalized( ) );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > CannonBallRadiationPressurePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Check if parameter dependency exists.
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        // Set function returning partial w.r.t. radiation pressure coefficient.
        case estimatable_parameters::radiation_pressure_coefficient:

            partialFunction = std::bind( &CannonBallRadiationPressurePartial::wrtRadiationPressureCoefficient,
                                           this, std::placeholders::_1 );
            numberOfRows = 1;

            break;
        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );
}


//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > CannonBallRadiationPressurePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Check if parameter dependency exists.
    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        // Set function returning partial w.r.t. radiation pressure coefficient.
        case estimatable_parameters::arc_wise_radiation_pressure_coefficient:

            if( std::dynamic_pointer_cast< estimatable_parameters::ArcWiseRadiationPressureCoefficient >( parameter ) != nullptr )
            {
                partialFunction = std::bind(
                            &CannonBallRadiationPressurePartial::wrtArcWiseRadiationPressureCoefficient, this, std::placeholders::_1,
                            std::dynamic_pointer_cast< estimatable_parameters::ArcWiseRadiationPressureCoefficient >( parameter ) );
            }
            else
            {
                throw std::runtime_error( "Error when making radiation pressure partial, arcwise radiation pressure parameter not consistent" );
            }
            numberOfRows = parameter->getParameterSize( );

            break;
        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );
}


}

}

