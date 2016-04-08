/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef NUMERICALACCELERATIONPARTIAL_H
#define NUMERICALACCELERATIONPARTIAL_H

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

void emptyFunction( );

void emptyTimeFunction( const double time );

Eigen::Matrix3d calculateAccelerationWrtStatePartials(
        boost::function< void( basic_mathematics::Vector6d ) > setBodyState,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        basic_mathematics::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        boost::function< void( ) > updateFunction = emptyFunction );


Eigen::Vector3d calculateAccelerationWrtParameterPartials(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double parameterPerturbation,
        boost::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        boost::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );


Eigen::Matrix< double, 3, Eigen::Dynamic > calculateAccelerationWrtParameterPartials(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::VectorXd parameterPerturbation,
        boost::function< void( ) > updateDependentVariables = emptyFunction,
        const double currentTime = 0.0,
        boost::function< void( const double ) > timeDependentUpdateDependentVariables = emptyTimeFunction );


}

}

}

#endif // NUMERICALACCELERATIONPARTIAL_H
