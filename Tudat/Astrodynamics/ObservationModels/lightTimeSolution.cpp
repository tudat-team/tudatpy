/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{
namespace observation_models
{


template< >
double getDefaultLightTimeTolerance< double, double >( )
{
    return 1.0E-12;
}

template< >
long double getDefaultLightTimeTolerance< long double, long double >( )
{
    return 1.0E-14L;
}

template< >
double getDefaultLightTimeTolerance< double, long double >( )
{
    return 1.0E-12;
}

template< >
long double getDefaultLightTimeTolerance< long double, double >( )
{
    return 1.0E-12;
}

} // namespace observation_models
} // namespace tudat
