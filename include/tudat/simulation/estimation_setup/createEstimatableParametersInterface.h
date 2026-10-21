/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEESTIMATABLEPARAMETERSINTERFACE_H
#define TUDAT_CREATEESTIMATABLEPARAMETERSINTERFACE_H

#include <memory>
#include <vector>

namespace tudat
{

namespace simulation_setup
{
class SystemOfBodies;
}

namespace estimatable_parameters
{
class EstimatableParameterSettings;

template< typename ParameterType >
class EstimatableParameter;

template< typename ParameterType >
class EstimatableParameterSet;

}  // namespace estimatable_parameters

namespace propagators
{

template< typename StateScalarType >
class PropagatorSettings;

template< typename StateScalarType, typename TimeType >
class SingleArcPropagatorSettings;

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_CREATEESTIMATABLEPARAMETERSINTERFACE_H
