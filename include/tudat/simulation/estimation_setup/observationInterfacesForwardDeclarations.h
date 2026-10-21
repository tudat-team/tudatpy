/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONINTERFACESFORWARDDECLARATIONS_H
#define TUDAT_OBSERVATIONINTERFACESFORWARDDECLARATIONS_H

#include <type_traits>

#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace observation_models
{

class ObservationModelSettings;
class ObservationBiasSettings;
class DopplerProperTimeRateSettings;
class LightTimeCorrectionSettings;
struct ObservationAncillarySimulationSettings;
class ObservationViabilitySettings;

template< typename ObservationScalarType, typename TimeType >
class ObservationManagerBase;

template< typename ObservationScalarType, typename TimeType >
class ObservationSimulatorBase;

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
class ObservationCollection;

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
class SingleObservationSet;

struct ObservationCollectionParser;

}  // namespace observation_models

namespace simulation_setup
{

class ObservationDependentVariableSettings;
class ObservationDependentVariableBookkeeping;
class ObservationDependentVariableCalculator;

template< typename TimeType >
class ObservationSimulationSettings;

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONINTERFACESFORWARDDECLARATIONS_H
