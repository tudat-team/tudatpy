/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/kernelExternTemplatesObservationSimulators.h"

namespace tudat
{

namespace observation_models
{

template std::vector< std::shared_ptr< ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >
createObservationSimulators< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies );

}  // namespace observation_models

}  // namespace tudat
