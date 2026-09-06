/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#include "tudat/io/serialization/file_io.h"

#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationPrintSettings.h"
#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"

TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::simulation_setup::AccelerationSettings )
TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::simulation_setup::TorqueSettings )
TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::propagators::VariableSettings )
TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::propagators::PropagationTerminationSettings )
TUDAT_IMPLEMENT_FILE_IO( tudat::propagators::PropagationPrintSettings )
TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::propagators::PropagationTerminationDetails )
TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::propagators::PropagatorProcessingSettings )
