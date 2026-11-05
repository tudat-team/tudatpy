/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#include "tudat/io/serialization/file_io.h"

#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"

TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( tudat::simulation_setup::ObservationDependentVariableSettings )
TUDAT_IMPLEMENT_BINARY_IO( tudat::simulation_setup::ObservationDependentVariableBookkeeping )
