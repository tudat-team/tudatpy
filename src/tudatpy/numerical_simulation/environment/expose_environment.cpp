/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_environment.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/aerodynamics.h>
#include <tudat/astro/ephemerides.h>
#include <tudat/astro/gravitation.h>
#include <tudat/basics/deprecationWarnings.h>

#include "scalarTypes.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

// namespace py = pybind11;
// namespace tba = tudat::basic_astrodynamics;
// namespace tss = tudat::simulation_setup;
// namespace tp = tudat::propagators;
// namespace tinterp = tudat::interpolators;
// namespace te = tudat::ephemerides;
// namespace tni = tudat::numerical_integrators;
// namespace trf = tudat::reference_frames;
// namespace tmrf = tudat::root_finders;

namespace py = pybind11;

namespace tba = tudat::basic_astrodynamics;
namespace ta = tudat::aerodynamics;
namespace tr = tudat::reference_frames;
namespace te = tudat::ephemerides;
namespace teo = tudat::earth_orientation;
namespace tgs = tudat::ground_stations;
namespace tr = tudat::reference_frames;
namespace tg = tudat::gravitation;
namespace trf = tudat::reference_frames;
namespace tss = tudat::simulation_setup;
namespace ti = tudat::interpolators;
namespace tsm = tudat::system_models;
namespace tom = tudat::observation_models;
namespace tem = tudat::electromagnetism;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment
{

void expose_environment( py::module &m )
{

}
}  // namespace environment
}  // namespace numerical_simulation
}  // namespace tudatpy
