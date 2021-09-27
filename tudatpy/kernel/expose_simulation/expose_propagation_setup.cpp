/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagation_setup.h"

#include "expose_propagation_setup/expose_dependent_variable_setup.h"
#include "expose_propagation_setup/expose_torque_setup.h"
#include "expose_propagation_setup/expose_acceleration_setup.h"
#include "expose_propagation_setup/expose_integrator_setup.h"
#include "expose_propagation_setup/expose_propagator_setup.h"
#include "expose_propagation_setup/expose_mass_rate_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy {
namespace simulation {
namespace propagation_setup {

void expose_propagation_setup(py::module &m) {


    auto acceleration_setup = m.def_submodule("acceleration");
    acceleration::expose_acceleration_setup(acceleration_setup);

    auto torque_setup = m.def_submodule("torque");
    torque::expose_torque_setup(torque_setup);

    auto integrator_setup = m.def_submodule("integrator");
    integrator::expose_integrator_setup(integrator_setup);

    auto propagator_setup = m.def_submodule("propagator");
    propagator::expose_propagator_setup(propagator_setup);

    auto mass_setup = m.def_submodule("mass");
    mass::expose_mass_rate_setup(mass_setup);

    auto dependent_variable_setup = m.def_submodule("dependent_variable");
    dependent_variable::expose_dependent_variable_setup(dependent_variable_setup);
}
}// namespace propagation_setup
}// namespace simulation
}// namespace tudatpy
