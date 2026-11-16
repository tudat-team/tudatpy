#include "expose_simulator_bindings.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "independent_arc_models.h"

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void expose_simulator_independent_arc_bindings( py::module& module )
{
    using namespace independent_arcs;
    py::class_< MarsGravityArcModel >( module, "MarsGravityArcModel", R"doc(
Native-only Mars/JGMRO120D gravity model for isolated Cowell arcs.

Load standard SPICE kernels before construction. No Python callbacks or custom
forces are accepted. Reference radius is in metres and mu in m^3/s^2.
)doc" )
            .def( py::init< int, int >( ), py::arg( "degree" ) = 100, py::arg( "order" ) = -1 )
            .def_readonly( "degree", &MarsGravityArcModel::degree )
            .def_readonly( "order", &MarsGravityArcModel::order )
            .def_readonly( "reference_radius", &MarsGravityArcModel::referenceRadius )
            .def_readonly( "gravitational_parameter", &MarsGravityArcModel::gravitationalParameter );

    // Intentionally retain the GIL. All workers are native-only; retaining it
    // prevents concurrent Python access to live native result containers.
    // Do not expose the underlying mutable solver or a Python factory here.
    py::class_< IndependentArcVariationalSimulator, std::shared_ptr< IndependentArcVariationalSimulator > >(
            module, "IndependentArcVariationalSimulator", R"doc(
Owned native multi-arc simulation with isolated Mars gravity models.

Initial states have shape (arcs, 6), in Mars/J2000 metres and m/s. Start epochs
are increasing TDB seconds since J2000; duration and step are seconds. Each arc
has six initial-state parameters and shares one Mars gravitational parameter.
Construction never integrates. One worker uses the sequential native route.
State-transfer configurations remain sequential even after NaNs are resolved.

Inputs and factory configuration are copied into C++; the model can be deleted
after construction. Native result types retain their usual live/update semantics.
The GIL remains held throughout integration. C++ workers require no GIL and
never execute Python callbacks. Arbitrary Python factories/custom models and
free-threaded Python are not supported by this first interface.
)doc" )
            .def( py::init< const MarsGravityArcModel&,
                            std::vector< double >,
                            const Eigen::MatrixXd&,
                            double,
                            double,
                            std::int64_t,
                            bool >( ),
                  py::arg( "model" ),
                  py::arg( "arc_start_times" ),
                  py::arg( "initial_states" ),
                  py::arg( "arc_duration" ),
                  py::arg( "step_size" ) = 30.0,
                  py::arg( "number_of_workers" ) = 1,
                  py::arg( "transfer_state_to_next_arc" ) = false )
            .def( "set_independent_arc_execution",
                  &IndependentArcVariationalSimulator::setIndependentArcExecution,
                  py::arg( "number_of_workers" ),
                  "Configure native arc workers; no Python factory is accepted." )
            .def( "clear_independent_arc_execution", &IndependentArcVariationalSimulator::clearIndependentArcExecution )
            .def( "integrate_full_equations",
                  &IndependentArcVariationalSimulator::integrateFullEquations,
                  py::arg( "initial_states" ) = py::none( ),
                  "Integrate dynamics and variations with optional copied concatenated initial states. Retains the GIL." )
            .def( "integrate_equations_of_motion_only",
                  &IndependentArcVariationalSimulator::integrateEquationsOfMotionOnly,
                  py::arg( "initial_states" ) = py::none( ),
                  "Integrate dynamics only; variational-result access is then invalid." )
            .def( "reset_parameter_estimate",
                  &IndependentArcVariationalSimulator::resetParameterEstimate,
                  py::arg( "parameter_vector" ),
                  py::arg( "reintegrate" ) = true,
                  "Update copied global parameters (6*arcs states, then shared Mars mu); optionally reintegrate." )
            .def_property_readonly( "number_of_workers", &IndependentArcVariationalSimulator::getNumberOfWorkers )
            .def_property_readonly( "last_integration_used_independent_arcs",
                                    &IndependentArcVariationalSimulator::lastIntegrationUsedIndependentArcs )
            .def_property_readonly( "initial_states", &IndependentArcVariationalSimulator::getInitialStates )
            .def_property_readonly( "parameter_vector", &IndependentArcVariationalSimulator::getParameterVector )
            .def_property_readonly( "propagation_results", &IndependentArcVariationalSimulator::getPropagationResults )
            .def_property_readonly( "variational_propagation_results", &IndependentArcVariationalSimulator::getVariationalResults )
            .def_property_readonly( "state_transition_interface", &IndependentArcVariationalSimulator::getStateTransitionInterface );
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
