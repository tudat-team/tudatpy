//
// Created by elmar on 13-07-20.
//

#include "expose_low_thrust.h"

namespace py = pybind11;
namespace tltt = tudat::low_thrust_trajectories;
namespace sbm = tudat::shape_based_methods;

namespace tudatpy {

	void expose_low_thrust( py::module &m )
	{

		py::enum_<tltt::LowThrustLegTypes>(m, "LowThrustLegTypes")
				.value("hodographic_shaping_leg", tltt::hodographic_shaping_leg)
				.value("spherical_shaping_leg", tltt::spherical_shaping_leg)
				.export_values();

		m.def("create_low_thrust_leg", &tltt::createLowThrustLeg, py::arg("low_thrust_leg_settings"),
			  py::arg("state_at_departure"), py::arg("state_at_arrival"), py::arg("time_of_flight"));

		py::class_<tltt::SphericalShapingLegSettings, std::shared_ptr<tltt::SphericalShapingLegSettings>>(m, "SphericalShapingLegSettings")
		.def(py::init<
		        const int,
		        const double,
		        const double,
				const std::shared_ptr< tudat::root_finders::RootFinderSettings >,
				const std::pair< double, double >
		        >(),
		        py::arg("number_of_revolutions"),
		        py::arg("central_body_gravitational_parameter"),
		        py::arg("initial_value_free_coefficient"),
		        py::arg("root_finder_settings"),
		        py::arg("bounds_free_coefficient") = std::make_pair( TUDAT_NAN, TUDAT_NAN )
			);

		py::class_<sbm::SphericalShaping, std::shared_ptr<sbm::SphericalShaping>>(m, "SphericalShaping")
		.def(py::init<
		        const Eigen::Vector6d&,
		        const Eigen::Vector6d&,
		        const double,
		        const double,
		        const int,
		        const double,
		        const std::shared_ptr<tudat::root_finders::RootFinderSettings>,
				const double,
				const double,
				const double
				>(),
				py::arg("initial_state"),
				py::arg("final_state"),
				py::arg("required_time_of_flight"),
				py::arg("central_body_gravitational_parameter"),
				py::arg("number_of_revolutions"),
				py::arg("initial_value_free_coefficient"),
				py::arg("root_finder_settings"),
				py::arg("lower_bound_free_coefficient") = TUDAT_NAN,
			 	py::arg("upper_bound_free_coefficient") = TUDAT_NAN,
			 	py::arg("initial_body_mass") = TUDAT_NAN
			)
		.def("compute_delta_V", &sbm::SphericalShaping::computeDeltaV)
		.def("compute_time_of_flight", &sbm::SphericalShaping::computeTimeOfFlight);

	}
}
