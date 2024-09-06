/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */



// #include <tudat/astro/low_thrust/lowThrustLeg.h>
// #include <tudat/astro/low_thrust/shape_based/shapeBasedMethod.h>
// #include <tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/getRecommendedBaseFunctionsHodographicShaping.h>

namespace py = pybind11;
namespace tsbm = tudat::shape_based_methods;
// namespace tltt = tudat::low_thrust_trajectories;

namespace tudatpy {
    namespace trajectory_design {
        namespace shape_based_thrust {


            PYBIND11_MODULE(expose_shape_based_thrust, m) {
                //     py::class_<
                //             tltt::LowThrustLeg,
                //             std::shared_ptr<tltt::LowThrustLeg> >(m,
                //             "LowThrustLeg",
                //             .def( "get_trajectory",
                //                   py::overload_cast<
                //                   std::vector< double >& >(
                //                   &tltt::LowThrustLeg::getTrajectory ),
                //                   py::arg("times"),
                //                   )
                //             .def( "get_state",
                //                   &tltt::LowThrustLeg::getStateAtEpoch,
                //                   py::arg("time"),
                //                   )
                //             .def( "compute_delta_v",
                //                   &tltt::LowThrustLeg::computeDeltaV,
                //                   );

                //     py::class_<
                //             tsbm::ShapeBasedMethod,
                //             std::shared_ptr<tsbm::ShapeBasedMethod>,
                //             tltt::LowThrustLeg
                //             >(m, "ShapeBasedMethod",

                //     py::class_<
                //             tsbm::HodographicShaping,
                //             std::shared_ptr<tsbm::HodographicShaping>,
                //             tsbm::ShapeBasedMethod
                //             >(m, "HodographicShaping",
                //             .def(py::init<
                //                  const Eigen::Vector6d&,
                //                  const Eigen::Vector6d&,
                //                  const double,
                //                  const double,
                //                  const int,
                //                  const std::vector< std::shared_ptr<
                //                  tsbm::BaseFunctionHodographicShaping > >&,
                //                  const std::vector< std::shared_ptr<
                //                  tsbm::BaseFunctionHodographicShaping > >&,
                //                  const std::vector< std::shared_ptr<
                //                  tsbm::BaseFunctionHodographicShaping > >&,
                //                  const Eigen::VectorXd&,
                //                  const Eigen::VectorXd&,
                //                  const Eigen::VectorXd& >(),
                //                  py::arg("initial_state"),
                //                  py::arg("final_state"),
                //                  py::arg("time_of_flight"),
                //                  py::arg("central_body_gravitational_parameter"),
                //                  py::arg("number_of_revolutions"),
                //                  py::arg("radial_velocity_functions"),
                //                  py::arg("normal_velocity_functions"),
                //                  py::arg("axial_velocity_functions"),
                //                  py::arg("radial_free_coefficients"),
                //                  py::arg("normal_free_coefficients"),
                //                  py::arg("axial_free_coefficients"),
                //             .def( "get_thrust",
                //                   py::overload_cast< double >(
                //                   &tsbm::HodographicShaping::computeCurrentThrustAcceleration
                //                   ), py::arg( "time_since_departure" ),


                py::class_<
                    tsbm::BaseFunctionHodographicShaping,
                    std::shared_ptr<tsbm::BaseFunctionHodographicShaping> >(
                    m, "BaseFunctionHodographicShaping",
R"doc(Base class for defining settings of the shape functions for hodographic shaping method.

	Base class for defining settings of the shape functions for Hodograph shaping method. Objects derived
	from this class are created by calling the dedicated factory functions in this module

)doc");


                m.def("recommended_radial_hodograph_functions",
                      py::overload_cast<const double>(
                          &tsbm::getRecommendedRadialVelocityBaseFunctions),
                      py::arg("time_of_flight"),
R"doc(Factory function for creating the default radial hodographic trajectory shaping functions.

	Factory function for creating the default radial hodographic trajectory shaping functions. This function
	(and its counterparts normal and axial components) provided three shaping functions that have been found in
	literature to work well for this method. For a given time-of-flight :math:`T`, this function returns a list of
	three shaping functions:

	* Constant term, see :func:`hodograph_constant`
	* Power function, see :func:`hodograph_power`, with exponent = 1.0, scale_factor = :math:`1/T`
	* Power function, see :func:`hodograph_power`, with exponent = 2.0, scale_factor = :math:`1/T`


	:param time_of_flight:
		Total time of flight (in seconds) of the trajectory that is to be generated.
	:return:
		List of default settings object for radial hodographic shaping
)doc");

                m.def("recommended_normal_hodograph_functions",
                      py::overload_cast<const double>(
                          &tsbm::getRecommendedNormalBaseFunctions),
                      py::arg("time_of_flight"),
R"doc(Factory function for creating the default normal hodographic trajectory shaping functions.

	Factory function for creating the default normal hodographic trajectory shaping functions. This function
	(and its counterparts radial and axial components) provided three shaping functions that have been found in
	literature to work well for this method. For a given time-of-flight :math:`T`, this function returns a list of
	three shaping functions:

	* Constant term, see :func:`hodograph_constant`
	* Power function, see :func:`hodograph_power`, with exponent = 1.0, scale_factor = :math:`1/T`
	* Power function, see :func:`hodograph_power`, with exponent = 2.0, scale_factor = :math:`1/T`


	:param time_of_flight:
		Total time of flight (in seconds) of the trajectory that is to be generated.
	:return:
		List of default settings object for axial hodographic shaping
)doc");

                m.def("recommended_axial_hodograph_functions",
                      py::overload_cast<const double, const int>(
                          &tsbm::getRecommendedAxialVelocityBaseFunctions),
                      py::arg("time_of_flight"),
                      py::arg("number_of_revolutions"),
R"doc(Factory function for creating the default axial hodograph	ic trajectory shaping functions.

	Factory function for creating the default axial hodographic trajectory shaping functions. This function
	(and its counterparts radial and normal components) provided three shaping functions that have been found in
	literature to work well for this method. For a given time-of-flight :math:`T` and number of revolutions :math:`N`, this function returns a list of
	three shaping functions:

	* Cosine term, see :func:`hodograph_cosine` with frequency = :math:`\frac{2\pi(N+1/2)}{T}`
	* Power cosine function term, see :func:`hodograph_power_cosine` with  exponent = 3.0, frequency = :math:`\frac{2\pi(N+1/2)}{T}`, scale_factor = :math:`1/T`
	* Power sine function term, see :func:`hodograph_power_sine` with  exponent = 3.0, frequency = :math:`\frac{2\pi(N+1/2)}{T}`, scale_factor = :math:`1/T`


	:param time_of_flight:
		Total time of flight (in seconds) of the trajectory that is to be generated.
	:param number_of_revolutions:
		Number of full revolutions around the central body that are to be used.
	:return:
		List of default settings object for axial hodographic shaping
)doc");


                m.def("hodograph_constant", &tsbm::hodographConstant,
R"doc(Factory function for creating a constant contribution to hodographic trajectory shaping.

	Factory function for creating a constant contribution to hodographic trajectory shaping. This adds a contribution
	:math:`K` to the selected velocity component, with :math:`K` a free parameter.

	:return:
		Settings object for a constant contribution to hodographic shaping.
)doc");

                m.def("hodograph_sine", &tsbm::hodographSine,
                      py::arg("frequency"),
R"doc(Factory function for creating a sine contribution to hodographic trajectory shaping.

	Factory function for creating a sine contribution to hodographic trajectory shaping. For a
	provided frequency :math:`f`, this adds a contribution :math:`K\sin(f\cdot t)` to the selected
	velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


	:param frequency:
		Frequency of the sine contribution to the shape function.
	:return:
		Settings object for a cosine contribution to hodographic shaping.
)doc");

                m.def("hodograph_cosine", &tsbm::hodographCosine,
                      py::arg("frequency"),
R"doc(Factory function for creating a cosine contribution to hodographic trajectory shaping.

	Factory function for creating a cosine contribution to hodographic trajectory shaping. For a
	provided frequency :math:`f`, this adds a contribution :math:`K\cos(f\cdot T)` to the selected
	velocity component, with :math:`T` the time since departure, and :math:`K` a free parameter.


	:param frequency:
		Frequency of the cosine contribution to the shape function.
	:return:
		Settings object for a cosine contribution to hodographic shaping.
)doc");

                m.def("hodograph_exponential", &tsbm::hodographExponential,
                      py::arg("exponent"));

                m.def("hodograph_scaled_exponential",
                      &tsbm::hodographScaledExponential, py::arg("exponent"),
                      py::arg("scale_factor") = 1.0,
"");

                //    m.def("hodograph_scaled_exponential",
                //          &tsbm::hodographScaledExponential,
                //          py::arg("exponent"),
                //          py::arg("scale_factor"));

                m.def("hodograph_exponential_sine",
                      &tsbm::hodographExponentialSine, py::arg("exponent"),
                      py::arg("frequency"));

                m.def(
                    "hodograph_scaled_exponential_sine",
                    &tsbm::hodographScaledExponentialSine, py::arg("exponent"),
                    py::arg("frequency"), py::arg("scale_factor") = 1.0,
"");

                //    m.def("hodograph_scaled_exponential_sine",
                //          &tsbm::hodographScaledExponentialSine,
                //          py::arg("exponent"),
                //          py::arg("frequency"),
                //          py::arg("scale_factor") );

                m.def("hodograph_exponential_cosine",
                      &tsbm::hodographExponentialCosine, py::arg("exponent"),
                      py::arg("frequency"));

                m.def("hodograph_scaled_exponential_cosine",
                      &tsbm::hodographScaledExponentialCosine,
                      py::arg("exponent"), py::arg("frequency"),
                      py::arg("scale_factor") = 1.0,
"");

                //    m.def("hodograph_scaled_exponential_cosine",
                //          &tsbm::hodographScaledExponentialCosine,
                //          py::arg("exponent"),
                //          py::arg("frequency"),
                //          py::arg("scale_factor") );

                m.def("hodograph_power", &tsbm::hodographPower,
                      py::arg("exponent"));

                m.def("hodograph_scaled_power", &tsbm::hodographScaledPower,
                      py::arg("exponent"), py::arg("scale_factor") = 1.0,
"");

                //    m.def("hodograph_scaled_power",
                //          &tsbm::hodographScaledPower,
                //          py::arg("exponent"),
                //          py::arg("scale_factor") );


                m.def("hodograph_power_sine", &tsbm::hodographScaledPowerSine,
                      py::arg("exponent"), py::arg("frequency"),
                      py::arg("scale_factor") = 1.0,
R"doc(Factory function for creating a power sine function contribution to hodographic trajectory shaping.

	Factory function for creating a power sine function contribution to hodographic trajectory shaping. For a
	provided exponent :math:`r`, (optional) scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\cdot c\sin(f\cdot t)\cdot t^{r}` to the selected
	velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


	:param frequency:
		Frequency of the sine contribution to the shape function.
	:param exponent:
		Exponent of the power function contribution to the shape function.
	:param scale_factor:
		Optional scale factor, which can be used to scale the physical meaning of the free parameter :math:`K`.
	:return:
		Settings object for a power sine function contribution to hodographic shaping.
)doc");

                m.def("hodograph_scaled_power_sine",
                      &tsbm::hodographScaledPowerSine, py::arg("exponent"),
                      py::arg("frequency"), py::arg("scale_factor"),
"");

                m.def("hodograph_power_cosine",
                      &tsbm::hodographScaledPowerCosine, py::arg("exponent"),
                      py::arg("frequency"), py::arg("scale_factor") = 1.0,
R"doc(Factory function for creating a power cosine function contribution to hodographic trajectory shaping.

	Factory function for creating a power cosine function contribution to hodographic trajectory shaping. For a
	provided exponent :math:`r`, (optional) scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\cdot c\cos(f\cdot t)\cdot t^{r}` to the selected
	velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


	:param frequency:
		Frequency of the cosine contribution to the shape function.
	:param exponent:
		Exponent of the power function contribution to the shape function.
	:param scale_factor:
		Optional scale factor, which can be used to scale the physical meaning of the free parameter :math:`K`.
	:return:
		Settings object for a power cosine function contribution to hodographic shaping.
)doc");

                m.def("hodograph_scaled_power_cosine",
                      &tsbm::hodographScaledPowerCosine, py::arg("exponent"),
                      py::arg("frequency"), py::arg("scale_factor"),
"");
            }


        }  // namespace shape_based_thrust
    }  // namespace trajectory_design
}  // namespace tudatpy
