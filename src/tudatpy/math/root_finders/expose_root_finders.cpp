/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/pybind11.h>

#include "tudat/math/root_finders.h"

namespace py = pybind11;
namespace trf = tudat::root_finders;

namespace tudatpy {
    namespace math {

        namespace root_finders {

            PYBIND11_MODULE(expose_root_finders, m) {
                /*
                 *
                 *
                 *  root_finders
                 *  ├── bisection.h
                 *  ├── createRootFinder.h
                 *  ├── halleyRootFinder.h
                 *  ├── newtonRaphson.h
                 *  ├── rootFinder.h
                 *  ├── secantRootFinder.h
                 *  └── terminationConditions.h
                 *
                 */

                py::enum_<trf::MaximumIterationHandling>(
                    m, "MaximumIterationHandling",
                    R"doc(

        Enumeration of types of behaviour to be used when the convergence criterion on maximum number of iterations is reached.





     )doc")
                    .value("accept_result",
                           trf::MaximumIterationHandling::accept_result,
                           R"doc(
The program will accept the root at the final iteration, without any additional output
     )doc")
                    .value("accept_result_with_warning",
                           trf::MaximumIterationHandling::
                               accept_result_with_warning,
                           R"doc(
The program will accept the root at the final iteration, but will print a warning to the terminal that the root finder may not have converged
     )doc")
                    .value("throw_exception",
                           trf::MaximumIterationHandling::throw_exception,
                           R"doc(
The program will not accept the root at the final iteration, and will throw an exception
     )doc")
                    .export_values();


                py::class_<trf::RootFinder<double>,
                           std::shared_ptr<trf::RootFinder<double>>>(
                    m, "RootFinderCore");

                py::class_<trf::NewtonRaphson<double>,
                           std::shared_ptr<trf::NewtonRaphson<double>>,
                           trf::RootFinder<double>>(m, "NewtonRaphsonCore")
                    .def(py::init<const double, const unsigned int>(),
                         py::arg("x_tol"), py::arg("max_iter"));

                py::class_<trf::RootFinderSettings,
                           std::shared_ptr<trf::RootFinderSettings>>(
                    m, "RootFinderSettings",
                    R"doc(

        Class to define settings for a root finder.





     )doc");

                m.def("bisection", &trf::bisectionRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
                      R"doc(

Function to create settings for a bisection root-finder.

Function to create settings for a bisection root finder. This root finder approximates the root by initializing with
two initial guesses :math:`x_{\downarrow,0}` and :math:`x_{\uparrow,0}`, for which it is required that
:math:`f(x_{\downarrow,0}) < 0` and :math:`f(x_{\uparrow,0}) > 0`. At each iteration :math:`i`, the current guess of
the root :math:`x_{i}` is:

.. math::
   x_{i}=\begin{cases}
   x_{\downarrow,i}, & |f(x_{\downarrow,i})|<|f(x_{\uparrow,i})|\\
	  x_{\uparrow,i}, & |f(x_{\downarrow,i})|\ge|f(x_{\uparrow,i})|
	       \end{cases}

The midpoint :math:`x_{m,i}` of :math:`x_{\downarrow,i}` and :math:`x_{\uparrow,i}` is then computed from :math:`x_{m,i}=(x_{\downarrow,i}-x_{\uparrow,i})/2`.
Depending on the sign of :math:`f(x_{m,i})`, it then replaces either :math:`x_{\downarrow,i}` or :math:`x_{\uparrow,i}` (depending on whether
its sign matches :math:`f(x_{\downarrow,i})` for iteration :math:`i+1` and :math:`f(x_{\uparrow,i})`), while the other point from iteration :math:`i` is retained.

Although slow, the algorithm is ensured to converge to a root, if the two initial guesses indeed have opposite signs (if not, an exception is thrown).


Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Bisection root-finding settings object






    )doc");


                m.def("newton_raphson", &trf::newtonRaphsonRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
                      R"doc(

Function to create settings for a Newton-Raphson root-finder.

Function to create settings for a bisection root finder. This root finder approximates the root by initializing with
a single initial guesses :math:`x_{0}` and requires an analytical formulation for :math:`f(x)` and :math:`f'(x)=\frac{d}{dx}f(x)`.
The algorithm uses the following equation to iterate:

.. math::
   x_{i+1}=x_{i}-\frac{f(x_{i})}{f'(x_{i})}


Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Newton-Raphson root-finding settings object






    )doc");

                m.def("halley", &trf::halleyRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
                      R"doc(

Function to create settings for a Halley root-finder.

Function to create settings for a Halley root finder. This root finder approximates the root by initializing with
a single initial guesses :math:`x_{0}` and requires an analytical formulation for :math:`f(x)`, :math:`f'(x)=\frac{d}{dx}f(x)` and :math:`f''(x)=\frac{d^{2}}{dx^{2}}f(x)`.
The algorithm uses the following equation to iterate:

.. math::
   x_{i+1}=x_{i}-\frac{2f(x_{i})f'(x_{i})}{2(f'(x_{i}))^{2}-f(x_{i})f''(x_{i})}




Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Halley root-finding settings object






    )doc");

                m.def("secant", &trf::secantRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
                      R"doc(

Function to create settings for a secant method root-finder.

Function to create settings for a root finder using the secant method. This root finder approximates the root by initializing with
two initial guesses :math:`x_{0}` and :math:`x_{1}`. The algorithm uses the following equation to iterate:

.. math::
   x_{i+1}=x_{i}-f(x_{i})\frac{x_{i}-x_{i-1}}{f(x_{i})-f(x_{i-1})}



Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Secant root-finding settings object






    )doc");
            }

        }  // namespace root_finders
    }  // namespace math

}  // namespace tudatpy
