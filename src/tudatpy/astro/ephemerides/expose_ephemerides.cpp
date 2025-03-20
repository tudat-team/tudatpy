/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/ephemerides.h>
#include <tudat/simulation/simulation.h>  // TODO: EphemerisType should be in <tudat/astro/ephemerides.h>

namespace py = pybind11;

namespace te = tudat::ephemerides;
namespace tss = tudat::simulation_setup;
namespace ti = tudat::interpolators;

namespace tudatpy
{

namespace trampoline
{

class PyEphemeris : public te::Ephemeris
{
//   public:
//    /* Inherit the constructors */
//    using te::Ephemeris::Ephemeris;
//
//    Eigen::Vector6d getCartesianState( const double seconds_since_epoch ) override
//    {
//        PYBIND11_OVERLOAD_PURE( Eigen::Vector6d,    /* Return type */
//                                te::Ephemeris,      /* Parent class */
//                                getCartesianState,  /* Name of function in C++ (must match
//                                                       Python name) */
//                                seconds_since_epoch /* Argument(s) */
//        );
//    }
//
//    //  Eigen::Matrix<long double, 6, 1> getCartesianLongState(const
    //  double seconds_since_epoch) override {
    //    PYBIND11_OVERLOAD_PURE(
    //        Eigen::Matrix<long double, 6, 1>, /* Return type */
    //        te::Ephemeris,                    /* Parent class */
    //        getCartesianLongState,            /* Name of function in
    //        C++ (must match Python name) */ seconds_since_epoch /*
    //        Argument(s) */
    //    );
    //  }

    //  Eigen::Matrix<long double, 6, 1>
    //  getCartesianStateFromExtendedTime() override {
    //    PYBIND11_OVERLOAD_PURE(
    //        Eigen::Matrix<long double, 6, 1>,  /* Return type */
    //        te::Ephemeris,                     /* Parent class */
    //        getCartesianStateFromExtendedTime, /* Name of function in
    //        C++ (must match Python name) */ const Time & /*
    //        Argument(s) */
    //    );
    //  }

    //  Eigen::Matrix<long double, 6, 1>
    //  getCartesianStateFromExtendedTime() override {
    //    PYBIND11_OVERLOAD_PURE(
    //        Eigen::Matrix<long double, 6, 1>,  /* Return type */
    //        te::Ephemeris,                     /* Parent class */
    //        getCartesianStateFromExtendedTime, /* Name of function in
    //        C++ (must match Python name) */ const double /*
    //        Argument(s) */
    //    );
};

}  // namespace trampoline

PYBIND11_MODULE( expose_ephemerides, m )
{

};
}  // namespace tudatpy
