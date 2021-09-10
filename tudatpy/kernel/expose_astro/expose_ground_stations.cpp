/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_ground_stations.h"

#include "tudat/astro/ground_stations/groundStation.h"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

namespace tgs = tudat::ground_stations;

namespace py = pybind11;

namespace tudatpy {

void expose_ground_stations(py::module &m) {


    py::class_<tgs::GroundStation,
               std::shared_ptr<tgs::GroundStation>>(m, "GroundStation")
            .def_property_readonly("pointing_angles_calculator", &tgs::GroundStation::getPointingAnglesCalculator );

    py::class_<tgs::PointingAnglesCalculator,
               std::shared_ptr<tgs::PointingAnglesCalculator>>(m, "PointingAnglesCalculator")
            .def("calculate_elevation_angle", &tgs::PointingAnglesCalculator::calculateElevationAngle,
                 py::arg( "inertial_vector_to_target" ),
                 py::arg( "time" ) )
            .def("calculate_azimuth_angle", &tgs::PointingAnglesCalculator::calculateAzimuthAngle,
                 py::arg( "inertial_vector_to_target" ),
                 py::arg( "time" ) )
            .def("convert_inertial_vector_to_topocentric",
                 &tgs::PointingAnglesCalculator::convertVectorFromInertialToTopocentricFrame,
                 py::arg( "inertial_vector" ),
                 py::arg( "time" ) );
}

}// namespace tudatpy
