/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_geometry.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/basicTypedefs.h>

#include "tudat/math/geometric/capsule.h"

namespace py = pybind11;

namespace tgs = tudat::geometric_shapes;

namespace tudatpy
{

void expose_geometry( py::module& m )
{
    py::class_< tudat::SurfaceGeometry, std::shared_ptr< tudat::SurfaceGeometry > >(
            m, "SurfaceGeometry", R"doc(Base representation of a parameterized surface used for vehicle geometry calculations.)doc" );

    py::class_< tgs::CompositeSurfaceGeometry, std::shared_ptr< tgs::CompositeSurfaceGeometry >, tudat::SurfaceGeometry >(
            m, "CompositeSurfaceGeometry", R"doc(Surface geometry assembled from multiple component surfaces.)doc" );

    py::class_< tgs::Capsule, std::shared_ptr< tgs::Capsule >, tgs::CompositeSurfaceGeometry >( m,
                                                                                                "Capsule",
                                                                                                R"doc(

         Composite capsule geometry defined by a rounded nose, middle radius, rear section and shoulder radius.

      )doc" )
            .def( py::init< const double, const double, const double, const double, const double >( ),
                  py::arg( "nose_radius" ),
                  py::arg( "middle_radius" ),
                  py::arg( "rear_length" ),
                  py::arg( "rear_angle" ),
                  py::arg( "side_radius" ),
                  R"doc(

         Create a capsule from nose, middle and side radii and rear length in metres, and the rear angle in radians.

         Parameters
         ----------
         nose_radius : float
             Radius of the capsule nose, in metres.
         middle_radius : float
             Radius of the capsule's middle section, in metres.
         rear_length : float
             Length of the rear section, in metres.
         rear_angle : float
             Angle of the rear conical section, in radians.
         side_radius : float
             Radius of curvature of the capsule side, in metres.

      )doc" )
            .def_property_readonly( "middle_radius",
                                    &tgs::Capsule::getMiddleRadius,
                                    R"doc(

         **read-only**

         Radius of the capsule middle section, in metres.

         :type: float

      )doc" )
            .def_property_readonly( "volume",
                                    &tgs::Capsule::getVolume,
                                    R"doc(

         **read-only**

         Volume enclosed by the capsule geometry, in cubic metres.

         :type: float

      )doc" )
            .def_property_readonly( "length",
                                    &tgs::Capsule::getLength,
                                    R"doc(

         **read-only**

         Total axial length of the capsule geometry, in metres.

         :type: float

      )doc" );
};

}  // namespace tudatpy
