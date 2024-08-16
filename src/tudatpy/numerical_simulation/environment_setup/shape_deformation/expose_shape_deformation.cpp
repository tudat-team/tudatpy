/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/environment_setup/createBodyDeformationModel.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace shape_deformation {

                PYBIND11_MODULE(expose_shape_deformation, m) {
                    py::class_<tss::BodyDeformationSettings,
                               std::shared_ptr<tss::BodyDeformationSettings>>(
                        m, "BodyDeformationSettings",
R"doc(Base class for providing settings for body shape deformation model.

	Functional (base) class for settings of body shape deformation models that require no information in addition to their type.
	Body shape deformation model settings requiring additional information must be defined using an object derived from this class.
	
)doc");

                    py::class_<
                        tss::BasicSolidBodyDeformationSettings,
                        std::shared_ptr<tss::BasicSolidBodyDeformationSettings>,
                        tss::BodyDeformationSettings>(
                        m, "BasicSolidBodyDeformationSettings",
R"doc(Class for defining model settings for simple tidal solid-body shape deformation.

	`BodyDeformationSettings` derived class for simple tidal solid-body shape deformation.
)doc");


                    m.def("basic_solid_body_tidal",
                          &tss::basicTidalBodyShapeDeformation,
                          py::arg("tide_raising_bodies"),
                          py::arg("displacement_love_numbers"),
                          py::arg("reference_radius") = TUDAT_NAN,
R"doc(Factory function for creating basic tidal solid-body shape deformation

	Factory function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{m}` and :math:`l_{m}` (with only :math:`m=2,3` currently supported). This function implements equations (7.5) and (7.6) of the IERS 2010 Conventions.
	

	:param tide_raising_bodies:
		List of bodies that raise a tide that induces the shape variation.
	:param displacement_love_numbers:
		Dictionary of pairs. The dictionary key the spherical harmonic degree :math:`l` of the tidal deformation (2 or 3 are currenty supported). The dictionary value is comprised of a pair of floats representing the :math:`h_{2}` and :math:`l_{2}` deformation Love numbers
	:param reference_radius:
		Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyDeformationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.BasicSolidBodyDeformationSettings` class
)doc");

                    m.def("degree_two_basic_solid_body_tidal",
                          &tss::degreeTwoBasicTidalBodyShapeDeformation,
                          py::arg("tide_raising_bodies"),
                          py::arg("love_number"), py::arg("shida_number"),
                          py::arg("reference_radius") = TUDAT_NAN,
R"doc(Factory function for creating degree 2 basic tidal solid-body shape deformation

	Factory function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{2}` and :math:`l_{2}`. This function implements equations (7.5) of the IERS 2010 Conventions, and provides a simplified interface (for degree 2 only) of :func:`basic_solid_body_tidal`.
	

	:param tide_raising_bodies:
		List of bodies that raise a tide that induces the shape variation.
	:param love_number:
		Value of :math:`h_{2}` deformation Love number`
	:param shida_number:
		Value of :math:`l_{2}` deformation Shida number`
	:param reference_radius:
		Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyDeformationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.BasicSolidBodyDeformationSettings` class
)doc");

                    m.def("iers_2010_solid_body_tidal",
                          &tss::iers2010TidalBodyShapeDeformation,
R"doc(Factory function for creating full IERS 2010 shape deformation model

	Factory function for creating full IERS 2010 shape deformation model, computing the tidal shape variation due to the full model defined in Section 7.1.1 of the 2010 IERS conventions, implementing Eqs. (7.5)-(7.13), including all terms from Tables 7.3a and 7.3b. At present, none of the input parameters of the model can be varid.
	
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyDeformationSettings` defining the IERS 2010 settings
)doc");
                }

            }  // namespace shape_deformation
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
