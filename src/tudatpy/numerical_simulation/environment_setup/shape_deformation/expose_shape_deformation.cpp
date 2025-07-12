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
#include "expose_shape_deformation.h"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/environment_setup/createBodyDeformationModel.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace shape_deformation
{

void expose_shape_deformation_setup( py::module &m )
{
    py::class_< tss::BodyDeformationSettings, std::shared_ptr< tss::BodyDeformationSettings > >(
            m,
            "BodyDeformationSettings",
            R"doc(

         Base class for providing settings for body shape deformation model.

         Functional (base) class for settings of body shape deformation models that require no information in addition to their type.
         Body shape deformation model settings requiring additional information must be defined using an object derived from this class.


      )doc" );

    py::class_< tss::BasicSolidBodyDeformationSettings,
                std::shared_ptr< tss::BasicSolidBodyDeformationSettings >,
                tss::BodyDeformationSettings >( m,
                                                "BasicSolidBodyDeformationSettings",
                                                R"doc(

         Class for defining model settings for simple tidal solid-body shape deformation.

         `BodyDeformationSettings` derived class for simple tidal solid-body shape deformation.




      )doc" );

    m.def( "basic_solid_body_tidal",
           &tss::basicTidalBodyShapeDeformation,
           py::arg( "tide_raising_bodies" ),
           py::arg( "displacement_love_numbers" ),
           py::arg( "reference_radius" ) = TUDAT_NAN,
           R"doc(

 Function for creating basic tidal solid-body shape deformation

 Function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{m}` and :math:`l_{m}` (with only :math:`m=2,3` currently supported). This function implements equations (7.5) and (7.6) of the `IERS Conventions 2010 <https://iers-conventions.obspm.fr/conventions_material.php>`_.


 Parameters
 ----------
 tide_raising_bodies : list[ string ]
     List of bodies that raise a tide that induces the shape variation.
 displacement_love_numbers : dict[ int, [float,float] ]
     Dictionary of pairs. The dictionary key the spherical harmonic degree :math:`l` of the tidal deformation (2 or 3 are currently supported). The dictionary value is comprised of a pair of floats representing the :math:`h_{2}` and :math:`l_{2}` deformation Love numbers
 reference_radius : float, default = NaN
     Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
 Returns
 -------
 BasicSolidBodyDeformationSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class





 Examples
 --------
 In this example, we create a settings for degree 2 tidal deformation of the Earth due to the Sun and Moon:

 .. code-block:: python

   # Create Love numbers
   love_numbers = dict()
   love_numbers[2] = (0.6, 0.08)

   # Create tide raising bodies
   tide_raising_bodies = ["Sun", "Moon"]

   # Append shape variation settings to existing (default is empty) list
   body_settings.get( "Earth" ).shape_deformation_settings.append(
       environment_setup.shape_deformation.basic_solid_body_tidal(
           tide_raising_bodies, love_numbers ) )


     )doc" );

    m.def( "degree_two_basic_solid_body_tidal",
           &tss::degreeTwoBasicTidalBodyShapeDeformation,
           py::arg( "tide_raising_bodies" ),
           py::arg( "love_number" ),
           py::arg( "shida_number" ),
           py::arg( "reference_radius" ) = TUDAT_NAN,
           R"doc(

 Function for creating degree 2 basic tidal solid-body shape deformation

 Function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{2}` and :math:`l_{2}`. This function implements equations (7.5) of the IERS Conventions 2010, and provides a simplified interface (for degree 2 only) of :func:`basic_solid_body_tidal`.


 Parameters
 ----------
 tide_raising_bodies : list[ string ]
     List of bodies that raise a tide that induces the shape variation.
 love_number : float
     Value of :math:`h_{2}` deformation Love number`
 shida_number : float
     Value of :math:`l_{2}` deformation Shida number`
 reference_radius : float, default = NaN
     Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
 Returns
 -------
 BasicSolidBodyDeformationSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class





 Examples
 --------
 In this example, we create a settings for degree 2 tidal deformation of the Earth due to the Sun and Moon:

 .. code-block:: python

   # Create Love numbers
   h2_love_number = 0.6
   l2_shida_number = 0.08

   # Create tide raising bodies
   tide_raising_bodies = ["Sun", "Moon"]

   # Append shape variation settings to existing (default is empty) list
   body_settings.get( "Earth" ).shape_deformation_settings.append(
       environment_setup.shape_deformation.degree_two_basic_solid_body_tidal(
           tide_raising_bodies, h2_love_number, l2_shida_number ) )


     )doc" );

    m.def( "iers_2010_solid_body_tidal",
           &tss::iers2010TidalBodyShapeDeformation,
           R"doc(

 Function for creating full IERS 2010 shape deformation model

 Function for creating full IERS 2010 shape deformation model, computing the tidal shape variation due to the full model defined in
 Section 7.1.1 of the IERS Conventions 2010, implementing Eqs. (7.5)-(7.13), including all terms from Tables 7.3a and 7.3b. At present, none of the input parameters of the model can be varied.

 Returns
 -------
 BodyDeformationSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BodyDeformationSettings` defining the IERS 2010 settings

     )doc" );

    m.def( "pole_tidal", &tss::poleTideBodyShapeDeformation, R"doc(No documentation found.)doc" );

    m.def( "ocean_tidal",
           &tss::oceanTideBodyShapeDeformation,
           py::arg( "blq_files" ),
           R"doc(No documentation found.)doc" );
}

}  // namespace shape_deformation
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
