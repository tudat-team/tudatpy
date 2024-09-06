/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>


// #include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tpc = tudat::physical_constants;

namespace tudat {

    namespace simulation_setup {

        inline std::shared_ptr<GravityFieldSettings>
        fromFileSphericalHarmonicsGravityFieldSettings(
            const std::string& filePath, const int maximumDegree,
            const int maximumOrder,
            const std::string& associatedReferenceFrame = "",
            const int gravitationalParameterIndex = 0,
            const int referenceRadiusIndex = 1) {
            return std::make_shared<
                FromFileSphericalHarmonicsGravityFieldSettings>(
                filePath, associatedReferenceFrame, maximumDegree, maximumOrder,
                gravitationalParameterIndex, referenceRadiusIndex);
        }


        inline std::shared_ptr<GravityFieldSettings>
        predefinedSphericalHarmonic(
            const SphericalHarmonicsModel sphericalHarmonicsModel,
            const int maximumDegree = -1) {
            return std::make_shared<
                FromFileSphericalHarmonicsGravityFieldSettings>(
                sphericalHarmonicsModel, maximumDegree);
        }

        inline std::shared_ptr<GravityFieldSettings>
        createHomogeneousTriAxialEllipsoidGravitySettingsDeprecated(
            const double axisA, const double axisB, const double axisC,
            const double ellipsoidDensity, const int maximumDegree,
            const int maximumOrder, const std::string& associatedReferenceFrame,
            const double gravitationalConstant = tpc::GRAVITATIONAL_CONSTANT) {
            static bool isWarningPrinted = false;
            if(isWarningPrinted == false) {
                tudat::utilities::printDeprecationWarning(
                    "tudatpy.numerical_simulation.environment_setup.gravity_"
                    "field.spherical_harmonic_triaxial_body",
                    "tudatpy.numerical_simulation.environment_setup.gravity_"
                    "field.sh_triaxial_ellipsoid_from_density");
                isWarningPrinted = true;
            }

            return createHomogeneousTriAxialEllipsoidGravitySettings(
                axisA, axisB, axisC, ellipsoidDensity, maximumDegree,
                maximumOrder, associatedReferenceFrame, gravitationalConstant);
        }

    }  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace gravity_field {

                PYBIND11_MODULE(expose_gravity_field, m) {
                    /////////////////////////////////////////////////////////////////////////////
                    // createGravityField.h
                    /////////////////////////////////////////////////////////////////////////////
                    py::enum_<tss::GravityFieldType>(
                        m, "GravityFieldType",
                        R"doc(Enumeration of gravity field types.

	Enumeration of gravity field types supported by tudat.


	:member polyhedron:
	:member central_gravity:
	:member central_spice_gravity:
	:member spherical_harmonic_gravity:
	:member polyhedron_gravity:
)doc")
                        .value("central_gravity",
                               tss::GravityFieldType::central, "")
                        .value("central_spice_gravity",
                               tss::GravityFieldType::central_spice, "")
                        .value("spherical_harmonic_gravity",
                               tss::GravityFieldType::spherical_harmonic, "")
                        .value("polyhedron_gravity",
                               tss::GravityFieldType::polyhedron, "")
                        .value("ring_gravity",
                               tss::GravityFieldType::one_dimensional_ring, "")
                        .export_values();


                    py::enum_<tss::SphericalHarmonicsModel>(
                        m, "PredefinedSphericalHarmonicsModel",
                        R"doc(Enumeration of predefined spherical harmonics models.

	Enumeration of predefined spherical harmonics models supported by tudat, for which thee coefficient files are automatically available (downloaded from
	`here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_). The directory where these files are stored can be
	extracted using the :func:`~tudatpy.io.get_gravity_models_path` function.


	:member egm96:
	:member ggm02c:
	:member ggm02s:
	:member goco05c:
	:member glgm3150:
	:member lpe200:
	:member gggrx1200:
	:member jgmro120d:
	:member jgmess160a:
	:member shgj180u:
)doc")
                        .value("egm96", tss::SphericalHarmonicsModel::egm96, "")
                        .value("ggm02c", tss::SphericalHarmonicsModel::ggm02c,
                               "")
                        .value("ggm02s", tss::SphericalHarmonicsModel::ggm02s,
                               "")
                        .value("goco05c", tss::SphericalHarmonicsModel::ggm02s,
                               "")
                        .value("glgm3150",
                               tss::SphericalHarmonicsModel::glgm3150, "")
                        .value("lpe200", tss::SphericalHarmonicsModel::lpe200,
                               "")
                        .value("gggrx1200",
                               tss::SphericalHarmonicsModel::lpe200, "")
                        .value("jgmro120d",
                               tss::SphericalHarmonicsModel::jgmro120d, "")
                        .value("jgmess160a",
                               tss::SphericalHarmonicsModel::jgmess160a, "")
                        .value("shgj180u",
                               tss::SphericalHarmonicsModel::shgj180u, "")
                        .export_values();


                    py::class_<tss::GravityFieldSettings,
                               std::shared_ptr<tss::GravityFieldSettings> >(
                        m, "GravityFieldSettings",
                        R"doc(Base class for providing settings for automatic gravity field model creation.

	This class is a functional base class for settings of gravity field models that require no information in addition to their type.
	Gravity field model classes requiring additional information must be created using an object derived from this class.

)doc")
                        .def_property_readonly(
                            "gravity_field_type",
                            &tss::GravityFieldSettings::getGravityFieldType,
                            R"doc(Type of gravity field model that is to be created.
	)doc");


                    py::class_<
                        tss::CentralGravityFieldSettings,
                        std::shared_ptr<tss::CentralGravityFieldSettings>,
                        tss::GravityFieldSettings>(
                        m, "CentralGravityFieldSettings",
                        R"doc(`GravityFieldSettings` derived class defining settings of point mass gravity field.

	Derived class of `GravityFieldSettings` for central gravity fields, which are defined by a single gravitational parameter.

)doc")
                        //            .def(py::init<double>(),
                        //            py::arg("gravitational_parameter") )
                        .def_property(
                            "gravitational_parameter",
                            &tss::CentralGravityFieldSettings::
                                getGravitationalParameter,
                            &tss::CentralGravityFieldSettings::
                                resetGravitationalParameter,
                            R"doc(Gravitational parameter of central gravity field.
	)doc");


                    py::class_<tss::SphericalHarmonicsGravityFieldSettings,
                               std::shared_ptr<
                                   tss::SphericalHarmonicsGravityFieldSettings>,
                               tss::GravityFieldSettings>(
                        m, "SphericalHarmonicsGravityFieldSettings",
                        R"doc(`GravityFieldSettings` derived class defining settings of spherical harmonic gravity field representation.

	Derived class of `GravityFieldSettings` for gravity fields, which are defined by a spherical harmonic gravity field representation.

)doc")
                        //            .def(py::init<const double, const double,
                        //            const Eigen::MatrixXd, const
                        //            Eigen::MatrixXd, const std::string&>(),
                        //                 py::arg("gravitational_parameter"),
                        //                 py::arg("reference_radius"),
                        //                 py::arg("cosine_coefficients"),
                        //                 py::arg("sine_coefficients"),
                        //                 py::arg("associated_reference_frame"))
                        .def_property(
                            "gravitational_parameter",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getGravitationalParameter,
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                resetGravitationalParameter,
                            R"doc(Gravitational parameter of gravity field.
	)doc")
                        .def_property(
                            "normalized_cosine_coefficients",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getCosineCoefficients,
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                resetCosineCoefficients,
                            R"doc(Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
	)doc")
                        .def_property(
                            "normalized_sine_coefficients",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getSineCoefficients,
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                resetSineCoefficients,
                            R"doc(Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
	)doc")
                        .def_property(
                            "associated_reference_frame",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getAssociatedReferenceFrame,
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                resetAssociatedReferenceFrame,
                            R"doc(Identifier for body-fixed reference frame with which the coefficients are associated.
	)doc")
                        .def_property(
                            "create_time_dependent_field",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getCreateTimeDependentField,
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                setCreateTimeDependentField,
                            R"doc(Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed initially).
	)doc")
                        .def_property(
                            "scaled_mean_moment_of_inertia",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getScaledMeanMomentOfInertia,
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                setScaledMeanMomentOfInertia,
                            "")
                        .def_property_readonly(
                            "reference_radius",
                            &tss::SphericalHarmonicsGravityFieldSettings::
                                getReferenceRadius,
                            R"doc(Reference radius of spherical harmonic field expansion.
	)doc");


                    py::class_<
                        tss::FromFileSphericalHarmonicsGravityFieldSettings,
                        std::shared_ptr<
                            tss::
                                FromFileSphericalHarmonicsGravityFieldSettings>,
                        tss::SphericalHarmonicsGravityFieldSettings>(
                        m, "FromFileSphericalHarmonicsGravityFieldSettings",
                        "");


                    py::class_<
                        tss::PolyhedronGravityFieldSettings,
                        std::shared_ptr<tss::PolyhedronGravityFieldSettings>,
                        tss::GravityFieldSettings>(
                        m, "PolyhedronGravityFieldSettings",
                        R"doc(`GravityFieldSettings` derived class defining settings of a polyhedron gravity field representation.

	Derived class of `GravityFieldSettings` for gravity fields, which are defined by a polyhedron gravity field representation.

)doc")
                        .def_property(
                            "gravitational_parameter",
                            &tss::PolyhedronGravityFieldSettings::
                                getGravitationalParameter,
                            &tss::PolyhedronGravityFieldSettings::
                                resetGravitationalParameter,
                            R"doc(Gravitational parameter of gravity field.
	)doc")
                        .def_property(
                            "density",
                            &tss::PolyhedronGravityFieldSettings::getDensity,
                            &tss::PolyhedronGravityFieldSettings::resetDensity,
                            R"doc(Density of the polyhedron.
	)doc")
                        .def_property(
                            "associated_reference_frame",
                            &tss::PolyhedronGravityFieldSettings::
                                getAssociatedReferenceFrame,
                            &tss::PolyhedronGravityFieldSettings::
                                resetAssociatedReferenceFrame,
                            R"doc(Identifier for body-fixed reference frame with which the vertices coordinates are associated.
	)doc")
                        .def_property_readonly(
                            "vertices_coordinates",
                            &tss::PolyhedronGravityFieldSettings::
                                getVerticesCoordinates,
                            R"doc(Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
row per vertex, 3 columns).

	)doc")
                        .def_property_readonly(
                            "vertices_defining_each_facet",
                            &tss::PolyhedronGravityFieldSettings::
                                getVerticesDefiningEachFacet,
                            R"doc(Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
when seen from the outside of the polyhedron.

	)doc");


                    m.def(
                        "central", &tss::centralGravitySettings,
                        py::arg("gravitational_parameter"),
                        R"doc(Factory function for central gravity field settings object.

	Factory function for settings object, defining a point-mass gravity field model with user-defined gravitational parameter :math:`\mu`. The gravitational potential is the defined as:

	.. math::
	   U(\mathbf{r})=\frac{\mu}{||\mathbf{r}||}

	with :math:`\mathbf{r}` the position vector measured from the body's center of mass.


	:param gravitational_parameter:
		Gravitational parameter defining the point-mass gravity field.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.CentralGravityFieldSettings` class
)doc");

                    m.def(
                        "central_spice", &tss::centralGravityFromSpiceSettings,
                        R"doc(Factory function to create central gravity field settings from Spice settings.

	Factory function for settings object, defining a point-mass gravity field model. This function provides the same model as :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.central`), but with gravitational parameter :math:`\mu` from Spice.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` class of gravity field type ``central_spice``"
)doc");

                    m.def(
                        "spherical_harmonic",
                        py::overload_cast<
                            const double, const double, const Eigen::MatrixXd,
                            const Eigen::MatrixXd, const std::string&>(
                            &tss::sphericalHarmonicsGravitySettings),
                        py::arg("gravitational_parameter"),
                        py::arg("reference_radius"),
                        py::arg("normalized_cosine_coefficients"),
                        py::arg("normalized_sine_coefficients"),
                        py::arg("associated_reference_frame"),
                        R"doc(Factory function for creating a spherical harmonics gravity field settings object.

	Factory function for settings object, defining a gravity field model through spherical harmonic expansion.
	The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
	It represents the frame in which the spherical harmonic field is defined.

	The gravitational potential is the defined as:

	.. math::
	   U(\mathbf{r})=\sum_{l=0}^{l_{max}}\sum_{m=0}^{l}\mu\left(\frac{{R}^{l}}{r^{l+1}}\right)\bar{P}_{lm}(\sin\phi)\left(\bar{C}_{lm}\cos m\theta+\bar{S}_{lm}\sin m\theta\right)

	with :math:`\mathbf{r}` the position vector of the evaluation point, measured from the body's center of mass. The angles :math:`\phi` and :math:`\theta` are the body-fixed latitude and longitude of the evaluation point, and :math:`\bar{P}_{lm}` is the associated Legendre polynomial (at degree/order :math`l/m`).

	Note: Spherical harmonic coefficients used for this environment model must *always* be fully normalized.
	To normalize un-normalized spherical harmonic coefficients, see :func:`~tudatpy.astro.gravitation.normalize_spherical_harmonic_coefficients`.


	:param gravitational_parameter:
		Gravitational parameter :math:`\mu` of gravity field.
	:param reference_radius:
		Reference radius :math:`R` of spherical harmonic field expansion.
	:param normalized_cosine_coefficients:
		Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\bar{C}_{ij}` at degree i and order j.
		As such, note that entry (0,0) of cosine coefficients should be equal to 1.

	:param normalized_sine_coefficients:
		Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\bar{S}_{ij}`
		at degree i and order j.

	:param associated_reference_frame:
		Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
)doc");

                    m.def(
                        "from_file_spherical_harmonic",
                        tss::fromFileSphericalHarmonicsGravityFieldSettings,
                        py::arg("file"), py::arg("maximum_degree"),
                        py::arg("maximum_order"),
                        py::arg("associated_reference_frame") = "",
                        py::arg("gravitational_parameter_index") = 0,
                        py::arg("reference_radius_index") = 1,
                        R"doc(Factory function to load a custom spherical harmonics gravity field settings from a file.

	Factory function to load a custom spherical harmonics gravity field settings from a file. The file should contain **fully normalized** spherical harmonic coefficients.
	The associated gravitational paramerer and reference radius should be given in m^3/s^2 and m, respectively. The file format should be the same as that used for the files
	in the directories `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_. Specifically, the file should contain

	- The first line should be a series of text blocks (typically numerical data). Two of these blocks (by default the first and second one) should be the gravitational parameter and reference radius, respectively. The text block should be separated by spaces, tabs and/or commas
	- Each subsequent line should contain a set of spherical harmonic coefficients (first ordered in ascending order by degree, then in ascending order by order), where the first, second, third and fourth value of the line should be: degree :math:`l`, order :math:`m`, normalized cosine coefficient :math:`\bar{C}_{lm}`, normalized sine coefficient :math:`\bar{S}_{lm}`. Additional entries (for instance with coefficient uncertainties) are ignored.


	:param file:
		Full file path and name where th gravity field file is located
	:param maximum_degree:
		Maximum degree of the coefficients that are to be loaded
	:param maximum_order:
		Maximum order of the coefficients that are to be loaded
	:param associated_reference_frame:
		Name of the body-fixed reference frame to which the gravity field is to be fixed. If left empty, this reference frame will automatically be set to the body-fixed frame defined by this body's rotation (see :ref:`\`\`rotation_model\`\`` for specifying rotation models).
	:param gravitational_parameter_index:
		Index of the values in the file header (first line of file) that contains the gravitational parameter
	:param reference_radius_index:
		Index of the values in the file header (first line of file) that contains the reference radius
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
)doc");

                    m.def(
                        "predefined_spherical_harmonic",
                        tss::predefinedSphericalHarmonic,
                        py::arg("predefined_model"),
                        py::arg("maximum_degree") = -1,
                        R"doc(Factory function for spherical harmonics gravity field settings of a predefined model.

	Factory function for spherical harmonics gravity field settings of a predefined model


	:param predefined_model:
		Identified for gravity field model that is to be loaded
	:param maximum_degree:
		Maximum degree and order to which the coefficients are to be loaded. If value is negative, all coefficients for the specified gravity field are loaded
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
)doc");


                    m.def(
                        "polyhedron_from_mu",
                        py::overload_cast<const double, const Eigen::MatrixXd,
                                          const Eigen::MatrixXi,
                                          const std::string&, const double>(
                            &tss::polyhedronGravitySettingsFromMu),
                        py::arg("gravitational_parameter"),
                        py::arg("vertices_coordinates"),
                        py::arg("vertices_defining_each_facet"),
                        py::arg("associated_reference_frame"),
                        py::arg("gravitational_constant") =
                            tpc::GRAVITATIONAL_CONSTANT,
                        R"doc(Factory function for creating a polyhedron gravity field settings object, using the gravitational parameter.

	Factory function for settings object, defining a gravity field model through a polyhedron.
	The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
	It represents the frame in which the polyhedron field is defined.

	The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
	to Werner and Scheeres [2]_.

	This function uses the gravitational parameter to define the gravity field. To instead use the density
	constant see :func:`~tudatpy.astro.gravitation.polyhedron_from_density`.


	:param gravitational_parameter:
		Gravitational parameter :math:`\mu` of gravity field.
	:param vertices_coordinates:
		Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
		row per vertex, 3 columns).

	:param vertices_defining_each_facet:
		Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
		the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
		when seen from the outside of the polyhedron.

	:param associated_reference_frame:
		Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_constant:
		Newton's gravitational constant G, used to computed the density

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class
)doc");

                    m.def(
                        "polyhedron_from_density",
                        py::overload_cast<const double, const Eigen::MatrixXd,
                                          const Eigen::MatrixXi,
                                          const std::string&, const double>(
                            &tss::polyhedronGravitySettings),
                        py::arg("density"), py::arg("vertices_coordinates"),
                        py::arg("vertices_defining_each_facet"),
                        py::arg("associated_reference_frame"),
                        py::arg("gravitational_constant") =
                            tpc::GRAVITATIONAL_CONSTANT,
                        R"doc(Factory function for creating a polyhedron gravity field settings object, using the density.

	Factory function for settings object, defining a gravity field model through a polyhedron.
	The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
	It represents the frame in which the polyhedron field is defined.

	The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
	to Werner and Scheeres [2]_.

	This function uses the density to define the gravity field. To instead use the
	gravitational parameter see :func:`~tudatpy.astro.gravitation.polyhedron_from_mu`.


	:param density:
		Density of the polyhedron.

	:param vertices_coordinates:
		Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
		row per vertex, 3 columns).

	:param vertices_defining_each_facet:
		Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
		the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
		when seen from the outside of the polyhedron.

	:param associated_reference_frame:
		Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_constant:
		Newton's gravitational constant G, used to computed the gravitational parameter

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class
)doc");

                    // Triaxial ellipsoid: overload 1
                    m.def(
                        "sh_triaxial_ellipsoid_from_density",
                        py::overload_cast<const double, const double,
                                          const double, const double, const int,
                                          const int, const std::string&,
                                          const double>(
                            &tss::
                                createHomogeneousTriAxialEllipsoidGravitySettings),
                        py::arg("axis_a"), py::arg("axis_b"), py::arg("axis_c"),
                        py::arg("density"), py::arg("maximum_degree"),
                        py::arg("maximum_order"),
                        py::arg("associated_reference_frame"),
                        py::arg("gravitational_constant") =
                            tudat::physical_constants::GRAVITATIONAL_CONSTANT,
                        R"doc(Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the density to define the mass distribution.

	Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
	The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
	Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
	The constant mass distribution is defined by the density and gravitational constant (optional).
	The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
	This function implements the models of (see Balmino [1]_).


	:param axis_a:
		Dimension of largest axis of triaxial ellipsoid.
	:param axis_b:
		Dimension of intermediate axis of triaxial ellipsoid.
	:param axis_c:
		Dimension of smallest axis of triaxial ellipsoid.
	:param density:
		Density of ellipsoid.
	:param maximum_degree:
		Maximum degree of spherical harmonics expansion.
	:param maximum_order:
		Maximum order of spherical harmonics expansion.
	:param associated_reference_frame:
		Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_constant:
		Gravitational constant G of the gravity field.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
)doc");

                    // Triaxial ellipsoid: overload 2
                    m.def(
                        "sh_triaxial_ellipsoid_from_gravitational_parameter",
                        py::overload_cast<const double, const double,
                                          const double, const int, const int,
                                          const std::string&, const double>(
                            &tss::
                                createHomogeneousTriAxialEllipsoidGravitySettings),
                        py::arg("axis_a"), py::arg("axis_b"), py::arg("axis_c"),
                        py::arg("maximum_degree"), py::arg("maximum_order"),
                        py::arg("associated_reference_frame"),
                        py::arg("gravitational_parameter"),
                        R"doc(Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the gravitational parameter to define the mass distribution..

	Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
	The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
	Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
	The constant mass distribution is defined by the gravitational parameter.
	The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
	This function implements the models of (see Balmino [1]_).


	:param axis_a:
		Dimension of largest axis of triaxial ellipsoid.
	:param axis_b:
		Dimension of intermediate axis of triaxial ellipsoid.
	:param axis_c:
		Dimension of smallest axis of triaxial ellipsoid.
	:param maximum_degree:
		Maximum degree of spherical harmonics expansion.
	:param maximum_order:
		Maximum order of spherical harmonics expansion.
	:param associated_reference_frame:
		Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_parameter:
		Gravitational parameter :math:`\mu` of gravity field.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
)doc");

                    m.def(
                        "spherical_harmonic_triaxial_body",
                        py::overload_cast<const double, const double,
                                          const double, const double, const int,
                                          const int, const std::string&,
                                          const double>(
                            &tss::
                                createHomogeneousTriAxialEllipsoidGravitySettingsDeprecated),
                        py::arg("axis_a"), py::arg("axis_b"), py::arg("axis_c"),
                        py::arg("density"), py::arg("maximum_degree"),
                        py::arg("maximum_order"),
                        py::arg("associated_reference_frame"),
                        py::arg("gravitational_constant") =
                            tudat::physical_constants::GRAVITATIONAL_CONSTANT);

                    m.def("ring_model", &tss::ringGravitySettings,
                          py::arg("gravitational_parameter"),
                          py::arg("ring_radius"),
                          py::arg("associated_reference_frame"),
                          py::arg("elliptic_integral_s_from_d_and_b"), "");
                }

            }  // namespace gravity_field
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
