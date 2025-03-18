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

namespace tudat
{

namespace simulation_setup
{

inline std::shared_ptr< GravityFieldSettings > fromFileSphericalHarmonicsGravityFieldSettings(
        const std::string& filePath,
        const int maximumDegree,
        const int maximumOrder,
        const std::string& associatedReferenceFrame = "",
        const int gravitationalParameterIndex = 0,
        const int referenceRadiusIndex = 1 )
{
    return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >(
            filePath,
            associatedReferenceFrame,
            maximumDegree,
            maximumOrder,
            gravitationalParameterIndex,
            referenceRadiusIndex );
}

inline std::shared_ptr< GravityFieldSettings > predefinedSphericalHarmonic(
        const SphericalHarmonicsModel sphericalHarmonicsModel,
        const int maximumDegree = -1 )
{
    return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >(
            sphericalHarmonicsModel, maximumDegree );
}

inline std::shared_ptr< GravityFieldSettings >
createHomogeneousTriAxialEllipsoidGravitySettingsDeprecated(
        const double axisA,
        const double axisB,
        const double axisC,
        const double ellipsoidDensity,
        const int maximumDegree,
        const int maximumOrder,
        const std::string& associatedReferenceFrame,
        const double gravitationalConstant = tpc::GRAVITATIONAL_CONSTANT )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning(
                "tudatpy.numerical_simulation.environment_setup.gravity_"
                "field.spherical_harmonic_triaxial_body",
                "tudatpy.numerical_simulation.environment_setup.gravity_"
                "field.sh_triaxial_ellipsoid_from_density" );
        isWarningPrinted = true;
    }

    return createHomogeneousTriAxialEllipsoidGravitySettings( axisA,
                                                              axisB,
                                                              axisC,
                                                              ellipsoidDensity,
                                                              maximumDegree,
                                                              maximumOrder,
                                                              associatedReferenceFrame,
                                                              gravitationalConstant );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace gravity_field
{

PYBIND11_MODULE( expose_gravity_field, m )
{
    /////////////////////////////////////////////////////////////////////////////
    // createGravityField.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_< tss::GravityFieldType >( m,
                                        "GravityFieldType",
                                        R"doc(

         Enumeration of gravity field types.

         Enumeration of gravity field types supported by tudat.





      )doc" )
            .value( "central_gravity", tss::GravityFieldType::central, R"doc(
      )doc" )
            .value( "central_spice_gravity", tss::GravityFieldType::central_spice, R"doc(
      )doc" )
            .value( "spherical_harmonic_gravity", tss::GravityFieldType::spherical_harmonic, R"doc(
      )doc" )
            .value( "polyhedron_gravity", tss::GravityFieldType::polyhedron, R"doc(
      )doc" )
            .value( "ring_gravity",
                    tss::GravityFieldType::one_dimensional_ring,
                    R"doc(No documentation found.)doc" )
            .export_values( );

    py::enum_< tss::SphericalHarmonicsModel >( m,
                                               "PredefinedSphericalHarmonicsModel",
                                               R"doc(

         Enumeration of predefined spherical harmonics models.

         Enumeration of predefined spherical harmonics models supported by tudat, for which thee coefficient files are automatically available (downloaded from
         `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_). The directory where these files are stored can be
         extracted using the :func:`~tudatpy.data.get_gravity_models_path` function.





      )doc" )
            .value( "egm96",
                    tss::SphericalHarmonicsModel::egm96,
                    R"doc(
 Coefficients for EGM96 Earth gravity field up to degree and order 200, (see `link <https://cddis.gsfc.nasa.gov/926/egm96/egm96.html>`_ )
      )doc" )
            .value( "ggm02c",
                    tss::SphericalHarmonicsModel::ggm02c,
                    R"doc(
 Coefficients for the combined GGM02 Earth gravity field up to degree and order 200, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ )
      )doc" )
            .value( "ggm02s",
                    tss::SphericalHarmonicsModel::ggm02s,
                    R"doc(
 Coefficients for the GRACE-only GGM02 Earth gravity field up to degree and order 160, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ )
      )doc" )
            .value( "goco05c",
                    tss::SphericalHarmonicsModel::goco05c,
                    R"doc(
 Coefficients for the GOCO05c combined Earth gravity field up to degree and order 719, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ )
      )doc" )
            .value( "glgm3150", tss::SphericalHarmonicsModel::glgm3150, R"doc(
 Coefficients for the GLGM3150 Moon gravity field up to degree and order 150, (see `link <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`_ )
      )doc" )
            .value( "lpe200",
                    tss::SphericalHarmonicsModel::lpe200,
                    R"doc(
 Coefficients for the LPE200 Moon gravity field up to degree and order 200, (see `link <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`_ )
      )doc" )
            .value( "gggrx1200", tss::SphericalHarmonicsModel::gggrx1200, R"doc(
 Coefficients for the GRGM1200A Moon gravity field up to degree and order 1199, (see `link <https://pgda.gsfc.nasa.gov/products/50>`_ )
      )doc" )
            .value( "jgmro120d", tss::SphericalHarmonicsModel::jgmro120d, R"doc(
 Coefficients for the MRO120D Moon gravity field up to degree and order 120, (see `link <https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/>`_ )
      )doc" )
            .value( "jgmess160a", tss::SphericalHarmonicsModel::jgmess160a, R"doc(
 Coefficients for the MESS160A Moon gravity field up to degree and order 160, (see `link <https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.lbl>`_ )
      )doc" )
            .value( "shgj180u", tss::SphericalHarmonicsModel::shgj180u, R"doc(
 Coefficients for the SHGJ180U Moon gravity field up to degree and order 180, (see `link <https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/shgj120u.lbl>`_ )
      )doc" )
            .export_values( );

    py::class_< tss::GravityFieldSettings, std::shared_ptr< tss::GravityFieldSettings > >(
            m,
            "GravityFieldSettings",
            R"doc(

         Base class for providing settings for automatic gravity field model creation.

         This class is a functional base class for settings of gravity field models that require no information in addition to their type.
         Gravity field model classes requiring additional information must be created using an object derived from this class.





      )doc" )
            //            .def(py::init<const
            //            tss::GravityFieldType>(),
            //                 py::arg("gravity_field_type"))
            .def_property_readonly( "gravity_field_type",
                                    &tss::GravityFieldSettings::getGravityFieldType,
                                    R"doc(

         **read-only**

         Type of gravity field model that is to be created.

         :type: GravityFieldType
      )doc" );

    py::class_< tss::CentralGravityFieldSettings,
                std::shared_ptr< tss::CentralGravityFieldSettings >,
                tss::GravityFieldSettings >( m,
                                             "CentralGravityFieldSettings",
                                             R"doc(

         `GravityFieldSettings` derived class defining settings of point mass gravity field.

         Derived class of `GravityFieldSettings` for central gravity fields, which are defined by a single gravitational parameter.





      )doc" )
            //            .def(py::init<double>(),
            //            py::arg("gravitational_parameter") )
            .def_property( "gravitational_parameter",
                           &tss::CentralGravityFieldSettings::getGravitationalParameter,
                           &tss::CentralGravityFieldSettings::resetGravitationalParameter,
                           R"doc(

         Gravitational parameter of central gravity field.

         :type: float
      )doc" );

    py::class_< tss::SphericalHarmonicsGravityFieldSettings,
                std::shared_ptr< tss::SphericalHarmonicsGravityFieldSettings >,
                tss::GravityFieldSettings >( m,
                                             "SphericalHarmonicsGravityFieldSettings",
                                             R"doc(

         `GravityFieldSettings` derived class defining settings of spherical harmonic gravity field representation.

         Derived class of `GravityFieldSettings` for gravity fields, which are defined by a spherical harmonic gravity field representation.





      )doc" )
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
                    &tss::SphericalHarmonicsGravityFieldSettings::getGravitationalParameter,
                    &tss::SphericalHarmonicsGravityFieldSettings::resetGravitationalParameter,
                    R"doc(

         Gravitational parameter of gravity field.

         :type: float
      )doc" )
            .def_property( "normalized_cosine_coefficients",
                           &tss::SphericalHarmonicsGravityFieldSettings::getCosineCoefficients,
                           &tss::SphericalHarmonicsGravityFieldSettings::resetCosineCoefficients,
                           R"doc(

         Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.

         :type: numpy.ndarray
      )doc" )
            .def_property( "normalized_sine_coefficients",
                           &tss::SphericalHarmonicsGravityFieldSettings::getSineCoefficients,
                           &tss::SphericalHarmonicsGravityFieldSettings::resetSineCoefficients,
                           R"doc(

         Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.

         :type: numpy.ndarray
      )doc" )
            .def_property(
                    "associated_reference_frame",
                    &tss::SphericalHarmonicsGravityFieldSettings::getAssociatedReferenceFrame,
                    &tss::SphericalHarmonicsGravityFieldSettings::resetAssociatedReferenceFrame,
                    R"doc(

         Identifier for body-fixed reference frame with which the coefficients are associated.

         :type: str
      )doc" )
            .def_property(
                    "create_time_dependent_field",
                    &tss::SphericalHarmonicsGravityFieldSettings::getCreateTimeDependentField,
                    &tss::SphericalHarmonicsGravityFieldSettings::setCreateTimeDependentField,
                    R"doc(

         Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed initially).

         :type: bool
      )doc" )
            .def_property(
                    "scaled_mean_moment_of_inertia",
                    &tss::SphericalHarmonicsGravityFieldSettings::getScaledMeanMomentOfInertia,
                    &tss::SphericalHarmonicsGravityFieldSettings::setScaledMeanMomentOfInertia,
                    R"doc(

         Value of the scaled mean moment of inertia :math:`I_{xx}+I_{yy}+I_{zz}/(MR^{2})`. This value does not influence the gravity field itself,
         but together with the degree 2 gravity field coefficients defines the body's inertia tensor.


         :type: float
      )doc" )
            .def_property_readonly(
                    "reference_radius",
                    &tss::SphericalHarmonicsGravityFieldSettings::getReferenceRadius,
                    R"doc(

         **read-only**

         Reference radius of spherical harmonic field expansion.

         :type: float
      )doc" );

    py::class_< tss::FromFileSphericalHarmonicsGravityFieldSettings,
                std::shared_ptr< tss::FromFileSphericalHarmonicsGravityFieldSettings >,
                tss::SphericalHarmonicsGravityFieldSettings >(
            m,
            "FromFileSphericalHarmonicsGravityFieldSettings",
            R"doc(No documentation found.)doc" );

    py::class_< tss::PolyhedronGravityFieldSettings,
                std::shared_ptr< tss::PolyhedronGravityFieldSettings >,
                tss::GravityFieldSettings >( m,
                                             "PolyhedronGravityFieldSettings",
                                             R"doc(

         `GravityFieldSettings` derived class defining settings of a polyhedron gravity field representation.

         Derived class of `GravityFieldSettings` for gravity fields, which are defined by a polyhedron gravity field representation.





      )doc" )
            .def_property( "gravitational_parameter",
                           &tss::PolyhedronGravityFieldSettings::getGravitationalParameter,
                           &tss::PolyhedronGravityFieldSettings::resetGravitationalParameter,
                           R"doc(

         Gravitational parameter of gravity field.

         :type: float
      )doc" )
            .def_property( "density",
                           &tss::PolyhedronGravityFieldSettings::getDensity,
                           &tss::PolyhedronGravityFieldSettings::resetDensity,
                           R"doc(

         Density of the polyhedron.

         :type: float
      )doc" )
            .def_property( "associated_reference_frame",
                           &tss::PolyhedronGravityFieldSettings::getAssociatedReferenceFrame,
                           &tss::PolyhedronGravityFieldSettings::resetAssociatedReferenceFrame,
                           R"doc(

         Identifier for body-fixed reference frame with which the vertices coordinates are associated.

         :type: str
      )doc" )
            .def_property_readonly( "vertices_coordinates",
                                    &tss::PolyhedronGravityFieldSettings::getVerticesCoordinates,
                                    R"doc(

         Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
         row per vertex, 3 columns).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly(
                    "vertices_defining_each_facet",
                    &tss::PolyhedronGravityFieldSettings::getVerticesDefiningEachFacet,
                    R"doc(

         Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
         the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
         when seen from the outside of the polyhedron.


         :type: numpy.ndarray
      )doc" );

    m.def( "central",
           &tss::centralGravitySettings,
           py::arg( "gravitational_parameter" ),
           R"doc(

 Function for central gravity field settings object.

 Function for settings object, defining a point-mass gravity field model with user-defined gravitational parameter :math:`\mu`. The gravitational potential is the defined as:

 .. math::
    U(\mathbf{r})=\frac{\mu}{||\mathbf{r}||}

 with :math:`\mathbf{r}` the position vector measured from the body's center of mass.


 Parameters
 ----------
 gravitational_parameter : float
     Gravitational parameter defining the point-mass gravity field.
 Returns
 -------
 CentralGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.CentralGravityFieldSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using a simple central gravity field model:

 .. code-block:: python

    # define parameters describing central gravity model
    gravitational_parameter = 3.986e14
    # create gravity field settings
    body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.central( gravitational_parameter )


     )doc" );

    m.def( "central_spice",
           &tss::centralGravityFromSpiceSettings,
           R"doc(

 Function to create central gravity field settings from Spice settings.

 Function for settings object, defining a point-mass gravity field model. This function provides the same model as :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.central`), but with gravitational parameter :math:`\mu` from Spice.

 Returns
 -------
 GravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` class of gravity field type ``central_spice``





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using a simple central gravity field model and data from Spice:

 .. code-block:: python

    # create gravity field settings
    body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.central_spice( )


     )doc" );

    m.def( "spherical_harmonic",
           py::overload_cast< const double,
                              const double,
                              const Eigen::MatrixXd,
                              const Eigen::MatrixXd,
                              const std::string& >( &tss::sphericalHarmonicsGravitySettings ),
           py::arg( "gravitational_parameter" ),
           py::arg( "reference_radius" ),
           py::arg( "normalized_cosine_coefficients" ),
           py::arg( "normalized_sine_coefficients" ),
           py::arg( "associated_reference_frame" ),
           R"doc(

 Function for creating a spherical harmonics gravity field settings object.

 Function for settings object, defining a gravity field model through spherical harmonic expansion.
 The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
 It represents the frame in which the spherical harmonic field is defined.

 The gravitational potential is the defined as:

 .. math::
    U(\mathbf{r})=\sum_{l=0}^{l_{max}}\sum_{m=0}^{l}\mu\left(\frac{{R}^{l}}{r^{l+1}}\right)\bar{P}_{lm}(\sin\phi)\left(\bar{C}_{lm}\cos m\theta+\bar{S}_{lm}\sin m\theta\right)

 with :math:`\mathbf{r}` the position vector of the evaluation point, measured from the body's center of mass. The angles :math:`\phi` and :math:`\theta` are the body-fixed latitude and longitude of the evaluation point, and :math:`\bar{P}_{lm}` is the associated Legendre polynomial (at degree/order :math:`l/m`).

 For the spherical harmonic gravity field (including other spherical harmonic functions), the normalized mean moment of inertia must be set by the user, to allow an inertia tensor to be computed. This is done using the :attr:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia` attribute of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class, as in the example below

 .. code-block:: python

   # Add gravity field model settings to body of spherical harmonic type
   body_settings.get( "Mars" ).gravity_field = ...

   # Add setting for moment of inertia
   body_settings.get( "Mars" ).gravity_field.scaled_mean_moment_of_inertia = 0.365

 This code snippet will automatically create a rigid body properties for Mars, with the inertia tensor computed from this value of 0.365 and the degree 2 gravity field coefficients. Note that, if gravity field variations are used for the body, time-variability of the degree 1- and 2- coefficients will be reflected in time-variability of the body's center of mass and inertia tensor.

 Note: Spherical harmonic coefficients used for this environment model must *always* be fully normalized.
 To normalize un-normalized spherical harmonic coefficients, see :func:`~tudatpy.astro.gravitation.normalize_spherical_harmonic_coefficients`.


 Parameters
 ----------
 gravitational_parameter : float
     Gravitational parameter :math:`\mu` of gravity field.
 reference_radius : float
     Reference radius :math:`R` of spherical harmonic field expansion.
 normalized_cosine_coefficients : numpy.ndarray
     Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\bar{C}_{ij}` at degree i and order j.
     As such, note that entry (0,0) of cosine coefficients should be equal to 1.

 normalized_sine_coefficients : numpy.ndarray
     Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\bar{S}_{ij}`
     at degree i and order j.

 associated_reference_frame : str
     Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
 Returns
 -------
 SphericalHarmonicsGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using a spherical harmonics gravity model:

 .. code-block:: python

   # Define the spherical harmonics gravity model
   gravitational_parameter = 3986004.415E+8
   reference_radius = 6378136.3
   # Normalized coefficients taken from https://cddis.nasa.gov/archive/egm96/general_info/egm96_to360.ascii
   # The above file is described in https://cddis.nasa.gov/archive/egm96/general_info/readme.egm96
   normalized_cosine_coefficients = [
       [1,                   0,                   0,                   0],
       [0,                   0,                   0,                   0],
       [-0.484165371736E-03, -0.186987635955E-09, 0.243914352398E-05,  0],
       [0.957254173792E-06,  0.202998882184E-05,  0.904627768605E-06,  0.721072657057E-06]
   ]
   normalized_sine_coefficients = [
       [0,                   0,                   0,                   0],
       [0,                   0,                   0,                   0],
       [0,                   0.119528012031E-08,  -0.140016683654E-05, 0],
       [0,                   0.248513158716E-06,  -0.619025944205E-06, 0.141435626958E-05]
   ]
   associated_reference_frame = "IAU_Earth"
   # Create the gravity field settings and add them to the body "Earth"
   body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.spherical_harmonic(
       gravitational_parameter,
       reference_radius,
       normalized_cosine_coefficients,
       normalized_sine_coefficients,
       associated_reference_frame )


     )doc" );

    m.def( "from_file_spherical_harmonic",
           tss::fromFileSphericalHarmonicsGravityFieldSettings,
           py::arg( "file" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           py::arg( "associated_reference_frame" ) = "",
           py::arg( "gravitational_parameter_index" ) = 0,
           py::arg( "reference_radius_index" ) = 1,
           R"doc(

 Function to load a custom spherical harmonics gravity field settings from a file.

 Function to load a custom spherical harmonics gravity field settings from a file. The file should contain **fully normalized** spherical harmonic coefficients.
 The associated gravitational parameter and reference radius should be given in m^3/s^2 and m, respectively. The file format should be the same as that used for the files
 in the directories `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_. Specifically, the file should contain

 - The first line should be a series of text blocks (typically numerical data). Two of these blocks (by default the first and second one) should be the gravitational parameter and reference radius, respectively. The text block should be separated by spaces, tabs and/or commas
 - Each subsequent line should contain a set of spherical harmonic coefficients (first ordered in ascending order by degree, then in ascending order by order), where the first, second, third and fourth value of the line should be: degree :math:`l`, order :math:`m`, normalized cosine coefficient :math:`\bar{C}_{lm}`, normalized sine coefficient :math:`\bar{S}_{lm}`. Additional entries (for instance with coefficient uncertainties) are ignored.


 Parameters
 ----------
 file : str
     Full file path and name where th gravity field file is located
 maximum_degree : int
     Maximum degree of the coefficients that are to be loaded
 maximum_order : int
     Maximum order of the coefficients that are to be loaded
 associated_reference_frame : str, default = ""
     Name of the body-fixed reference frame to which the gravity field is to be fixed. If left empty, this reference frame will automatically be set to the body-fixed frame defined by this body's rotation (see :ref:`\`\`rotation_model\`\`` for specifying rotation models).
 gravitational_parameter_index : int, default = 0
     Index of the values in the file header (first line of file) that contains the gravitational parameter
 reference_radius_index : int, default = 1
     Index of the values in the file header (first line of file) that contains the reference radius
 Returns
 -------
 SphericalHarmonicsGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for the Moon using GGGRX spherical harmonics gravity model, up to degree and order 300:

 .. code-block:: python

   # Create the gravity field settings for Moon with Spherical Harmonics loaded in from a Spherical Harmonics file
   body_settings.get("Moon").gravity_field_settings = environment_setup.gravity_field.from_file_spherical_harmonic(
       r"...\.tudat\resource\gravity_models\Moon\gggrx_1200l_sha.tab", 300, 300)


     )doc" );

    m.def( "predefined_spherical_harmonic",
           tss::predefinedSphericalHarmonic,
           py::arg( "predefined_model" ),
           py::arg( "maximum_degree" ) = -1,
           R"doc(

 Function for spherical harmonics gravity field settings of a predefined model.

 Function for spherical harmonics gravity field settings of a predefined model


 Parameters
 ----------
 predefined_model : PredefinedSphericalHarmonicsModel
     Identified for gravity field model that is to be loaded
 maximum_degree : int, default = -1
     Maximum degree and order to which the coefficients are to be loaded. If value is negative, all coefficients for the specified gravity field are loaded
 Returns
 -------
 SphericalHarmonicsGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using EGM96 spherical harmonics gravity model, up to degree and order 32:

 .. code-block:: python

   # Create the gravity field settings for Earth with Spherical Harmonics from a triaxial ellipsoid
   body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.predefined_spherical_harmonic(
       environment_setup.gravity_field.egm96, 32 )


     )doc" );

    m.def( "polyhedron_from_mu",
           py::overload_cast< const double,
                              const Eigen::MatrixXd,
                              const Eigen::MatrixXi,
                              const std::string&,
                              const double >( &tss::polyhedronGravitySettingsFromMu ),
           py::arg( "gravitational_parameter" ),
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           py::arg( "associated_reference_frame" ),
           py::arg( "gravitational_constant" ) = tpc::GRAVITATIONAL_CONSTANT,
           R"doc(

 Function for creating a polyhedron gravity field settings object, using the gravitational parameter.

 Function for settings object, defining a gravity field model through a polyhedron.
 The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
 It represents the frame in which the polyhedron field is defined.

 The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
 to Werner and Scheeres [2]_.

 This function uses the gravitational parameter to define the gravity field. To instead use the density
 constant see :func:`~tudatpy.astro.gravitation.polyhedron_from_density`. Since both models tend to be computationally intensive,
 it is recommended to use polyhedra with the lowest number of facets that allows meeting the desired accuracy. The number of facets of a polyhedron
 model can be reduced using any mesh processing software, for example `PyMeshLab <https://pymeshlab.readthedocs.io/en/latest/>`_.
 Additionally, different functions to process a polyhedron are available in `Polyhedron utilities <https://py.api.tudat.space/en/latest/polyhedron_utilities.html>`_.


 Parameters
 ----------
 gravitational_parameter : float
     Gravitational parameter :math:`\mu` of gravity field.
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 associated_reference_frame : str
     Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
 gravitational_constant : float, default=GRAVITATIONAL_CONSTANT
     Newton's gravitational constant G, used to computed the density

 Returns
 -------
 PolyhedronGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class






     )doc" );

    m.def( "polyhedron_from_density",
           py::overload_cast< const double,
                              const Eigen::MatrixXd,
                              const Eigen::MatrixXi,
                              const std::string&,
                              const double >( &tss::polyhedronGravitySettings ),
           py::arg( "density" ),
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           py::arg( "associated_reference_frame" ),
           py::arg( "gravitational_constant" ) = tpc::GRAVITATIONAL_CONSTANT,
           R"doc(

 Function for creating a polyhedron gravity field settings object, using the density.

 Function for settings object, defining a gravity field model through a polyhedron.
 The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
 It represents the frame in which the polyhedron field is defined.

 The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
 to Werner and Scheeres [2]_.

 This function uses the density to define the gravity field. To instead use the
 gravitational parameter see :func:`~tudatpy.astro.gravitation.polyhedron_from_mu`.


 Parameters
 ----------
 density : float, default=TUDAT_NAN
     Density of the polyhedron.

 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 associated_reference_frame : str
     Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
 gravitational_constant : float, default=GRAVITATIONAL_CONSTANT
     Newton's gravitational constant G, used to computed the gravitational parameter

 Returns
 -------
 PolyhedronGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class






     )doc" );

    // Triaxial ellipsoid: overload 1
    m.def( "sh_triaxial_ellipsoid_from_density",
           py::overload_cast< const double,
                              const double,
                              const double,
                              const double,
                              const int,
                              const int,
                              const std::string&,
                              const double >(
                   &tss::createHomogeneousTriAxialEllipsoidGravitySettings ),
           py::arg( "axis_a" ),
           py::arg( "axis_b" ),
           py::arg( "axis_c" ),
           py::arg( "density" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           py::arg( "associated_reference_frame" ),
           py::arg( "gravitational_constant" ) = tudat::physical_constants::GRAVITATIONAL_CONSTANT,
           R"doc(

 Function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the density to define the mass distribution.

 Function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
 The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
 Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
 The constant mass distribution is defined by the density and gravitational constant (optional).
 The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
 This function implements the models of (see Balmino [1]_).


 Parameters
 ----------
 axis_a : float
     Dimension of largest axis of triaxial ellipsoid.
 axis_b : float
     Dimension of intermediate axis of triaxial ellipsoid.
 axis_c : float
     Dimension of smallest axis of triaxial ellipsoid.
 density : float
     Density of ellipsoid.
 maximum_degree : int
     Maximum degree of spherical harmonics expansion.
 maximum_order : int
     Maximum order of spherical harmonics expansion.
 associated_reference_frame : str
     Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
 gravitational_constant : float, default=physical_constants::GRAVITATIONAL_CONSTANT
     Gravitational constant G of the gravity field.
 Returns
 -------
 SphericalHarmonicsGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





 Examples
 --------
 In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using the expansion of a homogeneous triaxial ellipsoid into a spherical harmonics gravity model:

 .. code-block:: python

   # Create the gravity field settings for Earth with Spherical Harmonics from a triaxial ellipsoid
   body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.spherical_harmonic_triaxial_ellipsoid_from_density(
       axis_a=6378171.88,
       axis_b=6378102.03,
       axis_c=6356752.24,
       density=5520,
       maximum_degree=50,
       maximum_order=50,
       associated_reference_frame="IAU_Earth" )


     )doc" );

    // Triaxial ellipsoid: overload 2
    m.def( "sh_triaxial_ellipsoid_from_gravitational_parameter",
           py::overload_cast< const double,
                              const double,
                              const double,
                              const int,
                              const int,
                              const std::string&,
                              const double >(
                   &tss::createHomogeneousTriAxialEllipsoidGravitySettings ),
           py::arg( "axis_a" ),
           py::arg( "axis_b" ),
           py::arg( "axis_c" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           py::arg( "associated_reference_frame" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the gravitational parameter to define the mass distribution..

 Function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
 The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
 Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
 The constant mass distribution is defined by the gravitational parameter.
 The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
 This function implements the models of (see Balmino [1]_).


 Parameters
 ----------
 axis_a : float
     Dimension of largest axis of triaxial ellipsoid.
 axis_b : float
     Dimension of intermediate axis of triaxial ellipsoid.
 axis_c : float
     Dimension of smallest axis of triaxial ellipsoid.
 maximum_degree : int
     Maximum degree of spherical harmonics expansion.
 maximum_order : int
     Maximum order of spherical harmonics expansion.
 associated_reference_frame : str
     Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
 gravitational_parameter : float
     Gravitational parameter :math:`\mu` of gravity field.
 Returns
 -------
 SphericalHarmonicsGravityFieldSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class






     )doc" );

    m.def( "spherical_harmonic_triaxial_body",
           py::overload_cast< const double,
                              const double,
                              const double,
                              const double,
                              const int,
                              const int,
                              const std::string&,
                              const double >(
                   &tss::createHomogeneousTriAxialEllipsoidGravitySettingsDeprecated ),
           py::arg( "axis_a" ),
           py::arg( "axis_b" ),
           py::arg( "axis_c" ),
           py::arg( "density" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           py::arg( "associated_reference_frame" ),
           py::arg( "gravitational_constant" ) =
                   tudat::physical_constants::GRAVITATIONAL_CONSTANT );

    m.def( "ring_model",
           &tss::ringGravitySettings,
           py::arg( "gravitational_parameter" ),
           py::arg( "ring_radius" ),
           py::arg( "associated_reference_frame" ),
           py::arg( "elliptic_integral_s_from_d_and_b" ),
           R"doc(No documentation found.)doc" );
}

}  // namespace gravity_field
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
