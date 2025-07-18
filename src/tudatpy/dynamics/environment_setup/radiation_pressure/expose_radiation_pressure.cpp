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
#include "expose_radiation_pressure.h"

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>

// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace radiation_pressure
{

void expose_radiation_pressure_setup( py::module& m )
{
    /////////////////////////////////////////////////////////////////////////////
    // createRadiationPressureInterface.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_< tss::RadiationPressureType >(
            m, "RadiationPressureType", R"doc(No documentation found.)doc" )
            .value( "cannonball_radiation_pressure_interface",
                    tss::RadiationPressureType::cannon_ball_radiation_pressure_interface,
                    R"doc(No documentation found.)doc" )
            //                .value("panelled_radiation_pressure_interface",
            //                       tss::RadiationPressureType::panelled_radiation_pressure_interface,
            //                       get_docstring("RadiationPressureType.panelled_radiation_pressure_interface").c_str())
            //                .value("solar_sailing_radiation_pressure_interface",
            //                       tss::RadiationPressureType::solar_sailing_radiation_pressure_interface,
            //                       get_docstring("RadiationPressureType.solar_sailing_radiation_pressure_interface").c_str())
            .export_values( );
    py::enum_< tss::RadiationPressureTargetModelType >(
            m, "RadiationPressureTargetModelType", R"doc(No documentation found.)doc" )
            .value( "cannonball_target",
                    tss::RadiationPressureTargetModelType::cannonball_target,
                    R"doc(No documentation found.)doc" )
            .value( "paneled_target",
                    tss::RadiationPressureTargetModelType::paneled_target,
                    R"doc(No documentation found.)doc" )
            .value( "multi_type_target",
                    tss::RadiationPressureTargetModelType::multi_type_target,
                    R"doc(No documentation found.)doc" )
            .value( "undefined_target",
                    tss::RadiationPressureTargetModelType::undefined_target,
                    R"doc(No documentation found.)doc" )
            .export_values( );

    py::class_< tss::RadiationPressureInterfaceSettings,
                std::shared_ptr< tss::RadiationPressureInterfaceSettings > >(
            m, "RadiationPressureInterfaceSettings", R"doc(No documentation found.)doc" );
    //            .def(py::init<const
    //            tss::RadiationPressureType, const std::string
    //            &,
    //                 const std::vector<std::string>>(),
    //                 py::arg("radiation_pressure_type"),
    //                 py::arg("source_body"),
    //                 py::arg("occulting_bodies") =
    //                 std::vector<std::string>());

    py::class_< tss::CannonBallRadiationPressureInterfaceSettings,
                std::shared_ptr< tss::CannonBallRadiationPressureInterfaceSettings >,
                tss::RadiationPressureInterfaceSettings >(
            m, "CannonBallRadiationPressureInterfaceSettings", R"doc(No documentation found.)doc" );

    //            .def(py::init<const std::string &, const
    //            double, const double,
    //                 const std::vector<std::string> &>(),
    //                 py::arg("source_body"), py::arg("area"),
    //                 py::arg("radiation_pressure_coefficient"),
    //                 py::arg("occulting_bodies") =
    //                 std::vector<std::string>())

    m.def( "cannonball",
           py::overload_cast< const std::string&,
                              const double,
                              const double,
                              const std::vector< std::string >& >(
                   &tss::cannonBallRadiationPressureSettings ),
           py::arg( "source_body" ),
           py::arg( "reference_area" ),
           py::arg( "radiation_pressure_coefficient" ),
           py::arg( "occulting_bodies" ) = std::vector< std::string >( ),
           R"doc(No documentation found.)doc" );

    //        m.def("panelled",
    //              &tss::panelledRadiationPressureInterfaceSettings,
    //              py::arg("source_body"),
    //              py::arg("emissivities"),
    //              py::arg("areas"),
    //              py::arg("diffusion_coefficients"),
    //              py::arg("surface_normals_in_body_fixed_frame"),
    //              py::arg("occulting_bodies") =
    //              std::vector<std::string>(),
    //              get_docstring("panelled").c_str()
    //              );

    ///////////////////////////////////////////////////////////
    ///////////   ENUMS
    ///////////////////////////////////////////////////////////

    py::enum_< tss::KnockeTypeSurfacePropertyDistributionModel >(
            m, "KnockeTypeSurfacePropertyDistributionModel" )
            .value( "custom", tss::KnockeTypeSurfacePropertyDistributionModel::custom )
            .value( "albedo_knocke",
                    tss::KnockeTypeSurfacePropertyDistributionModel::albedo_knocke )
            .value( "emissivity_knocke",
                    tss::KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke )
            .export_values( );

    py::enum_< tss::SphericalHarmonicsSurfacePropertyDistributionModel >(
            m, "SphericalHarmonicsSurfacePropertyDistributionModel" )
            .value( "albedo_dlam1",
                    tss::SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1 )
            .export_values( );

    enum class SphericalHarmonicsSurfacePropertyDistributionModel {
        custom,
        albedo_dlam1 /**< DLAM-1 lunar albedo model:
                        Floberghagen, R. et al. "Lunar
                        Albedo Force Modeling and its Effect
                        on Low Lunar Orbit and Gravity Field
                        Determination". ASR 23. 4(1999):
                        733-738. */
    };

    ///////////////////////////////////////////////////////////
    ///////////   LUMINOSITY MODELS
    ///////////////////////////////////////////////////////////

    py::class_< tss::LuminosityModelSettings, std::shared_ptr< tss::LuminosityModelSettings > >(
            m,
            "LuminosityModelSettings",
            R"doc(

         Base class for providing settings for body source luminosity settings, to be used (typically but not necessarily) for defining the Sun's luminosity.





      )doc" );

    m.def( "constant_luminosity",
           &tss::constantLuminosityModelSettings,
           py::arg( "luminosity" ),
           R"doc(

 Function for creating constant radiation source luminosity settings.

 Function for creating constant radiation source luminosity settings, defining the total
 radiated power (in Watts) of a given source. With this function, the source luminosity is constant,
 and is assumed to emit radiation isotropically.


 Parameters
 ----------
 luminosity : float
     Constant source luminosity (in Watt)
 Returns
 -------
 LuminosityModelSettings
     Object defining settings for source luminosity






     )doc" );

    m.def( "irradiance_based_constant_luminosity",
           &tss::irradianceBasedLuminosityModelSettings,
           py::arg( "constant_irradiance" ),
           py::arg( "reference_distance" ),
           R"doc(

 Function for creating source luminosity settings based on the irradiance at a reference distance.

 Function for creating source luminosity based on the irradiance at a reference distance. For instance,
 one can provide the solar irradiance at 1 AU, and this will be translated to the Sun's luminosity. With this function,
 the source luminosity is constant, and is assumed to emit radiation isotropically.


 Parameters
 ----------
 constant_irradiance : float
     Irradiance at reference distance from center of source (in :math:`W/m^{2}`)
 reference_distance : float
     Distance from center of source at which the irradiance is defined
 Returns
 -------
 LuminosityModelSettings
     Object defining settings for source luminosity

 Examples
 --------
 In this example, we create a constant luminosity model for the Sun, based on the solar constant, which defines the irradiance at 1 AU. Different values exist in literature for the solar constant, in this case we assume 1367 W/m^2 at 1 AU (Vallado 2013).
 Assuming an isotropic radiation source, we then create the radiation source settings for the Sun.

 .. code-block:: python

     ...

     from tudatpy import constants

     irradiance_at_1AU = 1367.0  # W/m^2, Vallado 2013

     luminosity_model_settings = (
         environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
             irradiance_at_1AU, constants.ASTRONOMICAL_UNIT
         )
     )
     radiation_source_settings_sun = (
         environment_setup.radiation_pressure.isotropic_radiation_source(
             luminosity_model_settings
         )
     )

     body_settings.get("Sun").radiation_source_settings = radiation_source_settings_sun




     )doc" );

    m.def( "time_variable_luminosity",
           &tss::timeVariableLuminosityModelSettings,
           py::arg( "luminosity_function" ),
           R"doc(

 Function for creating time-variable radiation source luminosity settings.

 Function for creating time-variable radiation source luminosity settings, defining the total
 radiated power (in Watts) of a given source as a function of time. With this function, the source
 is assumed to emit radiation isotropically.


 Parameters
 ----------
 luminosity_function : Callable[[float], float]
     Function returning source luminosity (in Watt) as a function of time
 Returns
 -------
 LuminosityModelSettings
     Object defining settings for source luminosity






     )doc" );

    m.def( "irradiance_based_time_variable_luminosity",
           &tss::timeVariableIrradianceBasedLuminosityModelSettings,
           py::arg( "irradiance_function" ),
           py::arg( "reference_distance" ),
           R"doc(

 Function for creating time-variable source luminosity settings based on the irradiance at a reference distance.

 Function for creating source time-variable luminosity based on the irradiance at a reference distance. For instance,
 one can provide the solar irradiance at 1 AU as a function of time, and this will be translated to the Sun's luminosity.
 With this function, the source is assumed to emit radiation isotropically.


 Parameters
 ----------
 irradiance_function : Callable[[float], float]
     Function returning irradiance at reference distance from center of source (in :math:`W/m^{2}`) as a function fo time
 reference_distance : float
     Distance from center of source at which the irradiance is defined
 Returns
 -------
 LuminosityModelSettings
     Object defining settings for source luminosity






     )doc" );

    ///////////////////////////////////////////////////////////
    ///////////   SURFACE PROPERTY MODELS
    ///////////////////////////////////////////////////////////

    py::class_< tss::SurfacePropertyDistributionSettings,
                std::shared_ptr< tss::SurfacePropertyDistributionSettings > >(
            m,
            "SurfacePropertyDistributionSettings",
            R"doc(

         Base class for providing settings for body surface property distribution settings, to be used (typically but not necessarily) for defining surface distribution of albedo and emissivity of solar system bodies for calculations of albedo and planetary radiation pressure.Note that not all albedo/emissivity models require this type of distribution model





      )doc" );

    m.def( "constant_surface_property_distribution",
           &tss::constantSurfacePropertyDistributionSettings,
           py::arg( "constant_value" ),
           R"doc(

 Function for creating constant radiative surface property distribution settings.

 Function for creating constant radiative surface property (e.g. albedo, emissivity, etc.) distribution settings.


 Parameters
 ----------
 constant_value : float
     Constant surface property value
 Returns
 -------
 SurfacePropertyDistributionSettings
     Object defining settings for surface property distribution






     )doc" );

    m.def( "spherical_harmonic_surface_property_distribution",
           py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd& >(
                   &tss::sphericalHarmonicsSurfacePropertyDistributionSettings ),
           py::arg( "cosine_coefficients" ),
           py::arg( "sine_coefficients" ),
           R"doc(

 Function for creating radiative surface property distribution settings according to a spherical harmonic model.

 Function for creating radiative surface property (e.g. albedo, emissivity, etc.) distribution settings
 according to a spherical harmonic model. The user provides unnormalized cosine and sine coefficients :math:`C_{lm}` and :math:`S_{lm}`,
 from which the surface property :math:`k` is computed from:

 .. math::
    k(\phi,\theta)=\sum_{l=0}^{l_{max}}\sum_{m=0}^{l}{P}_{lm}(\sin\phi)\left({C}_{lm}\cos m\theta+{S}_{lm}\sin m\theta\right)

 with the angles :math:`\phi` and :math:`\theta` the body-fixed latitude and longitude of the evaluation point.


 Parameters
 ----------
 cosine_coefficients : numpy.ndarray
     Cosine coefficients of surface distribution. Entry (i,j) denotes coefficient :math:`{C}_{ij}` at degree i and order j.
 sine_coefficients : numpy.ndarray
     Sine coefficients of surface distribution. Entry (i,j) denotes coefficient :math:`{C}_{ij}` at degree i and order j.
 Returns
 -------
 SurfacePropertyDistributionSettings
     Object defining settings for surface property distribution






     )doc" );

    m.def( "predefined_spherical_harmonic_surface_property_"
           "distribution",
           py::overload_cast< tss::SphericalHarmonicsSurfacePropertyDistributionModel >(
                   &tss::sphericalHarmonicsSurfacePropertyDistributionSettings ),
           py::arg( "predefined_model" ),
           R"doc(

 Function for creating radiative surface property distribution settings according to a predefined spherical harmonic model.

 As :func:`spherical_harmonic_surface_property_distribution`, but with a predefined spherical harmonic distribution.


 Parameters
 ----------
 predefined_model : SphericalHarmonicsSurfacePropertyDistributionModel
     Identifier for predefined spherical harmonic surface property model.
 Returns
 -------
 SurfacePropertyDistributionSettings
     Object defining settings for surface property distribution






     )doc" );

    m.def( "knocke_type_surface_property_distribution",
           &tss::manualSecondDegreeZonalPeriodicSurfacePropertyDistributionSettings,
           py::arg( "constant_contribution" ),
           py::arg( "constant_degree_one_contribution" ),
           py::arg( "cosine_periodic_degree_one_contribution" ),
           py::arg( "sine_periodic_degree_one_contribution" ),
           py::arg( "constant_degree_two_contribution" ),
           py::arg( "reference_epoch" ),
           py::arg( "period" ),
           R"doc(

 Function for creating radiative surface property distribution settings according to 'Knocke-type' model

 Function for creating radiative surface property (e.g. albedo, emissivity, etc.) distribution settings
 according to a model such as the one used by Knocke (1988). This model uses a degree two zonal spherical harmonic model, with
 a sinusoidal variation in the degree one coefficient. The surface property :math:`k` is computed from:

 .. math::
    k(\phi,t)=a_{0}+a_{1}P_{1}(\sin\phi)+a_{2}P_{2}(\sin\phi)
 .. math::
    a_{1}=c_{0}+c_{1}\cos\left(\frac{2\pi(t-t_{0})}{T}\right)+c_{2}\sin\left(\frac{2\pi(t-t_{0})}{T}\right)

 with the angle :math:`\phi` denotes the body-fixed latitude of the evaluation point, and :math:`t`, :math:`t_{0}` and :math:`T` define the current time,
 reference time and period of the variation, respectively. The coefficients :math:`a_{0}, a_{2}, c_{0}, c_{1}, c_{2}` are provided by the user.


 Parameters
 ----------
 constant_contribution : float
     Value of :math:`a_{0}` in above formulation.
 constant_degree_one_contribution : float
     Value of :math:`c_{0}` in above formulation.
 cosine_periodic_degree_one_contribution : float
     Value of :math:`c_{1}` in above formulation.
 sine_periodic_degree_one_contribution : float
     Value of :math:`c_{2}` in above formulation.
 constant_degree_two_contribution : float
     Value of :math:`a_{2}` in above formulation.
 reference_epoch : float
     Reference epoch :math:`t_{0}` of the periodic variation.
 period : float
     Period :math:`T` of the periodic variation.
 Returns
 -------
 SurfacePropertyDistributionSettings
     Object defining settings for surface property distribution






     )doc" );

    m.def( "predefined_knocke_type_surface_property_distribution",
           py::overload_cast< tss::KnockeTypeSurfacePropertyDistributionModel >(
                   &tss::secondDegreeZonalPeriodicSurfacePropertyDistributionSettings ),
           py::arg( "predefined_model" ),
           R"doc(

 Function for creating radiative surface property distribution settings according to a predefined 'Knocke-type` model.

 As :func:`spherical_harmonic_surface_property_distribution`, but with a predefined spherical harmonic distribution.


 Parameters
 ----------
 predefined_model : KnockeTypeSurfacePropertyDistributionModel
     Identifier for predefined Knocke-type surface property model.
 Returns
 -------
 SurfacePropertyDistributionSettings
     Object defining settings for surface property distribution






     )doc" );

    m.def( "custom_surface_property_distribution",
           &tss::customSurfacePropertyDistributionSettings,
           py::arg( "custom_function" ),
           R"doc(

 Function for creating radiative surface property distribution settings according to a custom user-defined model.

 Function for creating radiative surface property (e.g. albedo, emissivity, etc.) distribution settings
 according to a custom user-defined model, as a function of latitude, longitude and time.


 Parameters
 ----------
 custom_function : Callable[[float, float, float], float]
     Function providing surface property as a function of latitude, longitude and time (in that order).
 Returns
 -------
 SurfacePropertyDistributionSettings
     Object defining settings for surface property distribution






     )doc" );

    ///////////////////////////////////////////////////////////
    ///////////   PANEL RADIOSITY MODELS
    ///////////////////////////////////////////////////////////

    py::class_< tss::PanelRadiosityModelSettings,
                std::shared_ptr< tss::PanelRadiosityModelSettings > >(
            m,
            "PanelRadiosityModelSettings",
            R"doc(

         Base class for providing settings for body panel radiosity models, to be used (typically but not necessarily) for defining surface radiosity of a panelled solar system body as a result of albedo and/or planetary radiation pressure





      )doc" );

    m.def( "constant_radiosity",
           &tss::constantPanelRadiosityModelSettings,
           py::arg( "radiosity" ),
           R"doc(

 Function for creating settings for surface constant surface radiosity of an extended source

 Function for creating settings for surface radiosity of an extended source, using constant Lambertian radiosity :math:`J` (in :math:`W/m^{2}`).
 For a surface panel normal of :math:`\hat{\mathbf{n}}` and a vector :math:`\mathbf{r}` from the surface element to the target, the resulting
 irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if :math:`\theta>0`, or in other words if the panel is visible from the target):

 .. math::
    \Phi=J\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

 with :math:`A` the panel area, :math:`\theta` is the angle between :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


 Parameters
 ----------
 radiosity : float
     Constant Lambertian radiosity from surface in :math:`W/m^{2}`.
 Returns
 -------
 PanelRadiosityModelSettings
     Object defining settings for source panel radiosity






     )doc" );

    m.def( "constant_albedo_surface_radiosity",
           py::overload_cast< double, const std::string& >(
                   &tss::albedoPanelRadiosityModelSettings ),
           py::arg( "constant_albedo" ),
           py::arg( "original_source_name" ),
           R"doc(

 Function for creating settings for surface constant albedo surface radiosity of an extended source

 Function for creating settings for surface radiosity of an extended source, with surface radiation the result
 of albedo using a Lambertian scattering law, and a constant albedo value over the surface.
 For a surface panel normal of :math:`\hat{\mathbf{n}}`, a vector :math:`\mathbf{r}` from the surface element to the target, and a
 vector :math:`\mathbf{r}_{s}` from the surface element to the original source (typically the Sun),
 the resulting irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):

 .. math::
    \Phi=\cos\theta_{s}\Phi_{s}\frac{a}{\pi}\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

 with :math:`\theta_{s}` the angle between :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r_{s}}`, :math:`\Phi_{s}` the irradiance from the original source
 at the panel of the reflecting body, :math:`a` is the albedo coefficient, :math:`A` the panel area, :math:`\theta` is the angle between :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


 Parameters
 ----------
 constant_albedo : float
     Constant value of the albedo coefficient :math:`a`.
 original_source_name : str
     Name of the original source from which the radiation is reflection to the target.
 Returns
 -------
 PanelRadiosityModelSettings
     Object defining settings for source panel radiosity






     )doc" );

    m.def( "variable_albedo_surface_radiosity",
           &tss::albedoPanelRadiosityModelSettingsGeneric,
           py::arg( "albedo_distribution_settings" ),
           py::arg( "original_source_name" ),
           R"doc(

 Function for creating settings for surface variable albedo surface radiosity of an extended source

 As :func:`constant_albedo_surface_radiosity`, but with the surface albedo :math:`a` defined by a surface distribution model.


 Parameters
 ----------
 albedo_distribution_settings : SurfacePropertyDistributionSettings
     Model for the surface distribution of the albedo :math:`a`.
 original_source_name : str
     Name of the original source from which the radiation is reflection to the target.
 Returns
 -------
 PanelRadiosityModelSettings
     Object defining settings for source panel radiosity






     )doc" );

    m.def( "thermal_emission_blackbody_constant_emissivity",
           py::overload_cast< double, const std::string& >(
                   &tss::delayedThermalPanelRadiosityModelSettings ),
           py::arg( "constant_emissivity" ),
           py::arg( "original_source_name" ),
           R"doc(

 Function for creating settings for surface radiosity of an extended source from an isotropically heated body with constant emissivity

 Function for creating settings for surface radiosity of an extended source from an isotropically heated body (e.g. IR radiation) with constant surface
 emissivity,
 where the emitted power of the body is computed from the assumption that all heat absorbed from an original source is
 emitted isotropically by the body. For instance, for Earth with Sun as original source, this model is equivalent to
 assuming that a given fraction of all heat incident of the Sun on the Earth is absorbed and causes the full Earth surface to
 heat to a constant temperature, which then results in the body emitting infrared radiation from its surface.

 For a surface panel normal of :math:`\hat{\mathbf{n}}`, a vector :math:`\mathbf{r}` from the surface element to the target,
 the resulting irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):

 .. math::
    \Phi=\frac{\epsilon\Phi_{s}}{4}\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

 with :math:`\epsilon` the emissivity, :math:`\Phi_{s}` the irradiance from the original source,  :math:`A` the panel area, :math:`\theta` is the angle between
 :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


 Parameters
 ----------
 constant_emissivity : float
     Constant emissivity of the surface :math:`\epsilon`.
 original_source_name : str
     Name of the original source from which the radiation is reflection to the target.
 Returns
 -------
 PanelRadiosityModelSettings
     Object defining settings for source panel radiosity






     )doc" );

    m.def( "thermal_emission_blackbody_variable_emissivity",
           &tss::delayedThermalPanelRadiosityModelSettingsGeneric,
           py::arg( "emissivity_distribution_model" ),
           py::arg( "original_source_name" ),
           R"doc(

 Function for creating settings for surface radiosity of an extended source from an isotropically heated body with variable emissivity

 As :func:`thermal_emission_blackbody_constant_emissivity`, but with the surface emissivity :math:`\epsilon` defined by a surface distribution model.


 Parameters
 ----------
 emissivity_distribution_model : SurfacePropertyDistributionSettings
     Model for the surface distribution of the emissivity :math:`\epsilon`.
 original_source_name : str
     Name of the original source from which the radiation is reflection to the target.
 Returns
 -------
 PanelRadiosityModelSettings
     Object defining settings for source panel radiosity






     )doc" );

    m.def( "thermal_emission_angle_based_radiosity",
           &tss::angleBasedThermalPanelRadiosityModelSettings,
           py::arg( "minimum_temperature" ),
           py::arg( "maximum_temperature" ),
           py::arg( "constant_emissivity" ),
           py::arg( "original_source_name" ),
           R"doc(

 Function for creating settings for surface radiosity of an extended source with surface temperature from Lemoine (2013)

 Function for creating settings for surface radiosity of an extended source from an isotropically heated body (e.g. IR radiation)
 with surface temperature :math:`T` computed from the angle of the surface normal and the original source as follows:

 .. math::
    T=\max\left(T_{max}(\cos\phi_{s})^{1/4},T_{min} \right)

 with :math:`phi_{s}` the angle along a great circle arc from the panel to the subsolar (for the Sun as original source) point; for
 a circular body equivalent to the angle of the vector to the original source and the surface normal. The minimum and
 maximum temperatures are user parameters.

 For a surface panel normal of :math:`\hat{\mathbf{n}}`, a vector :math:`\mathbf{r}` from the surface element to the target,
 the resulting irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):

 .. math::
    \Phi=\epsilon kT^{4}\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

 with :math:`\epsilon` the emissivity, :math:`k` the Stefan-Boltzmann constant, :math:`A` the panel area, :math:`\theta` is the angle between
 :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


 Parameters
 ----------
 minimum_temperature : float
     Minimum surface temperature :math:`T_{min}`.
 maximum_temperature : float
     Maximum surface temperature :math:`T_{min}`.
 constant_emissivity : float
     Constant emissivity of the surface :math:`\epsilon`.
 original_source_name : str
     Name of the original source from which the radiation is reflection to the target.
 Returns
 -------
 PanelRadiosityModelSettings
     Object defining settings for source panel radiosity






     )doc" );

    ///////////////////////////////////////////////////////////
    ///////////   PANEL REFLECTION MODELS
    ///////////////////////////////////////////////////////////

    py::class_< tss::BodyPanelReflectionLawSettings,
                std::shared_ptr< tss::BodyPanelReflectionLawSettings > >(
            m,
            "BodyPanelReflectionLawSettings",
            R"doc(

         Base class for providing settings for body panel reflection law models, to be used for defining spacecraft surface properties relevant for the computation of radiation pressure acting on a macromodel.





      )doc" );

    m.def( "specular_diffuse_body_panel_reflection",
           &tss::specularDiffuseBodyPanelReflectionLawSettings,
           py::arg( "specular_reflectivity" ),
           py::arg( "diffuse_reflectivity" ),
           py::arg( "with_instantaneous_reradiation" ),
           R"doc(

 Function for creating settings for target panel reflection law using a specular-diffuse model

 Function for creating settings for target panel reflection law used for a radiation pressure target, with a
 specular diffuse model. The details of the implementation are given by Montenbruck et al. (2015). The reflection
 law is defined by the absorption coefficient :math:`\alpha`, diffuse reflectivity :math:`\delta` and specular reflectivity
 :math:`\rho`, which must meet the condition :math:`\alpha+\delta+\rho=1`. For the model definition, the user provides
 :math:`\alpha` and :math:`\delta` (and :math:`\rho` is calculated). The reaction vector :math:`\hat{\mathbf{f}}` for a panel
 with surface normal :math:`\hat{\mathbf{n}}`, and unit vector from panel surface to source :math:`\hat{\mathbf{r}}` then becomes:

 .. math::
    \hat{\mathbf{f}}=\cos\theta\left((\alpha+\delta)\hat{\mathbf{r}}+(\frac{2}{3}\delta+2\rho\cos\theta)\hat{\mathbf{n}} \right)

 In addition, it can be specified whether the absorbed radiation is also instantaneously retransmitted (according to Lambert's law), in
 which case the above is modified to become:

 .. math::
    \hat{\mathbf{f}}=\cos\theta\left((\alpha+\delta)\left(\hat{\mathbf{r}}+\frac{2}{3}\hat{\mathbf{n}}\right)+2\rho\cos\theta\hat{\mathbf{n}} \right)


 Parameters
 ----------
 specular_reflectivity : float
     Specular reflectivity :math:`\rho`.
 diffuse_reflectivity : float
     Diffuse reflectivity :math:`\delta`.
 with_instantaneous_reradiation : bool
     Boolean denoting whether absorbed radiation is instantaneously retransmitted (yes, if true).
 Returns
 -------
 BodyPanelReflectionLawSettings
     Object defining settings for target panel reflection law






     )doc" );

    m.def( "lambertian_body_panel_reflection",
           &tss::lambertainBodyPanelReflectionLawSettings,
           py::arg( "reflectivity" ),
           R"doc(

 Function for creating settings for target panel reflection law using a Lambertian model

 Function for creating settings for target panel reflection law used for a radiation pressure target, with a
 purely Lambertian model. The implementation is as :func:`specular_diffuse_body_panel_reflection`, with
 :math:`\rho=0` and no instantaneous reradiation. The only free parameter is the reflectivity :math:`\delta`, such that
 :math:`\alpha=1-\delta`.


 Parameters
 ----------
 reflectivity : float
     Reflectivity :math:`\delta`
 Returns
 -------
 BodyPanelReflectionLawSettings
     Object defining settings for target panel reflection law






     )doc" );

    ///////////////////////////////////////////////////////////
    ///////////   RADIATION SOURCE MODELS
    ///////////////////////////////////////////////////////////

    py::class_< tss::RadiationSourceModelSettings,
                std::shared_ptr< tss::RadiationSourceModelSettings > >(
            m,
            "RadiationSourceModelSettings",
            R"doc(

         Base class for providing settings for properties of a radiation source (e.g. Sun), to be used in the context of (for instance) calculation of radiation pressure on spacecraft





      )doc" );

    m.def( "isotropic_radiation_source",
           &tss::isotropicPointRadiationSourceModelSettings,
           py::arg( "luminosity_model" ),
           R"doc(

 Function for creating settings for an isotropic radiation source

 Function for creating settings for a radiation source that emits isotropically. The source is provided
 with a luminosity model :math:`L(t)` as a (possible) function of time :math:`t`. The irradiance :math:`\Phi` at a relative position
 :math:`\mathbf{r}` from the source's center is then computed from:

 .. math::
    \Phi=\frac{L}{4\pi||\mathbf{r}||^{2}}


 Parameters
 ----------
 luminosity_model : LuminosityModelSettings
     Settings for the luminosity model.
 Returns
 -------
 RadiationSourceModelSettings
     Object defining settings for source model irradiance
 Examples
 --------
 In this example, we create a constant luminosity model for the Sun, based on the solar constant, which defines the irradiance at 1 AU. Different values exist in literature for the solar constant, in this case we assume 1367 W/m^2 at 1 AU (Vallado 2013).
 Assuming an isotropic radiation source, we then create the radiation source settings for the Sun.

 .. code-block:: python

     ...

     from tudatpy import constants

     irradiance_at_1AU = 1367.0  # W/m^2, Vallado 2013

     luminosity_model_settings = (
         environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
             irradiance_at_1AU, constants.ASTRONOMICAL_UNIT
         )
     )
     radiation_source_settings_sun = (
         environment_setup.radiation_pressure.isotropic_radiation_source(
             luminosity_model_settings
         )
     )

     body_settings.get("Sun").radiation_source_settings = radiation_source_settings_sun






     )doc" );

    m.def( "panelled_extended_radiation_source",
           &tss::extendedRadiationSourceModelSettingsWithOccultationMap,
           py::arg( "panel_radiosity_settings" ),
           py::arg( "number_of_panels_per_ring" ),
           py::arg( "original_source_occulting_bodies" ) =
                   std::map< std::string, std::vector< std::string > >( ),
           R"doc(

 Function for creating settings for a dynamically panelled extended radiation source

 Function for creating settings for a radiation source that is modelled as an anisotropic extended source,
 such as a source due to albedo or planetary IR. The model can combined any number of superimposed surface panel radiosity models
 (e.g. albedo, direct radiation), each of which may or may not involve an 'original source' (e.g. the Sun for albedo).
 Each time the radiation at a given target is computed, the surface of the body is re-panelled, using the algorithm described by
 Knocke (1989). In short, a single panel is placed at the zenith of the evaluation point, with any number of rings around it, each of
 which has any number of (equispaced) panels on it. The width of each ring is such that all panels have the same projected, attenuated area.
 The panelling settings are defined by the user to this function. The
 The irradiance :math:`\Phi` at a relative position :math:`\mathbf{r}` from the source's center is then computed from all
 :math:`N` panels, each of which evaluated :math:`M` panel radiosity models

 .. math::
    \Phi=\sum_{i=1}^{N}\sum_{j=1}\Phi_{i,j}

 where :math:`\Phi_{i,j}` denotes the contribution to the total irradiance of panel radiosity model :math:`j` on panel :math:`i`.


 Parameters
 ----------
 luminosity_model : LuminosityModelSettings
     Settings for the luminosity model.
 Returns
 -------
 RadiationSourceModelSettings
     Object defining settings for source model irradiance






     )doc" );

    ///////////////////////////////////////////////////////////
    ///////////   RADIATION TARGET MODELS
    ///////////////////////////////////////////////////////////

    py::class_< tss::RadiationPressureTargetModelSettings,
                std::shared_ptr< tss::RadiationPressureTargetModelSettings > >(
            m,
            "RadiationPressureTargetModelSettings",
            R"doc(

         Base class for providing settings for properties of a radiation target (e.g. spacecraft), to be used in the context of (for instance) calculation of radiation pressure on spacecraft





      )doc" );

    m.def( "cannonball_radiation_target",
           &tss::cannonballRadiationPressureTargetModelSettingsWithOccultationMap,
           py::arg( "reference_area" ),
           py::arg( "radiation_pressure_coefficient" ),
           py::arg( "per_source_occulting_bodies" ) =
                   std::map< std::string, std::vector< std::string > >( ),
           R"doc(

 Function for cannonball radiation target


 Parameters
 ----------
 reference_area : float
     Cross-sectional area of cannonball [:math:`m^{2}`]
 radiation_pressure_coefficient : float
     Radiation pressure coefficient [-]
 per_source_occulting_bodies : Dict[str, List[str]]
     Names of bodies to occult the source as seen from this target
 Returns
 -------
 CannonballRadiationPressureTargetModelSettings
     Object defining settings for a cannonball radiation pressure target model






     )doc" );

    m.def( "panelled_radiation_target",
           &tss::paneledRadiationPressureTargetModelSettingsWithOccultationMap,
           py::arg( "source_to_target_occulting_bodies" ) = std::map< std::string, std::vector< std::string > >( ),
           py::arg( "maximum_number_of_pixels_per_source" ) = std::map< std::string, int >( ),
           R"doc(

 Function for creating settings for a paneled radiation pressure target model

 Function for creating settings for a paneled radiation pressure target model. Each source can have
 its own set of occulting bodies.
 This model requires the :attr:`~tudatpy.numerical_simulation.environment_setup.BodySettings.vehicle_shape_settings` of type :class:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.FullPanelledBodySettings` to be defined.
 The functions to define the panelled body settings are available in the :ref:`vehicle_systems` module.

 Parameters
 ----------
 source_to_target_occulting_bodies : dict[str, list[str]]
     Dictionary (source name -> list of occulting body names) of bodies to occult sources as seen from this target.
 maximum_number_of_pixels : Dict[str, int]
     Maximum number of pixels used in the self-shadowing algorithm per source, omitting a the value or setting it to zero equals to not considering self-shadowing for a given source.
 Returns
 -------
 RadiationPressureTargetModelSettings
     Object defining settings for a radiation pressure target model






     )doc" );
}

}  // namespace radiation_pressure
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
