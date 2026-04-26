/*    Copyright (c) 2010-2021, Delft University of Technology
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
#include "expose_light_time_corrections.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/ionosphereModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
// #include <pybind11/native_enum.h>

namespace tom = tudat::observation_models;
namespace tuc = tudat::unit_conversions;
namespace ti = tudat::interpolators;

namespace tudatpy
{
namespace estimation
{
namespace observable_models_setup
{

namespace light_time_corrections
{

void expose_light_time_corrections( py::module& m )
{
    py::enum_< tom::LightTimeCorrectionType >(
            m,
            "LightTimeCorrectionType",
            R"doc(Enum identifying each type of light-time correction registered on a link.

Used as a filter by :func:`~tudatpy.estimation.observations_setup.observations_dependent_variables.light_time_correction_components_dependent_variable`
to select which correction contributions are saved individually.)doc" )
            .value( "first_order_relativistic", tom::first_order_relativistic )
            .value( "function_wrapper_light_time_correction", tom::function_wrapper_light_time_correction )
            .value( "tabulated_tropospheric", tom::tabulated_tropospheric )
            .value( "saastamoinen_tropospheric", tom::saastamoinen_tropospheric )
            .value( "vmf3_tropospheric", tom::vmf3_tropospheric )
            .value( "vmf3o_tropospheric", tom::vmf3o_tropospheric )
            .value( "tabulated_ionospheric", tom::tabulated_ionospheric )
            .value( "jakowski_vtec_ionospheric", tom::jakowski_vtec_ionospheric )
            .value( "inverse_power_series_solar_corona", tom::inverse_power_series_solar_corona )
            .value( "ionex_vtec_ionospheric", tom::ionex_vtec_ionospheric )
            .value( "nequick2_ionospheric", tom::nequick2_ionospheric )
            .export_values( );

    py::class_< tom::LightTimeCalculatorBase, std::shared_ptr< tom::LightTimeCalculatorBase > >(
            m,
            "LightTimeCalculatorBase",
            R"doc(Non-templated view over a light-time calculator, primarily exposing access to the per-correction
breakdown from the last light-time evaluation.)doc" )
            .def( "get_current_light_time_correction_components",
                  &tom::LightTimeCalculatorBase::getCurrentLightTimeCorrectionComponents,
                  R"doc(Return the per-correction values cached during the last call to `setTotalLightTimeCorrection`.
The returned list's order matches `get_light_time_correction_list()`.)doc" )
            .def( "get_light_time_correction_list",
                  &tom::LightTimeCalculatorBase::getLightTimeCorrectionList,
                  R"doc(Return the list of light-time correction objects registered on this calculator.
Order corresponds to the values returned by `get_current_light_time_correction_components`.)doc" );

    py::class_< tom::LightTimeConvergenceCriteria, std::shared_ptr< tom::LightTimeConvergenceCriteria > >( m,
                                                                                                           "LightTimeConvergenceCriteria",
                                                                                                           R"doc(

         Base class to define criteria of light time convergence.

         Base class to define criteria of light time convergence.
         This class is not used for calculations of corrections, but is used for the purpose of defining the light time convergence criteria.
         Specific light time convergence criteria must be defined using an object derived from this class.
         Instances of this class are typically created via the :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings` function.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a LightTimeConvergenceCriteria object
             from tudatpy.estimation.observable_models_setup import light_time_corrections

             # Create Default Light Time Convergence Settings (no args specified = setting default arguments)
             light_time_convergence_settings = light_time_corrections.light_time_convergence_settings()

             # Show that it is an LightTimeConvergenceCriteria object.
             print(light_time_convergence_settings)




      )doc" );

    py::class_< tom::LightTimeCorrectionSettings, std::shared_ptr< tom::LightTimeCorrectionSettings > >( m,
                                                                                                         "LightTimeCorrectionSettings",
                                                                                                         R"doc(

         Base class to define light time correction settings.

         Base class to define light time correction settings.
         This class is not used for calculations of corrections, but is used for the purpose of defining the light time correction properties.
         Specific light time correction settings must be defined using an object derived from this class.

         Instances of this class are typically created via the
         :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction` function

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a LightTimeCorrectionSettings object
             from tudatpy.estimation.observable_models_setup import light_time_corrections, links

             # Create Link Ends dictionary
             link_ends = dict()
             link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
             link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

             # Create a Link Definition Object from link_ends dictionary
             Link_Definition_Object = links.LinkDefinition(link_ends)

             # Case 1: perturbing body (Earth) involved in the observations
             # In this case, Earth is a receiver, so the body’s state will be evaluated at the reception time.
             perturbing_body = ['Earth']
             doppler_observation_settings = light_time_corrections.first_order_relativistic_light_time_correction(perturbing_body)

             # Show that it is a LightTimeCorrectionSettings object.
             print(doppler_observation_settings)

             # Case 2: perturbing body (Sun) not involved in the observations
             # In this case, the body's state will be evaluated at the midpoint time between the transmission and reception events.
             perturbing_body = ['Sun']

             # Use: light_time_corrections.first_order_relativistic_light_time_correction to create a LightTimeCorrectionSettings object
             # Note: first_order_relativistic_light_time_correction only requires the perturbing list of bodies to be passed as arguments
             doppler_observation_settings = light_time_corrections.first_order_relativistic_light_time_correction(perturbing_body)

             # Show that it is an LightTimeCorrectionSettings object.
             print(doppler_observation_settings.transmitter_proper_time_rate_settings)
             print(dir(doppler_observation_settings))



      )doc" );

    py::enum_< tom::LightTimeFailureHandling >( m, "LightTimeFailureHandling", R"doc(

Enumeration of behaviour when failing to converge light-time with required settings.

Examples
--------
.. code-block:: python

    # Code snippet to print all available Light Time Failure Handling Types
    from tudatpy.estimation.observable_models_setup import light_time_corrections

    num_LightTimeFailureHandling_types = len(light_time_corrections.LightTimeFailureHandling.__members__)
    print(f'The length of all available Tudatpy Light Time Failure Handling Types is: {num_LightTimeFailureHandling_types}')

    # Print all available Observation Viability Types using the "name" property
    for i in range(num_LightTimeFailureHandling_types):
        print(i, light_time_corrections.LightTimeFailureHandling(i).name)



      )doc" )
            .value( "accept_without_warning", tom::LightTimeFailureHandling::accept_without_warning )
            .value( "print_warning_and_accept", tom::LightTimeFailureHandling::print_warning_and_accept )
            .value( "throw_exception", tom::LightTimeFailureHandling::throw_exception )
            .export_values( );

    m.def( "light_time_convergence_settings",
           &tom::lightTimeConvergenceCriteria,
           py::arg( "iterate_corrections" ) = false,
           py::arg( "maximum_number_of_iterations" ) = 50,
           py::arg( "absolute_tolerance" ) = TUDAT_NAN,
           py::arg( "failure_handling" ) = tom::accept_without_warning,
           R"doc(

 Function for creating convergence settings for solving the light-time equation.

 Function for creating convergence settings for solving the light-time equation. Computing the light time
 :math:`s=t_{R}-t_{T}` between two receiver :math:`R` and transmitter :math:`T` requires the iterative
 solution of the following equation:

 .. math::
     t_{R} - t_{T} = c\left(|\mathbf{r}_{R}(t_{R}) - \mathbf{r}_{T}(t_{T})| + \Delta s(t_{R}, t_{T}, \mathbf{r}_{R}(t_{R}), \mathbf{r}_{T}(t_{T}))\right)

 where either the reception time :math:`t_{R}` or the transmission time :math:`t_{T}` is kept fixed (reference link end time). The term :math:`\Delta s` contains any
 deviations in the light-time from straight-line propagation at speed of light (relativistic corrections, media corrections, etc.). The algorithm starts
 at :math:`t_{R}=t_{T}`, and uses this to evaluate the right-hand side of the above equation. This leads to a new value of :math:`t_{R}` or :math:`t_{T}` (depending on which is kept fixed)
 and the right-hand side is re-evaluated in a new iteration. The input to this function defines the settings for when the iteration will terminate.

 Parameters
 ----------
 iterate_corrections : bool, default = False
     Boolean denoting whether the terms :math:`\Delta s` are recomputed at each iteration or not. If false, the corrections are calculated only on the first iteration. Afterwards, the value
     is kept fixed until convergence. Once preliminarily converged, the algorithm recomputes :math:`\Delta s`, and continues the iteration (until proper convergence) while now recomputing
     :math:`\Delta s` each iteration. Setting this input to false is typically safe, and is computationally more efficient.

 maximum_number_of_iterations : int, default = 50
     Maximum number of iterations taken by the algorithm. If this number of iterations is reached without convergence (as defined by ``absolute_tolerance`` input),
     the behaviour of the algorithm is defined by the ``failure_handling`` input.

 absolute_tolerance : float, default = nan
     Difference in :math:`t_{R}-t_{T}` between two consecutive iterations below which the algorithm is considered to be converged. Default value is nan, which means the default value is taken.
     The default value depends on the time representation used (1 ps for float; 1 fs for Time class)

 failure_handling : LightTimeFailureHandling, default = accept_without_warning
     Input defines behaviour when failing to converge within the required number of iterations. NOTE: the default value should be overridden for high-accuracy applications

 Returns
 -------
 :class:`LightTimeConvergenceCriteria`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria` with the required settings.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the light_time_convergence_settings function
     from tudatpy.estimation.observable_models_setup import light_time_corrections

     # The light_time_convergence_settings function can be used with default inputs as just:
     light_time_convergence_settings = light_time_corrections.light_time_convergence_settings()
     # A LightTimeConvergenceCriteria object is returned
     print(light_time_convergence_settings)

     # Users can also specify the following input arguments:
     # iterate_corrections, maximum_number_of_iterations, absolute_tolerance, failure_handling.
     # Let's set the failure_handling argument to LightTimeFailureHandling.print_warning_and_accept (default was LightTimeFailureHandling.accept_without_warning)
     light_time_convergence_settings = light_time_corrections.light_time_convergence_settings(
         failure_handling = light_time_corrections.LightTimeFailureHandling.print_warning_and_accept
     )
     # Again, a LightTimeConvergenceCriteria object is returned
     print(light_time_convergence_settings)



     )doc" );

    m.def(
            "first_order_relativistic_light_time_correction",
            []( const std::vector< std::string >& perturbingBodies ) {
                // Force bending to always be false
                return tom::firstOrderRelativisticLightTimeCorrectionSettings( perturbingBodies, false );
            },
            py::arg( "perturbing_bodies" ),
            R"doc(

Function for creating settings for first-order relativistic light-time corrections.

Function for creating settings for first-order relativistic light-time corrections:  These corrections account for the delay in light travel time caused by stationary point masses, calculated up to
:math:`c^{-2}` according to general relativity (e.g., Moyer, 2000 Eq 8.55). A key consideration in the model is the time at which the states of the perturbing bodies are evaluated. This depends on their involvement in the observation link ends:

* 1. **Perturbing Body as a Link End:** If the perturbing body (e.g., Earth) is directly involved in the observation (e.g., as the location of a transmitter or receiver):

    - The body's state is evaluated at the **transmission time** if it acts as the **transmitter**.
    - The body's state is evaluated at the **reception time** if it acts as the **receiver**.

* 2. **Perturbing Body Not as a Link End:** If the perturbing body is not part of the observation link ends, its state is evaluated at the **midpoint time** between the transmission and reception events.

Parameters
----------
perturbing_bodies : List[str]
    A list containing the names of the bodies due to which the light-time correction is to be taken into account.

Returns
-------
:class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
    Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` configured to include
    first-order relativistic light-time corrections.

Examples
--------
.. code-block:: python

    # Code Snippet to showcase the use of the first_order_relativistic_light_time_correction function
    from tudatpy.estimation.observable_models_setup import light_time_corrections, links

    # Create Link Ends dictionary
    link_ends = dict()
    link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
    link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

    # Create a Link Definition Object from link_ends dictionary
    Link_Definition_Object = links.LinkDefinition(link_ends)

    # The function first_order_relativistic_light_time_correction() requires a list of strings (perturbing body/bodies) as input
    # and a boolean value for bending (default is True).
    perturbing_body = ['Earth']
    doppler_observation_settings = light_time_corrections.first_order_relativistic_light_time_correction(perturbing_body)

    # Show that it returns a LightTimeCorrectionSettings object.
    print(doppler_observation_settings)

     )doc" );

    m.def(
            "approximated_second_order_relativistic_light_time_correction",
            []( const std::vector< std::string >& perturbingBodies ) {
                // Force bending to always be true
                return tom::firstOrderRelativisticLightTimeCorrectionSettings( perturbingBodies, true );
            },
            py::arg( "perturbing_bodies" ),

            R"doc(

Function for creating settings for Moyer, 2000 Eq 8.55 approximated second-order relativistic light-time corrections.

Function for creating settings for approximated second-order relativistic light-time corrections:  These corrections account for the delay in light travel time caused by stationary point masses, calculated up to
:math:`c^{-2}` according to general relativity ( Moyer, 2000 Eq 8.55; correction term for Sun) and it includes the bending of light due to the perturbing body. A key consideration in the model is the time at which the states of the perturbing bodies are evaluated. This depends on their involvement in the observation link ends:

* 1. **Perturbing Body as a Link End:** If the perturbing body (e.g., Earth) is directly involved in the observation (e.g., as the location of a transmitter or receiver):

    - The body's state is evaluated at the **transmission time** if it acts as the **transmitter**.
    - The body's state is evaluated at the **reception time** if it acts as the **receiver**.

* 2. **Perturbing Body Not as a Link End:** If the perturbing body is not part of the observation link ends, its state is evaluated at the **midpoint time** between the transmission and reception events.

Parameters
----------
perturbing_bodies : List[str]
    A list containing the names of the bodies due to which the light-time correction is to be taken into account.

Returns
-------
:class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
    Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` configured to include
    approximated second-order relativistic light-time corrections.

Examples
--------
.. code-block:: python

    # Code Snippet to showcase the use of the first_order_relativistic_light_time_correction function
    from tudatpy.estimation.observable_models_setup import light_time_corrections, links

    # Create Link Ends dictionary
    link_ends = dict()
    link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
    link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

    # Create a Link Definition Object from link_ends dictionary
    Link_Definition_Object = links.LinkDefinition(link_ends)

    # The function first_order_relativistic_light_time_correction() requires a list of strings (perturbing body/bodies) as input
    perturbing_body = ['Earth']
    doppler_observation_settings = light_time_corrections.approximated_second_order_relativistic_light_time_correction(perturbing_body)

    # Show that it returns a LightTimeCorrectionSettings object.
    print(doppler_observation_settings)

     )doc" );

    py::enum_< tom::TroposphericMappingModel >( m, "TroposphericMappingModel", R"doc(No documentation found.)doc" )
            .value( "simplified_chao", tom::TroposphericMappingModel::simplified_chao )
            .value( "niell", tom::TroposphericMappingModel::niell )
            .value( "vmf3", tom::TroposphericMappingModel::vmf3 );

    py::enum_< tom::WaterVaporPartialPressureModel >(
            m, "WaterVaporPartialPressureModel", "enum.IntEnum", R"doc(No documentation found.)doc" )
            .value( "tabulated", tom::WaterVaporPartialPressureModel::tabulated )
            .value( "bean_and_dutton", tom::WaterVaporPartialPressureModel::bean_and_dutton );

    m.def( "dsn_tabulated_tropospheric_light_time_correction",
           &tom::tabulatedTroposphericCorrectionSettings,
           py::arg( "file_names" ),
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           py::arg( "mapping_model" ) = tom::TroposphericMappingModel::niell,
           R"doc(No documentation found.)doc" );

    m.def( "saastamoinen_tropospheric_light_time_correction",
           &tom::saastamoinenTroposphericCorrectionSettings,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           py::arg( "mapping_model" ) = tom::TroposphericMappingModel::niell,
           py::arg( "water_vapor_partial_pressure_model" ) = tom::WaterVaporPartialPressureModel::tabulated,
           R"doc(No documentation found.)doc" );

    m.def( "dsn_tabulated_ionospheric_light_time_correction",
           &tom::tabulatedIonosphericCorrectionSettings,
           py::arg( "file_names" ),
           py::arg( "spacecraft_name_per_id" ) = std::map< int, std::string >( ),
           py::arg( "quasar_name_per_id" ) = std::map< int, std::string >( ),
           py::arg( "reference_frequency" ) = 2295e6,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           R"doc(No documentation found.)doc" );

    m.def( "jakowski_ionospheric_light_time_correction",
           &tom::jakowskiIonosphericCorrectionSettings,
           py::arg( "ionosphere_height" ) = 400.0e3,
           py::arg( "first_order_delay_coefficient" ) = 40.3,
           py::arg( "solar_activity_data_path" ) = tudat::paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt",
           py::arg( "geomagnetic_pole_latitude" ) = tuc::convertDegreesToRadians( 80.9 ),
           py::arg( "geomagnetic_pole_longitude" ) = tuc::convertDegreesToRadians( -72.6 ),
           py::arg( "use_utc_for_local_time_computation" ) = false,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           R"doc(

Function for creating settings for Jakowski VTEC ionospheric light-time corrections.

Computes the ionospheric delay using the Jakowski et al. (2011) analytical VTEC model with a Modified Single Layer
Model (MSLM) mapping function. The model computes VTEC at the sub-ionospheric point as a function of geomagnetic
latitude, local time, and solar activity (F10.7 index), then maps to slant TEC via :math:`\text{STEC} = \text{VTEC} / \cos(z')`,
where :math:`z'` is the zenith angle at the ionospheric pierce point.

This correction is suitable for ground-to-GNSS geometries where the receiver is below the ionosphere.

Parameters
----------
ionosphere_height : float, default = 400.0e3
    Height of the ionospheric single-layer shell [m].
first_order_delay_coefficient : float, default = 40.3
    Coefficient relating STEC to range delay: :math:`\Delta\rho = c_1 \cdot \text{STEC} / f^2` [m :sup:`3` /s :sup:`2`].
solar_activity_data_path : str
    Path to space weather data file containing F10.7 solar radio flux values.
geomagnetic_pole_latitude : float, default = 80.9 deg (in radians)
    Geodetic latitude of the geomagnetic north pole [rad].
geomagnetic_pole_longitude : float, default = -72.6 deg (in radians)
    Geodetic longitude of the geomagnetic north pole [rad].
use_utc_for_local_time_computation : bool, default = False
    If True, use UTC for local time computation instead of TDB.
body_with_atmosphere_name : str, default = "Earth"
    Name of the body with the ionosphere.

Returns
-------
:class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
    Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
    configured for Jakowski VTEC ionospheric corrections.

           )doc" );

    // IONEX-based VTEC correction (DEPRECATED)
    m.def( "ionex_ionospheric_light_time_correction",
           &tom::ionexIonosphericCorrectionSettings,
           py::arg( "body_with_ionosphere_name" ),
           py::arg( "ionosphere_height" ),
           py::arg( "first_order_delay_coefficient" ) = 40.3,
           R"doc(

.. deprecated::
    Use :func:`nequick2_ionospheric_light_time_correction` instead, which provides physically
    correct path-integrated ionospheric corrections for all geometries including space-based receivers.

Function for creating settings for IONEX-based ionospheric light-time corrections using the
Modified Single Layer Model (MSLM). This approximation maps VTEC to STEC via
:math:`\text{STEC} = \text{VTEC} / \cos(z')` and is only valid when the receiver is below the ionosphere
and the transmitter is above it. For receivers inside the ionosphere (e.g., ISS at ~420 km), this
mapping is geometrically meaningless.

Parameters
----------
body_with_ionosphere_name : str
    Name of the body with the ionosphere (e.g., ``"Earth"``).
ionosphere_height : float
    Height of the ionospheric single-layer shell [m]. Typically 450e3 for IONEX maps.
first_order_delay_coefficient : float, default = 40.3
    Coefficient relating STEC to range delay [m :sup:`3` /s :sup:`2`].

           )doc" );

    // NeQuick-2 path-integrated ionospheric correction
    m.def( "nequick2_ionospheric_light_time_correction",
           []( bool useIonexRescaling,
               double firstOrderDelayCoefficient,
               int quadratureOrder,
               const std::string& ccirDataPath,
               const std::string& solarActivityDataPath,
               double ionexRmsBiasTecu ) {
               return tom::nequick2IonosphericCorrectionSettings(
                   "Earth", useIonexRescaling, firstOrderDelayCoefficient,
                   quadratureOrder, ccirDataPath, solarActivityDataPath, ionexRmsBiasTecu );
           },
           py::arg( "use_ionex_rescaling" ) = true,
           py::arg( "first_order_delay_coefficient" ) = 40.3,
           py::arg( "quadrature_order" ) = 50,
           py::arg( "ccir_data_path" ) = "",
           py::arg( "solar_activity_data_path" ) = "",
           py::arg( "ionex_rms_bias_tecu" ) = 0.0,
           R"doc(

Function for creating settings for NeQuick-2 path-integrated ionospheric light-time corrections.

Computes the ionospheric delay by numerically integrating the NeQuick-2 electron density profile (ITU-R P.531)
along the transmitter-receiver ray path using Gauss-Legendre quadrature. The delay is given by
:math:`\Delta\tau = c_1 \cdot \text{STEC} / f^2 / c`, where STEC is the path-integrated electron content.

Unlike the MSLM-based corrections (:func:`ionex_ionospheric_light_time_correction`,
:func:`jakowski_ionospheric_light_time_correction`), this correction works for **arbitrary geometries**:

- Ground-to-GNSS (classical)
- Ground-to-LEO (low elevation, receiver below ionosphere)
- Space-to-space (both endpoints above surface)
- **Receivers embedded in the ionosphere** (e.g., ISS at ~420 km, ACES MWL)

When ``use_ionex_rescaling`` is True and an IONEX ionosphere model has been loaded on the body, the NeQuick-2
F2-layer peak density (NmF2) is rescaled so that the NeQuick-2 vertical TEC column matches the IONEX VTEC at the
ionospheric pierce point. This preserves the physical profile shape while anchoring the total electron content to
the GNSS-derived IONEX measurement. When set to False, the model runs in free climatology mode using F10.7 only.

Parameters
----------
body_with_ionosphere_name : str, default = "Earth"
    Name of the body with the ionosphere.
use_ionex_rescaling : bool, default = True
    If True, rescale NmF2 so that the NeQuick-2 vertical TEC matches the IONEX VTEC at the ionospheric
    pierce point. Requires an ionosphere model (from IONEX) to be set on the body via
    :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.set_ionosphere_model_from_ionex`.
    If the body has no ionosphere model set, falls back to free-running NeQuick-2 without error.
first_order_delay_coefficient : float, default = 40.3
    Coefficient relating STEC to range delay: :math:`\Delta\rho = c_1 \cdot \text{STEC} / f^2` [m :sup:`3` /s :sup:`2`].
quadrature_order : int, default = 50
    Number of Gauss-Legendre quadrature nodes for numerical integration along the ray path.
    50 nodes provide sub-0.1% accuracy for smooth ionospheric profiles.
ccir_data_path : str, default = ""
    Path to directory containing the NeQuick-2 CCIR coefficient files (``ccir11.asc`` through ``ccir22.asc``)
    and the modified dip latitude grid (``modip.asc``). If empty, uses the default path
    (``~/.tudat/resource/ionosphere/nequick2/``).
solar_activity_data_path : str, default = ""
    Path to space weather data file containing F10.7 solar radio flux values.
    If empty, uses the default space weather file.
ionex_rms_bias_tecu : float, default = 0.0
    Additional bias term [TECU] added in quadrature to the IONEX RMS uncertainty.
    The total 1-sigma VTEC uncertainty is computed as
    :math:`\sigma_{\text{total}} = \sqrt{\text{RMS}_{\text{IONEX}}^2 + \text{bias}^2}`.
    A typical value of 2-5 TECU accounts for systematic errors in the IONEX maps
    not captured by the formal RMS. Set to 0 to use only the IONEX RMS.
    The resulting correction uncertainty is stored internally and can be retrieved
    as a dependent variable after simulation.

Returns
-------
:class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
    Instance of the :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings`
    configured for NeQuick-2 path-integrated ionospheric corrections.

Examples
--------
.. code-block:: python

    from tudatpy.estimation.observable_models_setup import light_time_corrections

    # Load IONEX maps onto Earth body (required for IONEX-constrained mode)
    light_time_corrections.set_ionosphere_model_from_ionex(
        data_files=["path/to/ionex_file.INX"],
        bodies=bodies
    )

    # Create NeQuick-2 correction with IONEX rescaling (default)
    nequick2_correction = light_time_corrections.nequick2_ionospheric_light_time_correction()

    # Or without IONEX rescaling (free-running NeQuick-2 climatology)
    nequick2_free = light_time_corrections.nequick2_ionospheric_light_time_correction(
        use_ionex_rescaling=False
    )

References
----------
- ITU-R Recommendation P.531-15 (2023). *Ionospheric propagation data and prediction methods
  required for the design of satellite networks and systems.*
- B. Nava, P. Coisson, and S.M. Radicella (2008). A new version of the NeQuick ionosphere
  electron density model. *J. Atmos. Sol.-Terr. Phys.*, 70(15), 1856-1862.
  `doi:10.1016/j.jastp.2008.01.015 <https://doi.org/10.1016/j.jastp.2008.01.015>`_
- M. Hernandez-Pajares et al. (2009). The IGS VTEC maps: a reliable source of ionospheric
  information since 1998. *J. Geodesy*, 83(3-4). `doi:10.1007/s00190-008-0266-1 <https://doi.org/10.1007/s00190-008-0266-1>`_

           )doc" );

    // VMF3 Tropospheric correction
    m.def( "vmf3_tropospheric_light_time_correction",
           &tom::vmf3TroposphericCorrectionSettings,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           py::arg( "use_gradient_correction" ) = true,
           py::arg( "tropospheric_mapping_model" ) = tom::TroposphericMappingModel::vmf3,
           R"doc(Create VMF3 tropospheric light time correction settings.)doc" );

    m.def( "vmf3o_tropospheric_light_time_correction",
           &tom::vmf3oTroposphericCorrectionSettings,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           py::arg( "use_gradient_correction" ) = true,
           py::arg( "tropospheric_mapping_model" ) = tom::TroposphericMappingModel::vmf3,
           py::arg( "observation_wavelength_nm" ) = 532.0,
           R"doc(
Create VMF3o (optical) tropospheric light time correction settings.

This uses VMF3 mapping with VMF3o-specific coefficient handling and wavelength-dependent scaling.
           )doc" );

    m.def( "inverse_power_series_solar_corona_light_time_correction",
           &tom::inversePowerSeriesSolarCoronaCorrectionSettings,
           py::arg( "coefficients" ) = std::vector< double >{ 1.3e14, 0.5e12 },
           py::arg( "positive_exponents" ) = std::vector< double >{ 6.0, 2.0 },
           py::arg( "delay_coefficient" ) = 40.3,
           py::arg( "sun_body_name" ) = "Sun",
           R"doc(No documentation found.)doc" );

    py::class_< tom::VtecCalculator, std::shared_ptr< tom::VtecCalculator > >( m, "VtecCalculator" );

    py::class_< tom::JakowskiVtecCalculator, std::shared_ptr< tom::JakowskiVtecCalculator >, tom::VtecCalculator >(
            m, "JakowskiVtecCalculator" )
            .def( py::init< std::function< double( double ) >, std::function< double( double ) >, bool >( ),
                  py::arg( "sun_declination_function" ),
                  py::arg( "f10p7_function" ),
                  py::arg( "use_utc_time_for_local_time" ) = false )
            .def( "calculate_vtec",
                  &tudat::observation_models::JakowskiVtecCalculator::calculateVtec,
                  py::arg( "time" ),
                  py::arg( "sub_ionospheric_point" ) );

    py::class_< tom::GlobalIonosphereModelVtecCalculator,
                std::shared_ptr< tom::GlobalIonosphereModelVtecCalculator >,
                tom::VtecCalculator >( m, "GlobalIonosphereModelVtecCalculator",
                R"doc(

VTEC calculator that wraps an IonosphereModel (e.g., from IONEX data).

Computes the vertical total electron content at a given sub-ionospheric point
by querying the underlying ionosphere model. This is used internally by the
IONEX ionospheric light-time correction, but can also be used directly to
sample VTEC on a grid for visualization (e.g., global VTEC maps).

The sub-ionospheric point is specified as a geodetic position vector
``[altitude_m, latitude_rad, longitude_rad]``.

Examples
--------
.. code-block:: python

    from tudatpy.estimation.observable_models_setup import light_time_corrections

    ionosphere_model = bodies.get_body("Earth").get_ionosphere_model()
    vtec_calc = light_time_corrections.GlobalIonosphereModelVtecCalculator(ionosphere_model)

    # Sample VTEC at 450 km, 45 deg N, 15 deg E
    import numpy as np
    geodetic = np.array([450e3, np.deg2rad(45.0), np.deg2rad(15.0)])
    vtec = vtec_calc.calculate_vtec(epoch_tdb, geodetic)  # returns VTEC in el/m^2

                )doc" )
            .def( py::init< std::shared_ptr< tudat::environment::IonosphereModel > >( ),
                  py::arg( "ionosphere_model" ),
                  R"doc(

Construct a VTEC calculator from an ionosphere model.

Parameters
----------
ionosphere_model : IonosphereModel
    The ionosphere model to wrap (e.g., from ``bodies.get_body("Earth").get_ionosphere_model()``).

                  )doc" )
            .def( "calculate_vtec",
                  &tudat::observation_models::GlobalIonosphereModelVtecCalculator::calculateVtec,
                  py::arg( "time" ),
                  py::arg( "sub_ionospheric_point" ),
                  R"doc(

Calculate the vertical total electron content at a sub-ionospheric point.

Parameters
----------
time : float
    Time [seconds since J2000 TDB].
sub_ionospheric_point : numpy.ndarray
    Geodetic position as [altitude_m, latitude_rad, longitude_rad].

Returns
-------
float
    Vertical TEC in electrons/m^2 (multiply by 1e-16 to convert to TECU).

                  )doc" );

    m.def( "set_vmf_troposphere_data",
           &tom::setVmfTroposphereCorrections,
           py::arg( "data_files" ),
           py::arg( "file_has_meteo" ),
           py::arg( "file_has_gradient" ),
           py::arg( "bodies" ),
           py::arg( "set_troposphere_data" ) = true,
           py::arg( "set_meteo_data" ) = true,
           py::arg( "interpolator_settings" ) = ti::cubicSplineInterpolation( ),
           py::arg( "retrieve_mapping_internally" ) = false,
           R"doc(
Set VMF/VMF3/VMF3o troposphere (and optional meteo) data in Earth ground stations.

If ``retrieve_mapping_internally`` is ``True``, station-name matching first attempts direct key matching and then
internally maps ILRS station code <-> DOMES identifiers using the default ILRS SINEX ``SITE/ID`` registry.
           )doc" );

    m.def( "set_ionosphere_model_from_ionex",
           &tom::setIonosphereModelFromIonex,
           py::arg( "data_files" ),
           py::arg( "bodies" ),
           py::arg( "interpolator_settings" ) = std::shared_ptr< ti::InterpolatorSettings >( ),
           py::arg( "station_subset" ) = std::vector< std::pair< std::string, std::string > >( ),
           py::arg( "subset_padding_deg" ) = 30.0,
           R"doc(

Load IONEX global ionosphere map data and set the resulting ionosphere model on the Earth body.

Reads one or more IONEX files (e.g., from IGS, CODE, JPL, ESA, UPC analysis centres), builds a 3D
interpolator over [time, latitude, longitude], and attaches the resulting
:class:`TabulatedIonosphereModel` to the ``"Earth"`` body. This model provides VTEC values at arbitrary
locations and times, which can then be used by the NeQuick-2 ionospheric light-time correction.

When ``station_subset`` is provided, the IONEX map is spatially cropped to a bounding box around the
specified ground stations (plus ``subset_padding_deg`` on each side). This significantly reduces memory
usage and interpolation time, especially for long arcs with high-resolution IONEX products. The station
geodetic positions are looked up from the ``bodies`` object.

Parameters
----------
data_files : list[str]
    List of paths to IONEX files (``*.INX`` or ``*.ionex``). Files from any IGS analysis centre
    (COD, EMR, ESA, IGS, JPL, UPC) and any temporal resolution (15 min to 2 hours) are supported.
bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    System of bodies. The ionosphere model will be set on the ``"Earth"`` body.
interpolator_settings : InterpolatorSettings, optional
    Settings for the 3D multi-linear interpolator. If not provided, defaults to multi-linear
    interpolation with hunting algorithm and boundary value extrapolation.
station_subset : list[tuple[str, str]], default = []
    List of ``(body_name, station_name)`` pairs defining the ground stations around which to crop
    the IONEX grid. When empty (default), the full global grid is loaded. Example:
    ``[("Earth", "OPMT"), ("Earth", "WTZR")]``.
subset_padding_deg : float, default = 30.0
    Padding in degrees of latitude and longitude around the bounding box of the station locations.

Examples
--------
.. code-block:: python

    # Load full global IONEX grid (default)
    light_time_corrections.set_ionosphere_model_from_ionex(
        data_files=["path/to/ionex.INX"], bodies=bodies)

    # Load only a regional subset around two stations (saves memory)
    light_time_corrections.set_ionosphere_model_from_ionex(
        data_files=["path/to/ionex.INX"], bodies=bodies,
        station_subset=[("Earth", "OPMT"), ("Earth", "WTZR")],
        subset_padding_deg=30.0)

           )doc" );
}

}  // namespace light_time_corrections
}  // namespace observable_models_setup
}  // namespace estimation
}  // namespace tudatpy
