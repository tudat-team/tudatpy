/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_light_time_corrections.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

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
            .export_values( );

    py::enum_< tom::WaterVaporPartialPressureModel >( m, "WaterVaporPartialPressureModel", R"doc(No documentation found.)doc" )
            .value( "tabulated", tom::WaterVaporPartialPressureModel::tabulated )
            .value( "bean_and_dutton", tom::WaterVaporPartialPressureModel::bean_and_dutton )
            .export_values( );

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
           R"doc(No documentation found.)doc" );
    
    // IONEX-based VTEC correction
    m.def( "ionex_ionospheric_light_time_correction",
        &tom::ionexIonosphericCorrectionSettings,
        py::arg( "body_with_ionosphere_name" ),
        py::arg( "ionosphere_height" ),
        py::arg( "first_order_delay_coefficient" ) = 40.3,
        R"doc(Create IONEX-based ionospheric light time correction settings.)doc" );

    // VMF3 Tropospheric correction
    m.def( "vmf3_tropospheric_light_time_correction",
        &tom::vmf3TroposphericCorrectionSettings,
        py::arg( "body_with_atmosphere_name" ) = "Earth",
        py::arg( "use_gradient_correction" ) = true,
        py::arg( "tropospheric_mapping_model" ) = tom::vmf3,
        R"doc(Create VMF3 tropospheric light time correction settings.)doc" );

    m.def( "inverse_power_series_solar_corona_light_time_correction",
           &tom::inversePowerSeriesSolarCoronaCorrectionSettings,
           py::arg( "coefficients" ) = std::vector< double >{ 1.3e14, 0.5e12 },
           py::arg( "positive_exponents" ) = std::vector< double >{ 6.0, 2.0 },
           py::arg( "delay_coefficient" ) = 40.3,
           py::arg( "sun_body_name" ) = "Sun",
           R"doc(No documentation found.)doc" );

    py::class_<tom::VtecCalculator, std::shared_ptr<tom::VtecCalculator>>(m, "VtecCalculator");

    py::class_<tom::JakowskiVtecCalculator, std::shared_ptr<tom::JakowskiVtecCalculator>, tom::VtecCalculator>(
        m, "JakowskiVtecCalculator")
        .def(py::init<
            std::function<double(double)>,
            std::function<double(double)>,
            bool>(),
            py::arg("sun_declination_function"),
            py::arg("f10p7_function"),
            py::arg("use_utc_time_for_local_time") = false)
        .def( "calculate_vtec",
            &tudat::observation_models::JakowskiVtecCalculator::calculateVtec,
            py::arg( "time" ),
            py::arg( "sub_ionospheric_point" ));

    py::class_<tom::GlobalIonosphereModelVtecCalculator,
            std::shared_ptr<tom::GlobalIonosphereModelVtecCalculator>,
            tom::VtecCalculator>(
        m, "GlobalIonosphereModelVtecCalculator")
        .def(py::init< std::shared_ptr<tudat::environment::IonosphereModel> >(),
            py::arg("ionosphere_model"))
        .def( "calculate_vtec",
          &tudat::observation_models::GlobalIonosphereModelVtecCalculator::calculateVtec,
          py::arg( "time" ),
          py::arg( "sub_ionospheric_point" ));

    m.def( "set_vmf_troposphere_data",
           &tom::setVmfTroposphereCorrections,
           py::arg( "data_files" ),
           py::arg( "file_has_meteo" ),
           py::arg( "file_has_gradient" ),
           py::arg( "bodies" ),
           py::arg( "set_troposphere_data" ) = true,
           py::arg( "set_meteo_data" ) = true,
           py::arg( "interpolator_settings" ) = ti::cubicSplineInterpolation( ) );

    m.def( "set_ionosphere_model_from_ionex",
       &tom::setIonosphereModelFromIonex,
       py::arg( "data_files" ),
       py::arg( "bodies" ),
       py::arg( "interpolator_settings" ) = std::shared_ptr< ti::InterpolatorSettings >( ) );
}

}
}
}
}