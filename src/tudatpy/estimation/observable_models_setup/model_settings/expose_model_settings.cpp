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
#include "expose_model_settings.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace estimation
{
namespace observable_models_setup
{

namespace model_settings
{

void expose_model_settings( py::module& m )
{
    py::enum_< tom::ObservableType >( m, "ObservableType", R"doc(

Enumeration of available observable types.

Examples
--------
.. code-block:: python

    # Code snippet to print all available Observable Types
    from tudatpy.estimation.observable_models_setup import model_settings

    num_observable_types = len(model_settings.ObservableType.__members__)
    print(f'The length of all available Tudatpy Observable Types is: {num_observable_types}')

    # Print all available Observable Types using the "name" property
    for i in range(num_observable_types):
        print(i, model_settings.ObservableType(i).name)




      )doc" )
            .value( "one_way_range_type", tom::ObservableType::one_way_range )
            .value( "n_way_range_type", tom::ObservableType::n_way_range )
            .value( "angular_position_type", tom::ObservableType::angular_position )
            .value( "relative_angular_position_type", tom::ObservableType::angular_position )
            .value( "position_observable_type", tom::ObservableType::position_observable )
            .value( "velocity_observable_type", tom::ObservableType::velocity_observable )
            .value( "relative_position_observable_type", tom::ObservableType::relative_position_observable )
            .value( "one_way_instantaneous_doppler_type", tom::ObservableType::one_way_doppler )
            .value( "one_way_averaged_doppler_type", tom::ObservableType::one_way_differenced_range )
            .value( "two_way_instantaneous_doppler_type", tom::ObservableType::two_way_doppler )
            .value( "n_way_averaged_doppler_type", tom::ObservableType::n_way_differenced_range )
            .value( "euler_angle_313_observable_type", tom::ObservableType::euler_angle_313_observable )
            .value( "dsn_one_way_averaged_doppler_type", tom::ObservableType::dsn_one_way_averaged_doppler )
            .value( "dsn_n_way_averaged_doppler_type", tom::ObservableType::dsn_n_way_averaged_doppler )
            .value( "doppler_measured_frequency_type", tom::ObservableType::doppler_measured_frequency )
            .value( "dsn_n_way_range_type", tom::ObservableType::dsn_n_way_range )
            .value( "differenced_time_of_arrival_type", tom::ObservableType::differenced_time_of_arrival )
            .export_values( );

    py::class_< tom::DopplerProperTimeRateSettings, std::shared_ptr< tom::DopplerProperTimeRateSettings > >(
            m,
            "DopplerProperTimeRateSettings",
            R"doc(

         Base class to define proper time rate settings.

         Base class to define proper time rate settings (at a single link end) for instantaneous Doppler observation model settings.
         Specific proper time rate settings must be defined using an object derived from this class.
         The derived classes are made accessible via dedicated functions.





      )doc" );

    py::class_< tom::ObservationModelSettings, std::shared_ptr< tom::ObservationModelSettings > >( m, "ObservationModelSettings", R"doc(

         Base class to define settings of observation models.

         Base class to define settings of observation models.
         Observation model settings define at least the type and geometry of a given observation.
         They can furthermore set observation biases and/or light-time corrections.
         Simple observation models settings that are fully characterised by these elements can be managed by this base class.
         Instances of this class are typically created via functions, such as
         :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`, :func:`~tudatpy.estimation.observable_models_setup.model_settings.cartesian_position`,
         :func:`~tudatpy.estimation.observable_models_setup.model_settings.angular_position`, etc.
         Model settings for specific observation models that require additional information such as integration time, retransmission time, etc. must be defined using an object derived from this class.
         The derived classes are made accessible through further functions.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationModelSettings object
             from tudatpy.estimation.observable_models_setup import links, model_settings

             # Create Link Ends dictionary
             link_ends = dict()
             link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
             link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

             # Create a Link Definition Object from link_ends dictionary
             Link_Definition_Object = links.LinkDefinition(link_ends)

             # Create minimal ObservationModelSettings object (only required Link_Definition_Object argument is passed)
             # Other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) are set by default
             observation_settings = model_settings.one_way_range(Link_Definition_Object)

             # Show that it is an ObservationModelSettings object.
             print(observation_settings)




      )doc" );

    py::class_< tom::OneWayDopplerObservationModelSettings,
                std::shared_ptr< tom::OneWayDopplerObservationModelSettings >,
                tom::ObservationModelSettings >( m,
                                                 "OneWayDopplerObservationModelSettings",
                                                 R"doc(

         Derived Class for defining the settings of one-way instantaneous Doppler observation models.

         Derived Class for defining the settings of one-way instantaneous Doppler observation models.
         Settings object can account for additional observation model aspects such as light time corrections and proper time rate settings.
         Instances of this class can be created via the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_instantaneous` function.
         Associated base class: :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a OneWayDopplerObservationModelSettings object
             from tudatpy.estimation.observable_models_setup import links, model_settings

             # Create Link Ends dictionary
             link_ends = dict()
             link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
             link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

             # Create a Link Definition Object from link_ends dictionary
             Link_Definition_Object = links.LinkDefinition(link_ends)

             # Use: model_settings.one_way_doppler_instantaneous to create a OneWayDopplerObservationModelSettings object (only required Link_Definition_Object argument is passed)
             # Other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings, proper time rate) are set by default
             doppler_observation_settings = model_settings.one_way_doppler_instantaneous(Link_Definition_Object)

             # Show that it is an OneWayDopplerObservationModelSettings object.
             print(doppler_observation_settings)




      )doc" );

    py::class_< tom::NWayRangeObservationModelSettings,
                std::shared_ptr< tom::NWayRangeObservationModelSettings >,
                tom::ObservationModelSettings >( m, "NWayRangeObservationModelSettings", R"doc(No documentation found.)doc" );

    m.def( "one_way_range",
           &tom::oneWayRangeSettings,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a one-way range observable.

 Function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
 a single-valued observable :math:`h_{_{\text{1-range}}}` as follows (in the unbiased case):

 .. math::
    h_{_{\text{1-range}}}(t_{R},t_{T})=|\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})| + \Delta s

 where :math:`\mathbf{r}_{R}`, :math:`\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
 of receiver and transmitter. The term :math:`\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
 The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
 As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the light-time equation (see :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings` for details)


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)

 light_time_convergence_settings : :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the one-way observable.


 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the one_way_range function
     from tudatpy.estimation.observable_models_setup import links, model_settings

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
     link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

     # Create a Link Definition Object from link_ends dictionary. This will be the input to the function.
     Link_Definition_Object = links.LinkDefinition(link_ends)

     # Create minimal ObservationModelSettings object (only required Link_Definition_Object argument is passed)
     # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
     observation_settings = model_settings.one_way_range(Link_Definition_Object)

     # Show that this returns an ObservationModelSettings object.
     print(observation_settings)



     )doc" );

    m.def( "two_way_range",
           &tom::twoWayRangeSimple,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a two-way range observable.

 Same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_range`, with :math:`n=2`. This function is provided
 for convenience.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter`, `retransmitter` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined

 light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 light_time_convergence_settings : :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the two_way_range function
     from tudatpy.estimation.observable_models_setup import links, model_settings

     # two_way_range() takes a Link Definition Object as input to the function.
     # Note: as for this case, transmitter, retransmitter and receiver are required to define the Link Ends dictionary
     link_ends = dict()
     link_ends[links.transmitter] = links.body_origin_link_end_id("Earth")
     link_ends[links.retransmitter] = links.body_origin_link_end_id("Delfi-C3")
     link_ends[links.receiver] = links.body_origin_link_end_id("Earth")

     # Create the LinkDefinition object
     Link_Definition_Object = links.LinkDefinition(link_ends)

     # Create minimal ObservationModelSettings object (only required Link_Definition_Object argument is passed)
     # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
     observation_settings = model_settings.two_way_range(Link_Definition_Object)

     # Show that two_way_range() returns an NWayRangeObservationModelSettings object.
     print(observation_settings)



     )doc" );

    m.def( "two_way_range_from_one_way_links",
           &tom::twoWayRange,
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a two-way range observable.

 Same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_range_from_one_way_links`, with :math:`n=2`. This function is provided
 for convenience.


 Parameters
 ----------
 one_way_range_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` ]
     List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
     The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
     ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
     ``transmitter`` and ``receiver`` of the second entry of this list.

 bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the two_way_range_from_one_way_links function
     from tudatpy.estimation.observable_models_setup import links, model_settings

     # two_way_range_from_one_way_links() takes a list of ObservationModelSettings objects
     # Note: as for this case, transmitter, retransmitter and receiver are required to define the Link Ends dictionary
     link_ends = dict()
     link_ends[links.transmitter] = links.body_origin_link_end_id("Earth")
     link_ends[links.retransmitter] = links.body_origin_link_end_id("Delfi-C3")
     link_ends[links.receiver] = links.body_origin_link_end_id("Earth")

     # Create the LinkDefinition object to be used as input
     Link_Definition_Object = links.LinkDefinition(link_ends) # define LinkDefinition object
     two_way_range_observation_settings_list = [model_settings.two_way_range(Link_Definition_Object)] # define (minimal) NWayRangeObservationModelSettings object
     two_way_range_one_way_link_settings = model_settings.two_way_range_from_one_way_links(two_way_range_observation_settings_list)

     # Show that two_way_range_from_one_way_links() returns an NWayRangeObservationModelSettings object.
     print(two_way_range_one_way_link_settings)



     )doc" );

    m.def( "n_way_range",
           &tom::nWayRangeSimple,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a n-way range observable.

 Function for creating observation model settings of n-way range type observables. The associated observation model creates
 a scalar observable :math:`h_{_{\text{N-range}}}` by combining together a series :math:`n` one-way range observations
 (see :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`).

 By default, the reception time of the :math:`i^{th}` one-way range in this n-way range is set as the
 transmission time of the :math:`(i+1)^{th}` one-way range. A retransmission delay may be defined by ancilliary settings
 (see :func:`~func.estimation.observations_setup.ancillary_settings.n_way_range_ancilliary_settings`) when creating observation
 simulation setings (see `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations/simulating-observations.html#defining-observation-simulation-settings>`_).

 For this function, the settings for each constituent one-way range (with the exception of the link end identifiers) are equal.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     ``transmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
     as a ``reflector1``, ``reflector2``, .... (with the number of reflectors to be defined by the user). For a two-way range, the
     ``retransmitter`` may be used instead of ``reflector1`` (which are identical in the role they play)

 light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 light_time_convergence_settings : :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings` class.


 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the n_way_range function
     from tudatpy.estimation.observable_models_setup import links, model_settings

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
     link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

     # n_way_range() takes a Link Definition Object as input to the function.
     # Let's create it from link_ends
     Link_Definition_Object = links.LinkDefinition(link_ends)

     # Create minimal ObservationModelSettings object (only required Link_Definition_Object argument is passed)
     # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
     observation_settings = model_settings.n_way_range(Link_Definition_Object)

     # Show that n_way_range() returns an NWayRangeObservationModelSettings object.
     print(observation_settings)





     )doc" );

    m.def( "n_way_range_from_one_way_links",
           &tom::nWayRange,
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a n-way range observable.

 Function for creating observation model settings of n-way range type observables, for a single link definition. The
 implementation is the same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_range`, with the difference
 that the input is provided as separate settings for the constituent one-way ranges.

 Parameters
 ----------
 one_way_range_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` ]
     List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way range observable.
     The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
     ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` (n-1) and ``receiver`` are defined by the
     ``transmitter`` and ``receiver`` of the :math:`\text{n}^{th}` entry of this list.

 bias_settings : :class:`~tudatpy.estimation.observable_models_setup.biases.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`~tudatpy.estimation.observable_models_setup.model_settings.NWayRangeObservationModelSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the n_way_range_from_one_way_links function
     from tudatpy.estimation.observable_models_setup import links, model_settings

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
     link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")

     # n_way_range_from_one_way_links() takes 1) a list of ObservationModelSettings objects and 2) bias as input (default is None)
     # Let's create it.
     Link_Definition_Object = links.LinkDefinition(link_ends) # define LinkDefinition object
     n_way_observation_settings_list = [model_settings.n_way_range(Link_Definition_Object)] # define (minimal) ObservationModelSettings object

     n_way_from_one_link_observation_settings = model_settings.n_way_range_from_one_way_links(n_way_observation_settings_list, bias_settings = None)

     # Show that n_way_range_from_one_way_links() returns an NWayRangeObservationModelSettings object.
     print(n_way_from_one_link_observation_settings)



     )doc" );

    m.def( "angular_position",
           &tom::angularPositionSettings,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an angular position observable.

 Function for creating observation model settings of angular position type observables (as right ascension :math:`\alpha` and declination :math:`\delta`),
 for a single link definition. The associated observation model creates an observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of type two as follows (in the unbiased case):

 .. math::
    \Delta\mathbf{r}&=\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})\\
    \tan\alpha&=\frac{\Delta r_{y}}{\Delta r_{x}}\\
    \delta&=\frac{\Delta r_{z}}{\sqrt{\Delta r_{x}^{2}+\Delta r_{y}^{2}}}\\
    \mathbf{h}_{_{\text{ang.pos.}}}&=[\alpha;\delta]

 The relative position vector :math:`\Delta\mathbf{r}` is computed identically as described for the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`
 The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
 environment.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the angular position observable.






     )doc" );

    m.def( "relative_angular_position",
           &tom::relativeAngularPositionSettings,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a relative angular position observable.

 Function for creating observation model settings of relative angular position type observables (as right ascension difference :math:`\Delta \alpha` and declination difference :math:`\Delta \delta`).
 It computes two :func:`~angular_position` observations :math:`[\alpha_{1};\delta_{1}]` and :math:`[\alpha_{2};\delta_{2}]` with the same receiver and different transmitters,
 and computes the observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of size two as follows (in the unbiased case):

 .. math::

    \mathbf{h}_{_{\text{rel.ang.pos.}}} = [\alpha_{2}-\alpha_{1};\delta_{2}-\delta_{1}]

 The relative angular position observable can be used for optical astrometry, optical navigation, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
 environment.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     ``transmitter``, `transmitter2` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the relative angular position observable.

     )doc" );

    m.def( "cartesian_position",
           &tom::positionObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a Cartesian position observable.

 Function for creating observation model settings of Cartesian position type observables.
 Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
 This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the ``observed_body`` defined by the `link_ends` input.
 The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     ``observed_body`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the cartesian position observable.






     )doc" );

    m.def( "relative_cartesian_position",
           &tom::relativePositionObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a relative Cartesian position observable.

 Function for creating observation model settings of relative Cartesian position type observables.
 Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
 This observable provides the inertial Cartesian position of the ``observed_body``, w.r.t. the ``observer`` defined by the `link_ends` input.
 The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     ``observed_body`` and ``observer`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the relative Cartesian position observable.






     )doc" );

    m.def( "cartesian_velocity",
           &tom::velocityObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a Cartesian velocity observable.

 Function for creating observation model settings of Cartesian position type observables.
 Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
 This observable provides the inertial (w.r.t. global frame origin) Cartesian velocity of the ``observed_body`` defined by the `link_ends` input.
 The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` velocity


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     ``observed_body`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the cartesian velocity observable.






     )doc" );

    m.def( "euler_angles_313",
           &tom::eulerAngle313ObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(No documentation found.)doc" );

    m.def( "one_way_doppler_instantaneous",
           &tom::oneWayOpenLoopDoppler,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "transmitter_proper_time_rate_settings" ) = nullptr,
           py::arg( "receiver_proper_time_rate_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false,
           R"doc(

 Function for creating settings for a one-way instantaneous Doppler observable.

 Function for creating settings for a one-way instantaneous Doppler observable for a single link definition. The associated observation model creates
 a scalar observable :math:`h_{_{\text{1-Dopp.}}}` as follows (in the unbiased case):

 .. math::
    h_{_{\text{1-Dopp.}}}=c\left(\frac{d\tau_{T}}{dt_{T}}\frac{dt_{T}}{dt_{R}}\frac{dt_{R}}{d\tau_{R}}-1\right)

 where :math:`t` and :math:`\tau` denote coordinate and proper time of the transmitter T and receiver R, respectively.
 The receiver and transmitter position and coordinate time are computed identically as described for the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`.

 The observable may be non-dimensionalized by the speed of light :math:`c` (using the ``normalized_with_speed_of_light`` input), which results in the observable being equal to the
 received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`. For the observation model provding a model for the received
 frequency, use :func:`~tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous_frequency`.

 This observable represents the 'instantaneous (non-integrated)' Doppler observable, such as those obtained from open-loop observations.
 It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the
 :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` should be used).

 The coordinate time derivative :math:`\frac{dt_{T}}{dt_{R}}` is always computed when generating this observable as:

 .. math::
   \frac{dt_{T}}{dt_{R}}=


 Settings for the proper time
 rates :math:`\frac{d\tau}{dt}` can be specified by the user through the ``transmitter_proper_time_rate_settings`` and ``receiver_proper_time_rate_settings``
 arguments (inputs, see Parameters). Whenever these are left empty, the proper time rates are omitted (set to 1.0).

 The above formulation is sensitive to numerical noise, since two almost equal quantities are subtracted. To mitigate this, we expand each
 term in the above equation into a Taylor series. By setting :math:`\Delta\tau=\tau-t`, we have:

 .. math::
    \frac{d\tau}{dt}&=1+\frac{d\Delta\tau}{dt}\\
    \frac{dt}{d\tau}&\approx 1+\sum_{i=0}^{N}(-1)^i\frac{d\Delta\tau}{dt}

 and we set:

 .. math::
   F_{\tau}&=\frac{d\tau_{T}}{dt_{T}}\frac{dt_{R}}{d\tau_{R}}-1\\
           &\approx \frac{d\Delta\tau_{T}}{dt_{T}}\cdot\left(1+\sum_{i=0}^{N}(-1)^i\frac{d\Delta\tau}{dt}\right)

 Similarly, we can expand the coordinate time derivative as, setting :math:`T=t_{R}-t_{T}`:

 .. math::

    s_{T}&=-\frac{dT}{d\mathbf{r}_{T}}\cdot\left(\frac{\mathbf{v_{T}}}{c}\right)\\
    s_{R}&=\frac{dT}{d\mathbf{r}_{R}}\cdot\left(\frac{\mathbf{v_{R}}}{c}\right)\\
    \frac{1}{1+s_{T}}&\approx 1+\sum_{i=0}^{N}s_{T}

 then:

 .. math::

    F_{t}&=\frac{dt_{T}}{dt_{R}}-1\\
         &\approx -s_{R}+\left(\sum_{i=0}^{N}s_{T}\right)(1-s_{R}\\

 Using these approximations, the complete observable is then computed as:

 .. math::
   h_{_{\text{1-Dopp.}}} = F_{t}F_{\tau}+F_{t}+F_{\tau}

 At present, we have hardcoded :math:`N=3`, which should be sufficient for all conceivable solar system cases.

 Note that, when calculating the partial derivatives of the observables for estimation, only the first-order Taylor series terms are retained.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     ``transmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 transmitter_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
     Settings for computing the transmitter proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

 receiver_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
     Settings for computing the receiver proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 normalized_with_speed_of_light : bool, default = false
     Option to non-dimensionalize the observable with speed of light :math:`c`

 Returns
 -------
 :class:`OneWayDopplerObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`OneWayDopplerObservationModelSettings` class defining the settings for the one-way open doppler observable observable.






     )doc" );

    m.def( "two_way_doppler_instantaneous_from_one_way_links",
           py::overload_cast< const std::shared_ptr< tom::OneWayDopplerObservationModelSettings >,
                              const std::shared_ptr< tom::OneWayDopplerObservationModelSettings >,
                              const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::twoWayOpenLoopDoppler ),
           py::arg( "uplink_doppler_settings" ),
           py::arg( "downlink_doppler_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

Function for creating settings for a two-way instantaneous Doppler observable.

Function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The associated observation model creates
a scalar observable :math:`h_{_{\text{2-Dopp.}}}`

.. math::

   h_{_{\text{2-Dopp.}}}=c\left(\frac{d\tau_{T}}{dt_{T}}\frac{dt_{T}}{dt_{R}}\frac{dt_{R}}{d\tau_{R}}-1\right)

with the link end :math:`T` and :math:`R` denoting the transmitter and receiver of a two-way link. The observable is computed
from a concatenation of a one-way Doppler uplink observation
:math:`h_{_{\text{1-Dopp.}\uparrow}}` and a one-way Doppler downlink obervation :math:`h_{_{\text{1-Dopp.}\downarrow}}` (see
:func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_instantaneous`) as (with the term :math:`c` omitted in the 1-way computations):

.. math::

   h_{_{\text{2-Dopp.}}}&=c\left(\left(h_{_{\text{1-Dopp.}\uparrow}}+1\right)\left(h_{_{\text{1-Dopp.}\downarrow}}+1\right)-1\right)
                        &=c\left(h_{_{\text{1-Dopp.}\uparrow}}h_{_{\text{1-Dopp.}\downarrow}}+h_{_{\text{1-Dopp.}\uparrow}}+h_{_{\text{1-Dopp.}\downarrow}}\right)



 Parameters
 ----------
 uplink_doppler_settings : :class:`OneWayDopplerObservationModelSettings`
     Settings for uplink leg of one-way observable, created using :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_instantaneous`

 downlink_doppler_settings : :class:`OneWayDopplerObservationModelSettings`
     Settings for downlink leg of one-way observable, created using :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_instantaneous`

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the full observation, default is none (unbiased observation). Note that,
     even if no bias is applied to the two-way observable, the constituent one-way observables may still be biased.

 Returns
 -------
 :class:`TwoWayDopplerObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`TwoWayDopplerObservationModelSettings` class defining the settings for the two-way open doppler observable.

     )doc" );

    m.def( "two_way_doppler_instantaneous",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                              const bool >( &tom::twoWayOpenLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false,
           R"doc(


Function for creating settings for a two-way instantaneous Doppler observable.

Implementation as in :func:`~tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous_from_one_way_links`,
but without full control over the up- and downlink settings (specifically, no control over proper time rate settings: :math:`d\tau/dt=1` is used
for this factory function)


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    ``transmitter``, ``retransmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
    Settings for convergence of the light-time

normalized_with_speed_of_light : bool, default = false
    Option to non-dimensionalize the observable with speed of light :math:`c`

 Returns
 -------
 :class:`OneWayDopplerObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`OneWayDopplerObservationModelSettings` class defining the settings for the one-way open doppler observable observable.




)doc" );

    m.def( "one_way_doppler_averaged",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >( &tom::oneWayClosedLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a one-way averaged Doppler observable.

 Function for creating observation model settings for one-way averaged Doppler observables, for a single link definition. The associated observation model creates
 a scalar observable :math:`h_{_{1-\bar{\text{Dopp}}}}` as follows (in the unbiased case):

 .. math::
    h_{_{1-\bar{\text{Dopp}}}}&=c\int_{t-\Delta t/2}^{t+\Delta t/2}\frac{dt_{T}}{dt_{R}}d\bar{t}\\
                              &=\frac{h_{_{\text{1-range}}}(t_{R}=t+\Delta t/2,t_{T})-h_{_{\text{1-range}}}(t_{R}=t+\Delta t/2,t_{T})}{\Delta t}

 where, in the latter formulation (which is the one that is implemented), the observable is referenced to the receiver time. This averaged Doppler observable
 is computed as the difference of two one-way range observables (see :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`),
 with the reference time shifted by :math:`\Delta t` (in TDB time, with the time tag in the center of this interval). As such, it is sensitive to numerical errors for small :math:`\Delta t`

 The integration time :math:`\Delta t` is defined in the ancilliary settings (see
 `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations/simulating-observations.html#defining-observation-simulation-settings>`_)
 when simulating the observations (with 60 s as default).

 Note that this observation model is a simplified version of the Doppler data as generated by the DSN/ESTRACK, which is suitable for simulation
 studies

 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     ``transmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived `OneWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "two_way_doppler_averaged",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::twoWayDifferencedRangeObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implementation is
 analogous to the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` observable. But, in the present case
 the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

 The integration time :math:`\Delta t` is defined in the ancilliary settings (see
 `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/creating-observations/simulating-observations.html#defining-observation-simulation-settings>`_)
 when simulating the observations (with 60 s as default).

 Note that this observation model is a simplified version of the Doppler data as generated by the DSN/ESTRACK, which is suitable for simulation
 studies

 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     ``transmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
     as a `retransmitter1`, ``retransmitter2``, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "two_way_doppler_averaged_from_one_way_links",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::twoWayDifferencedRangeObservationSettings ),
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
 analogous to the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` observable. But, in the present case
 the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

 The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).

 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     ``transmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
     as a ``retransmitter1``, ``retransmitter2``, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "n_way_doppler_averaged",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >( &tom::nWayDifferencedRangeObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implementation is
 analogous to the :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_doppler_averaged` observable. But, in the present case
 the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

 The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     ``transmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined, as well
     as a ``retransmitter1``, ``retransmitter2``, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.


     )doc" );

    m.def( "n_way_doppler_averaged_from_one_way_links",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >( &tom::nWayDifferencedRangeObservationSettings ),
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition.
 The implementation is the same as :func:`~tudatpy.estimation.observable_models_setup.model_settings.n_way_doppler_averaged`, with the difference
 that the constituent one-way range observables may have different settings.


 Parameters
 ----------
 one_way_range_settings : List[ :class:`ObservationModelSettings` ]
     List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way averaged range rate observable.
     The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
     ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` (n-1) and ``receiver`` are defined by the
     ``transmitter`` and ``receiver`` of the :math:`\text{n}^{th}` entry of this list.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived `~tudatpy.estimation.observable_models_setup.model_settings.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.




     )doc" );

    m.def( "dsn_n_way_doppler_averaged",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                              const bool >( &tom::dsnNWayAveragedDopplerObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "subtract_doppler_signature" ) = true,
           R"doc(

Function for creating settings for a two- or three-way Doppler observation model, as generated by DSN/ESTRACK.

Function for creating settings for a two- or three-way Doppler observation model, as generated by DSN/ESTRACK. This
model is used for processing of real radio tracking data, for instance from ODF, TNF or TFMS files (see TODO-ADDLINK).
It requires a frequency calculator to be set for the transmitter.

The Doppler observable at time tag :math:`t_{R}` (in TDB) is computed as follows, based on :cite:p:`moyer2005`:

* First the time tag is computed from TDB to UTC (denoted here :math:`\bar{t}_{R}`).
* The receiver start time :math:`\bar{t}_{R,s}` and end time :math:`\bar{t}_{R,e}` (in UTC) of the observation are then computed as :math:`\bar{t}_{R,s}=\bar{t}_{R}-\Delta t/2` and :math:`\bar{t}_{R,e}=\bar{t}_{R}+\Delta t/2`.
* These times are then converted back to TDB times :math:`t_{R,s}` and :math:`t_{R,e}`.
* The two-way light time from the receiver :math:`R` to the transmitter :math:`T` is then computed, for both :math:`t_{R,s}` and :math:`t_{R,e}` as reception times, to obtain the associated signal transmission times (in TDB) :math:`t_{T,s}` and :math:`t_{T,e}`. Note that, if frequency-dependent light-time corrections (such as ionosphere) are used, the transmission frequency :math:`f_{T}` is evaluated using the transmission epoch from the first iteration of the light-time solution (this incurs a very minor error on the frequencies used in the light-time solution).
* These transmission times are converted to UTC to obtain :math:`\bar{t}_{T,s}` and :math:`\bar{t}_{T,e}`.
* The transmitter frequency integral is then computed from the following:

.. math::

   F = \int_{\bar{t}_{T,s}}^{\bar{t}_{T,e}} f_{T}(t)\, dt

* The final observable :math:`h` is then computed from:

.. math::

   h = M_{2_{R}}f_{\text{ref}} \pm \frac{M_{2}}{\Delta t} F

where the sign is defined by the ``subtract_doppler_signature`` input (negative if ``true``; positive if ``false``).

Computing this observable requires the following ancillary settings to be set:

* ``doppler_integration_time`` — Defines the observable integration time at the receiver (in UTC seconds).
* ``doppler_reference_frequency`` — The Doppler reference frequency :math:`f_{\text{ref}}`. If not used, must be set to 0.
* ``frequency_bands`` — The frequency band for the up- and downlink, used to compute the turnaround ratio :math:`M_{2}` on the spacecraft.
* ``reception_reference_frequency_band`` — The receiver reference frequency band used to compute :math:`M_{2_{R}}`.

Typically, these quantities are loaded automatically from the tracking data files.

Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    ``transmitter``, ``retransmitter1`` and ``receiver`` and  :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
    Settings for convergence of the light-time

subtract_doppler_signature : bool, default = true
    Boolean defining if scaled frequency integral is to be subtracted or added (see formulation above)

Returns
-------
:class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
    Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived class defining the settings for doppler observable.

.)doc" );

    m.def( "two_way_doppler_instantaneous_frequency",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::dopplerMeasuredFrequencyObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

Function for creating settings for a two-way instantaneous Doppler frequency observable.

Implementation as in :func:`~tudatpy.estimation.observable_models_setup.model_settings.two_way_doppler_instantaneous`, but
with the observable in Hz rather than m/s. It requires a frequency calculator to be set for the transmitter.

Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    ``transmitter``, ``retransmitter`` and ``receiver`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
    Settings for convergence of the light-time

 Returns
 -------
 :class:`OneWayDopplerObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived :class:`OneWayDopplerObservationModelSettings` class defining the settings for the one-way open doppler observable observable.



)doc" );

    m.def( "dsn_n_way_range",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >( &tom::dsnNWayRangeObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(
Function for creating settings for a two-way range observation model, as generated by DSN/ESTRACK.

Function for creating settings for a two-way range observation model, as generated by DSN/ESTRACK. This
model is used for processing of real radio tracking data, for instance from ODF, TNF or TFMS files (see TODO-ADDLINK).
It requires a frequency calculator to be set for the transmitter.

The range observable at time tag :math:`t_{R}` (in TDB) is computed as follows, based on :cite:p:`moyer2005`:

* First the time tag is computed from TDB to UTC (denoted here :math:`\bar{t}_{R}`).
* The two-way light time from the receiver :math:`R` to the transmitter :math:`T` is then computed to obtain :math:`t_{T}`. Note that, if frequency-dependent light-time corrections (such as ionosphere) are used, the transmission frequency :math:`f_{T}` is evaluated using the transmission epoch from the first iteration of the light-time solution (this incurs a very minor error on the frequencies used in the light-time solution).
* This transmission time is converted to UTC to obtain :math:`\bar{t}_{T}`.
* The transmitter frequency integral is then computed from the following, with :math:`K` a conversion factor that depends on the up- and downlink frequency (1/2 for S up and X down; :math:`221/(2\cdot 749)` for X up and X down; :math:`221/(3599 \cdot 2.0)` for X up and Ka down).:

.. math::

   F = K \int_{\bar{t}_{T}}^{\bar{t}_{R}} f_{T}(t)\, dt



* The final observable :math:`h` is then computed from, where :math:`n` is the lowest ranging component (defined in the ancillary settings, see below).

.. math::

   h = \mod\left(F, 2^{n+6}\right)

This observable gives the range in “range units” (RU), modulo the factor :math:`2^{n+6}`. The conversion factor from
range units to meters depends on the transmitted frequency and can therefore not (if the signal is ramped) be computed
a priori. The conversion factor is equal to:

.. math::

  \frac{c}{K \cdot f_{T}(\bar{t}_{T})}

Computing this observable requires the following ancillary settings to be set:

* ``sequential_range_lowest_ranging_component`` — The lowest sequential ranging component.
* ``frequency_bands`` — The frequency band for the up- and downlink, used to compute the turnaround ratio :math:`M_{2}` on the spacecraft.

Typically, these quantities are loaded automatically from the tracking data files.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    ``transmitter``, ``retransmitter1`` and ``receiver`` and  :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
    Settings for convergence of the light-time

Returns
-------
:class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
    Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` derived class defining the settings for this range observable.

)doc" );

    m.def( "differenced_time_of_arrival",
           &tom::differencedTimeOfArrivalObservationSettings,
           py::arg( "link_ends" ),
           py::arg( "time_difference_time_scale" ) = tba::tdb_scale,
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a time difference of arrival observation model

 Function for creating settings for a time difference of arrival observation model, computed from the time difference :math:`\Delta t`
 at which the same transmitted signal is received by two different link ends. For a signal transmitted at :math:`t_{T}` and received by two
 receivers at :math:`t_{R,1}` and :math:`t_{R,2}`, respectively, the observable :math:`h` is computed from:

 .. math::

    h = t_{R,1} - t_{R,2}

 From the reference link (typically :math:`t_{R,1}`), the times at the other two link ends are computed as described for the
 :func:`~tudatpy.estimation.observable_models_setup.model_settings.one_way_range`. Using :math:`t_{R,1}` as reference, the time
 :math:`t_{T}` is computed from which :math:`t_{R,2}` is subsequently calculated, combining two light-time calculations (one from
 receiver to transmitter, and one from transmitter to receiver2).

 While the observation model performs all operations in TDB scale (including the time tag of teh observation), the model
 can convert :math:`t_{R,1}` and :math:`t_{R,2}` to a different time scale (typically UTC) to compute :math:`h` by using the ``time_difference_time_scale`` input.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     ``transmitter``, ``receiver`` and ``receiver2`` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` entries to be defined.

 link_ends : TimeScales, default = tdb_scale
     Time scale in which the epochs :math:`t_{R,1}` and :math:`t_{R,2}` are expressed before differencing to compute the observation


 light_time_correction_settings : List[ :class:`~tudatpy.estimation.observable_models_setup.light_time_corrections.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationModelSettings`
     Instance of the :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` class defining the settings for the differenced time of arrival model.


)doc" );

    //////////////////////////////////////////// DEPRECATED
    ///////////////////////////////////////////////

    m.def( "one_way_open_loop_doppler",
           &tom::oneWayOpenLoopDoppler,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "transmitter_proper_time_rate_settings" ) = nullptr,
           py::arg( "receiver_proper_time_rate_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false );

    m.def( "two_way_open_loop_doppler_from_one_way_links",
           py::overload_cast< const std::shared_ptr< tom::OneWayDopplerObservationModelSettings >,
                              const std::shared_ptr< tom::OneWayDopplerObservationModelSettings >,
                              const std::shared_ptr< tom::ObservationBiasSettings > >( &tom::twoWayOpenLoopDoppler ),
           py::arg( "uplink_doppler_settings" ),
           py::arg( "downlink_doppler_settings" ),
           py::arg( "bias_settings" ) = nullptr );

    m.def( "two_way_open_loop_doppler",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                              const bool >( &tom::twoWayOpenLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false );

    m.def( "one_way_closed_loop_doppler",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >( &tom::oneWayClosedLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ) );

    m.def( "one_way_closed_loop_doppler",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >( &tom::oneWayClosedLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ) );

    m.def( "doppler_measured_frequency",
           py::overload_cast< const tom::LinkDefinition&,
                              const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::dopplerMeasuredFrequencyObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) = std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(No documentation found.)doc" );
}

}  // namespace model_settings
}  // namespace observable_models_setup
}  // namespace estimation
}  // namespace tudatpy