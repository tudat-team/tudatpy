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
#include "expose_biases.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"

#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tom = tudat::observation_models;

namespace tudatpy
{
namespace estimation
{
namespace observable_models_setup
{

namespace biases
{

void expose_biases( py::module& m )
{
    
    py::class_< tom::ObservationBiasSettings, std::shared_ptr< tom::ObservationBiasSettings > >( m, "ObservationBiasSettings", R"doc(

         Base class to defining observation bias settings.

         Base class to defining observation bias settings.
         Specific observation bias settings must be defined using an object derived from this class.
         Instances of this class are typically created via the
         :func:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` or :func:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` function.


         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationBiasSettings object
             # using absolute and relative bias settings
             from tudatpy.numerical_simulation.estimation_setup import observation
             import numpy as np

             bias_array = np.array([1e-2])

             # Use absolute_bias function
             absolute_bias_settings = observation.absolute_bias(bias_array)
             # Show that it is an ObservationBiasSettings object.
             print(absolute_bias_settings)

             # Use relative_bias function
             relative_bias_settings = observation.relative_bias(bias_array)
             # Show that it is an ObservationBiasSettings object.
             print(relative_bias_settings)
      )doc" );

      m.def( "clock_induced_bias",
           &tom::clockInducedBias,
           py::arg( "body_name" ),
           py::arg( "station_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "absolute_bias",
           &tom::constantAbsoluteBias,
           py::arg( "bias_value" ),
           R"doc(

 Function for creating settings for an absolute observation bias.

 Function for creating settings for an absolute observation bias. When calculating the observable value, applying this setting
 will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

 .. math::
    \tilde{h}=h+K

 where :math:`K` is the `bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the addition is component-wise.


 Parameters
 ----------
 bias_value : numpy.ndarray
     A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` defining the settings for a constant, absolute observation bias.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the absolute_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function absolute_bias() requires a numpy.array of bias values in input
     bias_array = np.array([1e-2])
     absolute_bias_settings = observation.absolute_bias(bias_array)

     # Show that it returns an ObservationBiasSettings object.
     print(absolute_bias_settings)



     )doc" );

    m.def( "relative_bias",
           &tom::constantRelativeBias,
           py::arg( "bias_value" ),
           R"doc(

 Function for creating settings for a relative observation bias.

 Function for creating settings for a relative observation bias. When calculating the observable value, applying this setting
 will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

 .. math::
    \tilde{h}=h(1+K)

 where :math:`K` is the`bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the multiplication is component-wise.


 Parameters
 ----------
 bias_value : numpy.ndarray
     A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 Returns
 -------
 :class:`ConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` class,
     defining the settings for a constant, relative observation bias.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the relative_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function relative_bias() requires a numpy.array of bias values in input
     bias_array = np.array([1e-2])
     relative_bias_settings_settings = observation.relative_bias(bias_array)

     # Show that it returns an ObservationBiasSettings object.
     print(relative_bias_settings_settings)




     )doc" );

    m.def( "arcwise_absolute_bias",
           py::overload_cast< const std::vector< double >&, const std::vector< Eigen::VectorXd >&, const tom::LinkEndType >(
                   &tom::arcWiseAbsoluteBias ),
           py::arg( "arc_start_times" ),
           py::arg( "bias_values" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise absolute observation biases.

 Function for creating settings for arc-wise absolute observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 arc_start_times : List[ float ]
     List containing starting times for each arc.

 bias_values : List[ numpy.ndarray ]
     List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_absolute_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_absolute_bias() requires:
     # 1) an arc_start_times list ,2) a numpy.array of bias values in input, 3) reference_link_end_type
     # Let's simulate two arcs
     arc_start_times = [0, 60] # define start time in seconds
     arcwise_bias_array = [np.array([1e-2]), np.array([2e-2])] # set arc bias
     reference_link_end_type = observation.receiver # set bias at receiving link end
     arcwise_absolute_bias_settings = observation.arcwise_absolute_bias(arc_start_times, arcwise_bias_array, observation.receiver)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_absolute_bias_settings)




     )doc" );

    m.def( "arcwise_absolute_bias_per_time",
           py::overload_cast< const std::map< double, Eigen::VectorXd >&, const tom::LinkEndType >( &tom::arcWiseAbsoluteBias ),
           py::arg( "bias_values_per_start_time" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise absolute observation biases.

 Function for creating settings for arc-wise absolute observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
     Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
     The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

 # Code Snippet to showcase the use of the arcwise_absolute_bias function
 from tudatpy.numerical_simulation.estimation_setup import observation
 import numpy as np

     # The function arcwise_absolute_bias_settings_per_time() requires:
     # 1) a dictionary with times as keys and bias values as values ,2) a reference_link_end_type
     # Let's simulate two arcs
     bias_value_per_start_time = dict()
     bias_value_per_start_time[0] = np.array([1e-2])
     bias_value_per_start_time[60] = np.array([2e-2])
     reference_link_end_type = observation.receiver # set bias at receiving link end

     arcwise_absolute_bias_settings_per_time = observation.arcwise_absolute_bias_per_time(bias_value_per_start_time, reference_link_end_type)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_absolute_bias_settings_per_time)



     )doc" );

    m.def( "arcwise_relative_bias",
           py::overload_cast< const std::vector< double >&, const std::vector< Eigen::VectorXd >&, const tom::LinkEndType >(
                   &tom::arcWiseRelativeBias ),
           py::arg( "arc_start_times" ),
           py::arg( "bias_values" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise relative observation biases.

 Function for creating settings for arc-wise relative observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 arc_start_times : List[ float ]
     List containing starting times for each arc.

 bias_values : List[ numpy.ndarray ]
     List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_relative_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_relative_bias() requires:
     # 1) an arc_start_times list ,2) a numpy.array of bias values in input, 3) reference_link_end_type
     # Let's simulate two arcs
     arc_start_times = [0, 60] # define start time in seconds
     arcwise_bias_array = [np.array([1e-2]), np.array([2e-2])] # set arc bias
     reference_link_end_type = observation.receiver # set bias at receiving link end
     arcwise_relative_bias_settings = observation.arcwise_relative_bias(arc_start_times, arcwise_bias_array, reference_link_end_type)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_relative_bias_settings)




     )doc" );

    m.def( "arcwise_relative_bias_per_time",
           py::overload_cast< const std::map< double, Eigen::VectorXd >&, const tom::LinkEndType >( &tom::arcWiseRelativeBias ),
           py::arg( "bias_values_per_start_time" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise relative observation biases.

 Function for creating settings for arc-wise relative observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
     Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
     The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_relative_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_relative_bias_per_time() requires:
     # 1) a dictionary with times as keys and bias values as values ,2) a reference_link_end_type
     # Let's simulate two arcs
     bias_value_per_start_time = dict()
     bias_value_per_start_time[0] = np.array([1e-2])
     bias_value_per_start_time[60] = np.array([2e-2])
     reference_link_end_type = observation.receiver # set bias at receiving link end

     arcwise_relative_bias_settings_per_time = observation.arcwise_relative_bias_per_time(bias_value_per_start_time, reference_link_end_type)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_relative_bias_settings_per_time)





     )doc" );

    m.def( "time_drift_bias",
           &tom::constantTimeDriftBias,
           py::arg( "bias_value" ),
           py::arg( "time_link_end" ),
           py::arg( "ref_epoch" ),
           R"doc(

 Function for creating settings for a time-drift bias.

 Function for creating settings for a time-drift bias.
 This bias setting generates the configuration for applying a constant time-drift bias to an observation model.

 Parameters
 ----------
 bias_value : numpy.ndarray
     Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
     assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 time_link_end : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used the current time.

 ref_epoch : float
     Defines the reference epoch at which the effect of the time drift is initialised.

 Returns
 -------
 :class:`constantTimeDriftBias`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.constantTimeDriftBias` class,
     defining the settings for a constant, relative observation bias.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the time_drift_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function time_drift_bias() requires a numpy.array of bias value, time_link_end and ref_epoch as inputs
     bias_array = np.array([1e-2])
     time_link_end = observation.receiver
     ref_epoch = 0
     time_drift_bias_settings = observation.time_drift_bias(bias_array, time_link_end, ref_epoch)

     # Show that it returns an ObservationBiasSettings object.
     print(time_drift_bias_settings)






     )doc" );

    m.def( "arc_wise_time_drift_bias",
           py::overload_cast< const std::vector< Eigen::VectorXd >&,
                              const std::vector< double >&,
                              const tom::LinkEndType,
                              const std::vector< double >& >( &tom::arcWiseTimeDriftBias ),
           py::arg( "bias_value" ),
           py::arg( "arc_start_times" ),
           py::arg( "time_link_end" ),
           py::arg( "ref_epochs" ),
           R"doc(

 Function for creating settings for arc-wise time-drift biases.

 Function for creating settings for arc-wise time-drift biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.time_drift_bias` setting only through the option of setting the `bias_value` (time drift bias) to a different values for each arc.


 Parameters
 ----------
 bias_value : numpy.ndarray
     Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
     assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 arc_start_times : List[ float ]
     List containing starting times for each arc.

 time_link_end : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used the current time.

 ref_epochs : List[ float ]
     List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_time_drift_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_time_drift_bias() requires:
     # 1) an arc_start_times list ,2) a numpy.array of bias values in input, 3) reference_link_end_type, 4) list of reference epochs
     # Let's simulate two arcs
     arc_start_times = [0, 60] # define start time in seconds
     arcwise_bias_array = [np.array([1e-2]), np.array([2e-2])] # set arc bias
     reference_link_end_type = observation.receiver # set bias at receiving link end
     ref_epochs = [0,60]
     arcwise_time_drift_bias_settings = observation.arc_wise_time_drift_bias(arcwise_bias_array, arc_start_times, observation.receiver, ref_epochs)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_time_drift_bias_settings)




     )doc" );

    m.def( "arc_wise_time_drift_bias_per_time",
           py::overload_cast< const std::map< double, Eigen::VectorXd >&, const tom::LinkEndType, const std::vector< double > >(
                   &tom::arcWiseTimeDriftBias ),
           py::arg( "bias_value_per_start_time" ),
           py::arg( "time_link_end" ),
           py::arg( "ref_epochs" ),
           R"doc(

 Function for creating settings for arc-wise time-drift biases.

 Function for creating settings for arc-wise time-drift biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.time_drift_bias` setting only through the option of setting the `bias_value` (time drift bias) to a different values for each arc.

 Parameters
 ----------
 bias_value_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
     Dictionary, in which the time bias value vectors for each arc are directly mapped to the starting times of the respective arc.
     The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 time_link_end : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used the current time.

 ref_epochs : List[ float ]
     List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

 Returns
 -------
 :class:`ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_time_drift_bias_per_time function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_time_drift_bias_per_time() requires:
     # 1) an arc_start_times list ,2)  a dictionary with times as keys and bias values as values, 2) a reference_link_end_type, 3) reference_link_end_type, 4) list of reference epochs
     # Let's simulate two arcs
     bias_value_per_start_time = dict()
     bias_value_per_start_time[0] = np.array([1e-2])
     bias_value_per_start_time[60] = np.array([2e-2])
     reference_link_end_type = observation.receiver # set bias at receiving link end
     ref_epochs = [0,60]
     arcwise_time_drift_bias_settings = observation.arc_wise_time_drift_bias(bias_value_per_start_time, reference_link_end_type, ref_epochs)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_time_drift_bias_settings)




     )doc" );

    m.def( "time_bias", &tom::constantTimeBias, py::arg( "time_bias" ), py::arg( "associated_link_end" ) );

    m.def( "arcwise_time_bias",
           py::overload_cast< const std::map< double, double >&, const tom::LinkEndType >( &tom::arcWiseTimeBias ),
           py::arg( "time_bias_per_arc_start_time" ),
           py::arg( "associated_link_end" ) );

    m.def( "combined_bias",
           &tom::multipleObservationBiasSettings,
           py::arg( "bias_list" ),
           R"doc(

 Function for creating settings for a combined observation bias.

 Function for creating settings for a combined observation bias, calculating by combining any number of bias types.
 Each contribution of the combined bias is computed from the unbiased observable, so when applying both a relative and absolute bias, we get:

 .. math::
    \tilde{h}=h+K_{a}+hK_{r}

 And, crucially:

 .. math::
    \tilde{h}\neq (h+K_{a})(1+K_{r})

 where :math:`K_{r}` and :math:`K_{a}` is the relative and absolute bias, respectively.


 Parameters
 ----------
 bias_list : List[:class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings`]
     A list containing the bias settings that are to be applied to the observable.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.multipleObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.multipleObservationBiasSettings` class, combining the settings for multiple observation biases.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the combined_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function combined_bias() allows to combine multiple ObservationBiasSettings objects.
     # Let's combine absolute, relative and time_drift biases.
     bias_array = np.array([1e-2])

     # Define absolute and relative bias settings
     absolute_bias_settings = observation.absolute_bias(bias_array)
     relative_bias_settings = observation.absolute_bias(bias_array)

     # Define Time Drift Bias
     time_link_end = observation.receiver
     ref_epoch = 0
     time_drift_bias_settings = observation.time_drift_bias(bias_array, time_link_end, ref_epoch)

     # combined_bias takes a list of ObservationBiasSettings objects as input
     bias_list = [absolute_bias_settings, relative_bias_settings, time_drift_bias_settings]
     combined_bias = observation.combined_bias(bias_list)

     # Show that it returns an ObservationBiasSettings object.
     print(combined_bias)




     )doc" );

    m.def( "two_way_time_scale_range_bias", &tom::twoWayTimeScaleRangeBias, R"doc(No documentation found.)doc" );

}

}
}
}
}