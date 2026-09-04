.. _tracking_data_optical_utilities:

``optical_utilities``
=====================

.. automodule:: tudatpy.data_input.tracking_data.optical_utilities

This submodule contains functionality to load tracking data from optical
astrometry tables that already exist as pandas or astropy tables. The
:func:`read_optical_data` function is the main interface for loading the data
and converting it to objects that Tudat can process further; see also
:ref:`tracking_data`. The optical-data conversion can assign the default
optical weighting scheme of :cite:t:`veres2017` and star-catalog bias
corrections following :cite:t:`eggl2020`. All other functionality in this
module is reserved for better understanding what data is being loaded, and in
some cases manipulating it, before it is processed into Tudat-compatible
objects.

.. currentmodule:: tudatpy.data_input.tracking_data.optical_utilities

.. autofunction:: read_optical_data

Supporting API
--------------

The functions below implement the shared optical-data conversion pipeline used
by the MPC, 80-column and table-based readers. They are useful for validating,
augmenting, filtering or inspecting optical astrometry tables before they are
converted to Tudat-compatible tracking-data objects.

.. autofunction:: read_astropy_optical_data
.. autofunction:: read_pandas_optical_data
.. autofunction:: optical_table_to_tracking_data
.. autofunction:: create_augmented_optical_table
.. autofunction:: filter_augmented_optical_table
.. autofunction:: standardize_optical_dataframe
.. autofunction:: validate_optical_table
.. autofunction:: load_bias_file
.. autofunction:: get_biases_EFCC18
