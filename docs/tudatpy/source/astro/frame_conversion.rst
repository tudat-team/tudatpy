.. _frame_conversion:

``frame_conversion``
====================
Conversions between different reference frames.


This module provide a variety of functions and classes to convert
between different reference frames. Functionality to convert between
different state representations is provided in the
:ref:`element_conversion` module. Note that the functionality
here may be used independent of the Tudat models in :ref:`dynamics` or :ref:`estimation`. 
For more details on the use of frames in the context of these models, see 
`this page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html>`_






Notes
-----
- All reference frames used should be assumed as right handed:
  :math:`\mathbf{X}\times\mathbf{Y}=\mathbf{Z}`.

- We distinguish between two different uses of the term 'inertial':

  * An *inertial origin*: the frame has a non-accelerating origin. On solar
    system scales, the solar system barycenter (SSB) is the most 'typical' inertial origin. 
  * An *inertial orientation*: the unit axes of the frame are non-rotating
    with respect to the celestial background. This module is concerned 
    primarily with conversions between different orientations.

- Examples of an inertial orientation include J2000 (a.k.a EME2000), as
  well as the SPICE-defined ECLIPJ2000 frame (see `this description <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#predefined-inertial-frames>`_
  on our user guide). The use of the ICRF frame (difference from J2000
  is <0.1 arcseconds) in Tudat is presently limited to the
  :ref:`dynamics` model, see
  :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`.

.. raw:: html

    <object
    data="_static/J2000_.svg"
    type="image/svg+xml"
    class="center invertible">
    </object>




.. References
.. ----------
.. .. [1] Archinal, B.A., Acton, C.H., A’Hearn, M.F. et al. Report of
..        the IAU Working Group on Cartographic Coordinates and
..        Rotational Elements: 2015. Celest Mech Dyn Astr 130, 22
..        (2018). https://doi.org/10.1007/s10569-017-9805-5






Functions
---------
.. currentmodule:: tudatpy.astro.frame_conversion

.. autosummary::

   rotate_state_to_frame
   
   inertial_to_rsw_rotation_matrix

   rsw_to_inertial_rotation_matrix

   inertial_to_tnw_rotation_matrix

   tnw_to_inertial_rotation_matrix

   inertial_to_body_fixed_rotation_matrix

   body_fixed_to_inertial_rotation_matrix


.. autofunction:: rotate_state_to_frame

.. autofunction:: tudatpy.astro.frame_conversion.inertial_to_rsw_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.rsw_to_inertial_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.inertial_to_tnw_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.tnw_to_inertial_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.inertial_to_body_fixed_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.body_fixed_to_inertial_rotation_matrix

.. _terrestrial_frame_transformations:

Terrestrial affine transformations
---------------------------------

Transform Cartesian positions and moving states between terrestrial-frame
realizations using published linearized Helmert parameters. Load a catalogue
through :ref:`terrestrial_frame_data`, then construct a
``TerrestrialFrameTransformationModel``. Each catalogue retains its published
reference hub and source provenance.

Position evaluation accepts a decimal year. State evaluation also requires its
derivative with respect to the input state's seconds. The following example uses
a uniform 365.25-day year coordinate; the model itself performs no calendar or
time-scale conversion:

.. code-block:: python

   from tudatpy.data_input.environment_data import terrestrial_frames as data
   from tudatpy.astro import frame_conversion as frames

   model = frames.TerrestrialFrameTransformationModel(data.official_itrf2020_parameters())
   source_state = [4075578.385, 931852.890, 4801570.154, 0.012345, -0.023456, 0.034567]
   target_state = model.transform_state(
       source_state, "ITRF2020", "ITRF2014", 2025.0, 1.0 / 31557600.0)

``transform_position`` needs only the decimal year. To transform a full station
history with its calendar coordinate and source-rate policy, use
``ground_station.transform_station_history``; see :ref:`ground_station`.

The affine matrix includes scale and linearized signed rotations:

.. math::

   y=T+A x, \qquad \dot y=\dot T+\dot A x+A\dot x.

``inverse`` uses the actual matrix inverse and its derivative. ``compose_after``
includes all translation and matrix derivative terms. Model composition follows
the declared catalogue hub. Equal nonempty realization labels return identity
before catalogue lookup, including unknown identical labels. Aligned catalogue
release updates retain their realization identity and do not create additional
datum transformations.

The existing
:func:`~tudatpy.dynamics.environment_setup.convert_ground_station_state_between_itrf_frames`
retains the IGN ITRF2014 catalogue and its uniform-Julian-year epoch convention.
Station histories instead retain their calendar-decimal-year convention.

.. autoclass:: TerrestrialFrameAffineTransformation
   :members:

.. autoclass:: TerrestrialFrameTransformationModel
   :members:

.. autofunction:: evaluate_linearized_helmert







