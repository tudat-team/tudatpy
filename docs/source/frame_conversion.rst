``frame_conversion``
====================
Conversions between different reference frames.


This module provide a variety of functions and classes to convert
between different reference frames. Functionality to convert between
different state representations is provided in the
:ref:`\`\`element_conversion\`\`` module. Note that the functionality
here may be used independent of the Tudat models in :ref:`\`\`numerical_simulation\`\``. 
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
  :ref:`\`\`numerical_simulation\`\`` model, see
  :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.gcrs_to_itrs`.

.. raw:: html

    <object
    data="_static/J2000_.svg"
    type="image/svg+xml"
    class="center invertible">
    </object>




References
----------
.. [1] Archinal, B.A., Acton, C.H., A’Hearn, M.F. et al. Report of
       the IAU Working Group on Cartographic Coordinates and
       Rotational Elements: 2015. Celest Mech Dyn Astr 130, 22
       (2018). https://doi.org/10.1007/s10569-017-9805-5






Functions
---------
.. currentmodule:: tudatpy.astro.frame_conversion

.. autosummary::

   inertial_to_rsw_rotation_matrix

   rsw_to_inertial_rotation_matrix

   inertial_to_tnw_rotation_matrix

   tnw_to_inertial_rotation_matrix

   inertial_to_body_fixed_rotation_matrix

   body_fixed_to_inertial_rotation_matrix



.. autofunction:: tudatpy.astro.frame_conversion.inertial_to_rsw_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.rsw_to_inertial_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.inertial_to_tnw_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.tnw_to_inertial_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.inertial_to_body_fixed_rotation_matrix

.. autofunction:: tudatpy.astro.frame_conversion.body_fixed_to_inertial_rotation_matrix







