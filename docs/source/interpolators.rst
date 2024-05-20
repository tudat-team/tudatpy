``interpolators``
=================
Functions for performing numerical interpolation using various algorithms. More details on the usage of these functions is given in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/https://docs.tudat.space/en/latest/_src_user_guide/mathematics/interpolators.html>`_












Functions
---------
.. currentmodule:: tudatpy.math.interpolators

.. autosummary::

   linear_interpolation

   piecewise_constant_interpolation

   cubic_spline_interpolation

   hermite_spline_interpolation

   lagrange_interpolation

   create_one_dimensional_scalar_interpolator

   create_one_dimensional_vector_interpolator

   create_one_dimensional_matrix_interpolator



.. autofunction:: tudatpy.math.interpolators.linear_interpolation

.. autofunction:: tudatpy.math.interpolators.piecewise_constant_interpolation

.. autofunction:: tudatpy.math.interpolators.cubic_spline_interpolation

.. autofunction:: tudatpy.math.interpolators.hermite_spline_interpolation

.. autofunction:: tudatpy.math.interpolators.lagrange_interpolation

.. autofunction:: tudatpy.math.interpolators.create_one_dimensional_scalar_interpolator

.. autofunction:: tudatpy.math.interpolators.create_one_dimensional_vector_interpolator

.. autofunction:: tudatpy.math.interpolators.create_one_dimensional_matrix_interpolator




Enumerations
------------
.. currentmodule:: tudatpy.math.interpolators

.. autosummary::

   BoundaryInterpolationType

   AvailableLookupScheme

   LagrangeInterpolatorBoundaryHandling



.. autoclass:: tudatpy.math.interpolators.BoundaryInterpolationType
   :members:

.. autoclass:: tudatpy.math.interpolators.AvailableLookupScheme
   :members:

.. autoclass:: tudatpy.math.interpolators.LagrangeInterpolatorBoundaryHandling
   :members:




Classes
-------
.. currentmodule:: tudatpy.math.interpolators

.. autosummary::

   InterpolatorSettings

   LagrangeInterpolatorSettings

   OneDimensionalInterpolatorScalar

   OneDimensionalInterpolatorVector

   OneDimensionalInterpolatorMatrix



.. autoclass:: tudatpy.math.interpolators.InterpolatorSettings
   :members:

.. autoclass:: tudatpy.math.interpolators.LagrangeInterpolatorSettings
   :members:

.. autoclass:: tudatpy.math.interpolators.OneDimensionalInterpolatorScalar
   :members:

.. autoclass:: tudatpy.math.interpolators.OneDimensionalInterpolatorVector
   :members:

.. autoclass:: tudatpy.math.interpolators.OneDimensionalInterpolatorMatrix
   :members:



