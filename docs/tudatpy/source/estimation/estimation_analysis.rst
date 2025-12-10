.. _estimation_analysis:

``estimation_analysis``
========================

This module contains the top-level functionality for performing a covariance analysis or state/parameter estimation in Tudat. The class through which these procedures are performed is the :class:`~tudatpy.estimation.estimation_analysis.Estimator` class, with dedicacted classes for input and output for both a covariance analysis and estimation (:class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput`, :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput`, :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`, :class:`~tudatpy.estimation.estimation_analysis.EstimationOutput`). The procedure for state estimation is described in detail on a `series of user guide pages <https://docs.tudat.space/en/latest/user-guide/state-estimation.html>`_, with the final step of using the ``Estimator`` class described `here <https://docs.tudat.space/en/latest/user-guide/state-estimation/estimation-settings.html>`_.

In addition to functionality to perform and supprt the estimation, this modue also contains functionality to propagate the covariance that is obtained at the initial epoch :math:`t_{0}` to any other epoch (see :func:`tudatpy.estimation.estimation_analysis.propagate_covariance`).

Finally, this module also contains functions of convenience to perform the estimation of state and parameters based on pseudo-observations using the :func:`~tudatpy.estimation.estimation_analysis.create_best_fit_to_ephemeris` function.

Functions
---------
.. currentmodule:: tudatpy.estimation.estimation_analysis

.. autosummary::

   propagate_covariance

   propagate_covariance_from_analysis_objects

   propagate_formal_errors

   estimation_convergence_checker

   create_best_fit_to_ephemeris



.. autofunction:: tudatpy.estimation.estimation_analysis.propagate_covariance

.. autofunction:: tudatpy.estimation.estimation_analysis.propagate_covariance_from_analysis_objects

.. autofunction:: tudatpy.estimation.estimation_analysis.propagate_formal_errors

.. autofunction:: tudatpy.estimation.estimation_analysis.estimation_convergence_checker

.. autofunction:: tudatpy.estimation.estimation_analysis.create_best_fit_to_ephemeris





Classes
-------
.. currentmodule:: tudatpy.estimation.estimation_analysis

.. autosummary::

   CovarianceAnalysisInput

   EstimationInput

   CovarianceAnalysisOutput

   EstimationOutput

   Estimator

   EstimationConvergenceChecker



.. autoclass:: tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput
   :members:
   :special-members: __init__



.. autoclass:: tudatpy.estimation.estimation_analysis.EstimationInput
   :members:
   :special-members: __init__



.. autoclass:: tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.EstimationOutput
   :members:


.. autoclass:: tudatpy.estimation.estimation_analysis.Estimator
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker
   :members:

