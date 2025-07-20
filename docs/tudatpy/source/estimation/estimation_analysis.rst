.. _estimation_analysis:

``estimation_analysis``
========================
This module contains the functionality for managing the inputs and outputs of an estimation.












Functions
---------
.. currentmodule:: tudatpy.estimation.estimation_analysis

.. autosummary::

   simulate_observations

   compute_target_angles_and_range

   propagate_covariance

   propagate_formal_errors

   estimation_convergence_checker



.. autofunction:: tudatpy.estimation.estimation_analysis.simulate_observations

.. autofunction:: tudatpy.estimation.estimation_analysis.compute_target_angles_and_range

.. autofunction:: tudatpy.estimation.estimation_analysis.propagate_covariance

.. autofunction:: tudatpy.estimation.estimation_analysis.propagate_formal_errors

.. autofunction:: tudatpy.estimation.estimation_analysis.estimation_convergence_checker






Classes
-------
.. currentmodule:: tudatpy.estimation.estimation_analysis

.. autosummary::

   EstimatableParameterSet

   ObservationViabilityCalculator

   ObservationSimulator

   ObservationCollection

   SingleObservationSet

   CombinedStateTransitionAndSensitivityMatrixInterface

   EstimationConvergenceChecker

   CovarianceAnalysisInput

   EstimationInput

   CovarianceAnalysisOutput

   EstimationOutput

   Estimator



.. autoclass:: tudatpy.estimation.estimation_analysis.EstimatableParameterSet
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.ObservationViabilityCalculator
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.ObservationSimulator
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.ObservationCollection
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.SingleObservationSet
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.CombinedStateTransitionAndSensitivityMatrixInterface
   :members:

.. autoclass:: tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker
   :members:

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


