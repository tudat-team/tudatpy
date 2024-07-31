``root_finders``
================
Functions for root-finding of algebraic univariate functions using various algorithms. 
In this submodule, various methods can be selected to:

.. math::
   \text{find }x\text{, such that: } f(x)=0

The methods that can be selected here are iterative methods to find approximate solutions to the above problem,
starting from an initial guess :math:`x_{0}`, for which :math:`f_{0}=f(x_{0})`. The various methods implement different
iterative algorithms to compute :math:`x_{i}\rightarrow x_{i+1}`. That is, they (attempt to) compute an improved guess of the root
at iteration :math:`i+1` from the guess at iteration :math:`i`, and continue iterating until convergence has been reached.

Depending on the method that is use, the root-finder may need several initial guesses, or have a formulation for one or
more derivatives of the function :math:`f(x)`. If the required information is not available when performing the root-finding,
and exception will be thrown.

There are several convergence criteria that can be defined for this

* When an absolute tolerance :math:`\epsilon_{a}` is met, such that :math:`|x_{i}-x_{i-1}|<\epsilon_{a}`
* When a relative tolerance :math:`\epsilon_{r}` is met, such that :math:`|(x_{i}-x_{i-1})/x_{i}|<\epsilon_{r}`
* When the root function gets within :math:`\epsilon_{f}` of the true root :math:`|f(x_{i}|<\epsilon_{f}`)
* When the number if iterations exceeds some threshold :math:`N`, such that :math:`i=N`

The root finder algorithm continues until as single one of the required convergence criteria is met. 

When meeting the convergence criterion on number of iterations :math:`N`, a user can choose to deal with this in one of several manners (see below).












Functions
---------
.. currentmodule:: tudatpy.math.root_finders

.. autosummary::

   bisection

   newton_raphson

   secant

   halley



.. autofunction:: tudatpy.math.root_finders.bisection

.. autofunction:: tudatpy.math.root_finders.newton_raphson

.. autofunction:: tudatpy.math.root_finders.secant

.. autofunction:: tudatpy.math.root_finders.halley




Enumerations
------------
.. currentmodule:: tudatpy.math.root_finders

.. autosummary::

   MaximumIterationHandling



.. autoclass:: tudatpy.math.root_finders.MaximumIterationHandling
   :members:




Classes
-------
.. currentmodule:: tudatpy.math.root_finders

.. autosummary::

   RootFinderSettings



.. autoclass:: tudatpy.math.root_finders.RootFinderSettings
   :members:



