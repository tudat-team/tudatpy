.. tudatpy documentation master file, created by
   sphinx-quickstart on Wed Apr 29 16:32:15 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TudatPy
===================================

TU Delft Astrodynamics Toolbox in Python, or **TudatPy**, is a library that primarily exposes the powerful set of C++
libraries, `Tudat`_. TudatPy aims at accelerating the implementation of Tudat simulations,
providing an interface between Tudat and popular machine learning frameworks and establishing a platform to provide
quality education in the field of astrodynamics.

.. _`Tudat` : https://tudat.tudelft.nl/

.. toctree::
   :maxdepth: 2
   :caption: Modules
   :hidden:

   _src_modules/elements
   _src_modules/constants
   _src_modules/spice_interface
   _src_modules/basic_astrodynamics
   _src_modules/numerical_integrators
   _src_modules/propagators
   _src_modules/simulation_setup

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   _src_tutorials/tudat_tutorials
   _src_tutorials/python_ecosystem
   _src_tutorials/machine_learning

.. toctree::
   :maxdepth: 2
   :caption: Developers
   :hidden:

   _src_dev/general
   _src_dev/documentation
   _src_dev/pull_requests
   _src_dev/exposing_cpp
   _src_dev/changelog
   _src_dev/FAQ

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
