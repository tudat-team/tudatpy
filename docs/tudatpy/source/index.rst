API Reference
=============

The TU Delft Astrodynamics Toolbox (Tudat) is a powerful set of libraries that support astrodynamics and space research and education. It has been used for a wide variety of purposes, ranging from the study of reentry dynamics to interplanetary mission orbit estimation. The core functionality of Tudat is implemented in C++ and exposed to Python. Tudat(Py) is entirely open-source (under BSD 3-clause license) and disseminated using conda package. Installation instructions can be found `here <https://docs.tudat.space/en/latest/getting-started/installation.html>`_.

To get started with using Tudat, we recommend starting on our user guide where a top-level overview of core functionality on `state propagation <https://docs.tudat.space/en/latest/user-guide/state-propagation.html>`_ and `state estimation https://docs.tudat.space/en/latest/user-guide/state-estimation.html <https://docs.tudat.space/en/latest/user-guide/state-estimation.html>`_ can be found. A list of worked out `example scripts <https://docs.tudat.space/en/latest/index-examples.html>`_, ranging from simple starting points to detailed analyses is also provided

On this page, we provide a comprehensive overview of the functionality available throug the Tudat Python interface. The functionalitv is divided into the following submodules


* :ref:`dynamics`: This submodule contains the interfaces for one of core application of Tudatpy: numerical state propagation. The functionality in this submodule consists of a large number of interconnected elements that work together as a whole. It facilitates applications ranging from low-fidelity orbit modelling in the context of a global optimization to high-fidelity modelling for applications in precise orbit determination. The framework supports single- an multi-arc propagation, translational dynamics, rotational dynamics, as well as combinations of these.
* :ref:`estimation`: This submodule contains the interfaces for the other core application of Tudatpy: numerical state and parameter estimation. The functionality in this module builds on the ``dynamics`` module, and allows for the estimation of initial states, environmental parameters, observation model parameters from both simulated data and real observations. The framework supports range, Doppler and angular position data among others, and has interfaces load data from ESTRACK, DSN and MPC.
* :ref:`astro`: This submodule contains various (semi-)standalone functions for astrodynamics applications, which can be used very well outside of a Tudat application/propagation. Submodules contain lists of frame conversions, element conversion, elementary orbit calculations, *etc.*.
* :ref:`math`:  This submodule contains various functions and classes for purely mathematical operations, such as numerical interpolation, quadrature *etc.*.
* :ref:`trajectory_design`: This submodule contains (basic) functionality for the preliminary design of a full (transfer) orbit, using for instance a Multiple Gravity Assist (MGA) or a low-thrust system. It relies on functionality in the ``astro`` submodule. It is largely independent of the ``dynamics`` and ``estimation`` submodule, but does contain interface functions to allow the preliminary design to be used as an initial guess for a full numerical propagation.
* :ref:`interface`: This submodule contains functionality to interface with various external packages which Tudat uses, such as `SPICE <https://naif.jpl.nasa.gov/naif/toolkit.html>`_
* :ref:`data`: This submodule contains various pieces of functionality for file input-output in Tudatpy, including interfaces to retrieve data from `JPL Horizons <https://ssd.jpl.nasa.gov/horizons/>`_ and the `Minor Planet Center <https://www.minorplanetcenter.net/>`_ among others. Unlike most of the main Tudatpy submodules (which are written in C++, and exposed to Python), this submodule is written partially in Python
* :ref:`plotting`: This submodule contains various pieces of functionality to support the easy plotting of results generated with Tudatpy. Unlike most of the main Tudatpy submodules (which are written in C++, and exposed to Python), this submodule is written in Python
* :ref:`util`: This submodule contains various small pieces of functionality to support the easy post-processing of results generated with Tudatpy. Unlike most of the main Tudatpy submodules (which are written in C++, and exposed to Python), this submodule is written in Python
* :ref:`constants`: This submodule contains various numerical constants used inside Tudat


.. toctree::
   :maxdepth: 2
   :caption: Modules

   astro
   constants
   data
   dynamics
   estimation
   exceptions
   interface
   math
   trajectory_design
   plotting
   util




Bibliography
------------

.. bibliography::
   