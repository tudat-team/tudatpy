.. _migration-guide:

******************************
Migration Guides
******************************

TudatPy v1.0
-----------------

This section provides a guide for users migrating from TudatPy v0.9 and older to v1.0.

Have questions or feedback? As always, let us know in our `Github Discussion forum <https://github.com/orgs/tudat-team/discussions?discussions_q=>`_!

Installation
^^^^^^^^^^^^

How to obtain the new version: the installation of tudatpy for users has not changed, you can simply use our installation using the conda packages
as described on our `installation guide <https://docs.tudat.space/en/latest/getting-started/installation.html#installing-tudat-py>`_. The only difference is that installing tudatpy will now no longer result in a separate ``tudat`` and ``tudatpy`` conda packages being installed, since the two have been merged into a single ``tudatpy`` package.

New module structure
^^^^^^^^^^^^^^^^^^^^

As of v1.0, we have thoroughly re-organized the module structure of tudatpy. We have implemented this such that old import statements will still work (so v1.0 is backwards compatible with code written for older version), but will phase out this support in later versions. A comprehensive overview of the changes to the module structure is given below.

Why was this needed?
====================

There are two driving reasons behind our reorganization of the module structure. The primary reason is the poor organization of functionality in the ``numerical_simulation.estimation`` and ``numerical_simulation.estimation_setup`` modules, which led to a large amount of functionality being combined into a small number of modules
(such as ``numerical_simulation.estimation_setup.observation``), making tudatpy functionality poorly findable, and hampering transparent further developments. An additional reason was that the bulk of the functionality was in the ``numerical_simulation`` module. It was decided to split this into two top-level modules ``dynamics`` and ``estimation``, where the former contains functionality for propagation of dynamics, and the latter contains all functionality required for estimation of states and parameters from (tracking) data.

In addition to these two drivers, we also took the opportunity to re-organize some other minor aspects of the modules, moving some classes and functions into a more suitable module.


How does it effect users?
=========================

The impact on existing scripts should be none, since all reorganization has been done in a backwards compatible manner. If you find any functionality fot which this does not seem to be the case, let us know in our `Github Discussion forum <https://github.com/orgs/tudat-team/discussions?discussions_q=>`_.

Using the pre-v1.0 imports will result in a deprecation warning being printed. We expect to drop support for the older imports after two years, in the course of 2027.

All newly developed functionality will only be provided in the new module structure. In addition, all updates to documentation will only be available in the new module setup, so we recommend all users to update their code to the new import statements.

How to migrate?
====================

Migrating to the new module setup can be done by modifying top-level import statements of modules, and module specifications in the code when using a specific piece of functionality. Most of the propagation-related functionality in the ``environment``, ``environment_setup``, ``propagation`` and ``propagation_setup`` submodules is simply moved from the old ``numerical_simulation`` module to the new ``dynamics`` module. Estimation-related functionality in the ``estimation`` and ``estimation_setup`` submodules has been refactored in multiple submodules, which requires more detailed treatment. Below we provide an overview of all changes that have been made. Implementing these modifications in your code should fully migrate you to v1.0!

+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Old module                                            | New module                                                                                                                                                                |
+=======================================================+===========================================================================================================================================================================+
| ``astro.time_conversion``                             | ``astro.time_representation``                                                                                                                                             |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Top-level ``numerical_simulation``                    | Split in ``dynamics.simulator`` and ``estimation.estimation_analysis``                                                                                                    |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.environment``                  | ``dynamics.environment``                                                                                                                                                  |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.environment_setup``            | ``dynamics.environment_setup``                                                                                                                                            |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.propagation``                  | ``dynamics.propagation``                                                                                                                                                  |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.propagation_setup``            | ``dynamics.propagation_setup``                                                                                                                                            |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.estimation``                   | Distributed in ``estimation.estimation_analysis`` and ``estimation.observations.observations.geometry``                                                                   |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.estimation_setup``             | Distributed in ``dynamics.parameters``, ``dynamics.parameters_setup``, ``estimation.observations_setup.observations_simulation_settings``                                 |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.estimation_setup.parameter``   | Split in ``dynamics.parameters`` and ``dynamics.parameters_setup``                                                                                                        |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``numerical_simulation.estimation_setup.observation`` | Split in ``estimation.observable_models``, ``estimation.observable_models_setup``, ``estimation.observations`` and ``estimation.observations_setup`` and their submodules |
+-------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

New ``Time`` type
^^^^^^^^^^^^^^^^^

Up until v0.9, tudatpy used ``float`` variables to denote both epochs and durations of time. As of v1.0, we have moved to a new setup where time is represented by a dedicated :class:`~tudatpy.astro.time_conversion.Time` class. Also in v0.9 and earlier, it was possible to use the ``Time`` type internally, but this required manual recompilation with specific settings to trigger this behaviour. As of v1.0, we have chosen to make this the default behaviour in our packages.


Why was this needed?
====================

A ``float`` variable has a numerical resolution of about :math:`2\cdot 10^{-16}`, meaning that a relative change below this level cannot be represented. In Tudat, we use seconds since epoch J2000 as time representation. Using a ``float`` for this means that for epochs durther away from J2000, the resolution to which time can be represented degrades. For either 1950 or 2050 (about :math:`1.6\cdot 10^{9}` seconds from J2000) this imposes a hard limit of 0.35 microseconds in resolution of time.

There are several concrete examples of cases where this poor timing resolution limits the performance of analyses. For numerical integration with small time steps, rounding errors in the representation of time have been known to lead to confusing results in (for instance) benchmarking. As an additional example, in Doppler data analysis of planetary missions, the observable is computed by the difference of two light times. Due to the limited resolution in representing epochs, using ``Time`` is required to get state-of-the-art performance. This required Tudat to be manually compiled to use this functionality. With the data analysis framework of Tudat taking an ever more prominent place, it has become important to provide this functionality in the 'normal' package.

How does it affect users?
=========================

The modifications we have made to change to a different time representation have all been made in a backwards compatible manner. An implicit conversion between ``Time`` and ``float`` has been implemented, so that any function that requires a ``Time`` object as input can also take a ``float`` as input. In doing so, the value is 'upconverted` to the higher resolution representation, allowing all later computations to be done at the high resolution.

Various output structures in tudat are provided as dictionaries with time as the independent variable (key). By default, the output a user extracts, for instance from the :attr:`~tudatpy.dynamics.propagation.SingleArcSimulationResults.state_history` attribute for the state history from a numerical propagation will (as in v0.9 and earlier) provide this state history with ``float`` as independent variable, since for most post-processing purposes this is both sufficient and more convenient. However, we now also provide the option to retrieve the state history with time at the native resolution in which the propagation was performed by using the :attr:`~tudatpy.dynamics.propagation.SingleArcSimulationResults.state_history_time_object` attribute. A similar setup has been introduced in various other tudat output options, facilitating backwards compatibility, permitting use of the more typical ``float`` for post-processing, and providing access to the native resolution when required. Some more information on Tudat time representations can be found on our page for `using time representations <https://docs.tudat.space/en/latest/user-guide/state-propagation/environment-setup/times-and-dates.html#using-time-representations>`_.

How to migrate?
====================

No action is required to migrate for this modification. All v0.9 interfaces remain valid and are not deprecated. For various applications, it will not be relevant whether the ``float`` or ``Time`` representation is used internally, and inputs and outputs using ``float`` continue to be valid as they were before. Even for applications where the use of the high-accuracy internal time representation improves numerical results, it will often still be sufficient to provide the input and output at the original ``float`` representation.

Merging of ``tudatpy`` repositories and conda packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The previous structure of the project, with separate conda packages and code repositories for tudat (underlying C++ models) and tudatpy (Python exposure), as well as a tudat-bundle repositories for developers to compile both tudat and tudatpy concurrently, has been a source of various complications and inconsistencies. The codebase from the original tudat repository is now included as a subdirectory within the tudatpy repository (with some reorganization), which contains both the C++ source code and the Python bindings. The tudatpy-examples repository still exists, now as a submodule within tudatpy. The tudat-bundle repository has been deprecated, developers now compile the tudatpy repository directly.

The CMake configuration from tudat has been merged into the main CMakeLists.txt of tudatpy, resulting in a unified build system. The tudatpy repository now follows a mirrored structure: each component has its own tudat (for C++) and tudatpy (for Python) subdirectories. In general, the core logic is located in the tudat folders, while Python bindings and Python-only functionality are placed under the tudatpy folders. The build logic is now largely identical for developing, testing and deploying.

Where some files were previously duplicated or mirrored between tudat and tudatpy, they have now been consolidated—typically by keeping the version from tudatpy when the content was identical.

The tudat conda package and tudat-feedstock repository are now longer used with this change. The tudatpy conda package now contains both the underlying C++ models and the Python exposure.


Why was this needed?
====================
The decision to merge the tudat and tudatpy repositories was driven by the need to simplify development, testing, and packaging workflows. Maintaining them separately had become increasingly cumbersome, and the original motivation for the split no longer reflects how the project is used today.

Here are the main reasons behind the merge:

Simpler build process:
Previously, the tudat and tudatpy conda packages had to be built in sequence, which made the setup more fragile and added unnecessary overhead—especially for contributors and automated packaging.

More effective testing:
With separate repositories, the continuous integration pipelines were also split. This made it difficult to test changes to tudatpy until corresponding changes in tudat were merged, since the PRs weren’t linked. A single repo ensures all components can be tested together.

Reduced complexity:
Developers were expected to build from the tudat-bundle repository, even though the actual source code lived in tudat and tudatpy. This indirection often caused confusion, particularly for new contributors. The merge removes this extra layer.

Consistent configuration:
Maintaining separate build systems (CMake and conda feedstock) for two repositories sometimes led to inconsistencies or duplication of effort. A unified repo makes it easier to keep things aligned.

Changing usage patterns:
The repositories were originally split to support C++-only users. However, most users now rely on the Python interface. With the merged setup, C++-only workflows are still fully supported, but there's no longer a strong reason to keep the two codebases apart.

In short, the merge makes Tudat easier to work with, more robust to maintain, and better aligned with how it's actually used by the community.


How does it effect users?
=========================

The repository restructuring introduces a cleaner and more unified layout, but all core functionalities remain unchanged. Most users will continue working with tudatpy as before.

However, developers who were actively working on branches in the old tudat repository will need to migrate their work to the new combined repository. This typically involves rebasing or transplanting their changes into the appropriate location within the new structure (e.g. moving C++ code to the tudat/ subdirectory within tudatpy).

To preserve commit history and ensure smooth integration, we recommend carefully following the steps outlined in the "How to Migrate?" section below. If you encounter any issues or are unsure how to proceed, feel free to reach out to the core development team on our `discussion forum <https://github.com/orgs/tudat-team/discussions/>`_

How to migrate?
====================

For users, simply creating a new conda environment for tudatpy (as per out `installation guide <https://docs.tudat.space/en/latest/getting-started/installation.html#installing-tudat-py>`_) will migrate to the new setup, without any changes on the user side.

Developers without any active development branches on either tudat or tudatpy (pre-v1.0) should clone the new (v1.0) tudatpy ``develop`` branch, and work with this in the exact same manner as they interacted with the old tudat-bundle repository.

Developers with active development branches on either tudat or tudatpy that have diverged from the ``develop`` branch should contact the tudatpy development team. We can assist in migrating your code to the new repository setup.


