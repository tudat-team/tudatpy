==================
TudatPy Change Log
==================

.. current developments

v0.6.2
====================

**Authors:**




v0.6.1
====================

**Added:**

* Documentation for :module:`tudatpy.plotting` module added, with example plots.
* Copy to build tree for ``docs`` folder during CMake configuration.
* Rever tooling support.

**Changed:**

* ``tudatpy::spice_interface`` changed to ``tudatpy::spice``.
* ``tudatpy/docs`` changed to comply with implementation of :module:`multidoc`.
* ``tudatpy/include/tudatpy/docstrings.h`` to comply with :module:`multidoc`.
* Removed ``.conda`` build recipe from source.
* Removed ``wheel-setup.py``, as PyPi is not currently a goal.

**Deprecated:**

* Deprecated `tudatpy.kernel``, ``kernel`` layer has been collapsed. All
  code still using ``kernel`` will function. No deprecation warning has been
  added yet.
* Deprecated ``tudatpy.interface.spice_interface``, changed to
  :module:`tudatpy.interface.spice`. All code still using ``spice_interface``
  will function. No deprecation warning has been added yet.

**Removed:**

* ``tudatpy/.conda`` directory.
* ``tudatpy/include/docstrings`` directory (compliance with ``multidoc``
  in ``tudatpy/include/tudatpy/docstrings.h``).
* Removed test and tutorials from main ``tudatpy`` source.
* Removed old dev python layering ``_layer_propagation_setup.py``.
* Removed ``tudatpy/tools``.

**Fixed:**

* Fixed how version is handled from ``tudatpy/version``, so that labels
  are not passed to CMake project version.
* All example applications in ``tudatpy/examples`` for latest API following
  ``feature/po_updates`` development.

**Authors:**

* Geoffrey H. Garrett
* Dominic Dirkx
* foggionni
* Jonas Hener
* gaffarelj
* MiguelAvillez
* kimonito98
* transferorbit
* Filippo Oggionni
* jgaffarel



v0.5.27
====================

**Authors:**




v0.5.26
====================

**Authors:**




v0.5.25
====================

**Authors:**




v0.5.24
====================

**Authors:**

* Geoffrey H. Garrett



v0.5.23-1
====================

**Authors:**

* Geoffrey H. Garrett



v0.5.23
====================

**Authors:**

* Geoffrey H. Garrett



v0.5.23-rc1
====================

**Authors:**




v0.5.23-rc
====================

**Added:**

* Added ``rever`` support for the repository.
    - Manually configured aliases for Dominic Dirkx, Jonas Hener and Geoffrey
      Garrett.

**Authors:**

* Geoffrey H. Garrett
* Dominic Dirkx
* Elmar Puts
* Jonas Hener
* The Gitter Badger


