**Added:**

* Documentation for :module:`tudatpy.plotting` module added, with example plots.
* Copy to build tree for ``docs`` folder during CMake configuration.

**Changed:**

* ``tudatpy::spice_interface`` changed to ``tudatpy::spice``.
* ``tudatpy/docs`` changed to comply with implementation of :module:`multidoc`.
* ``tudatpy/include/tudatpy/docstrings.h`` to comply with :module:`multidoc`.

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

**Fixed:**

* Fixed how version is handled from ``tudatpy/version``, so that labels
  are not passed to CMake project version.
* All example applications in ``tudatpy/examples`` for latest API following
  ``feature/po_updates`` development.

**Security:**

* <news item>
