**Added:**

* Newly exposed :module:`tudatpy.astro.time_conversion`. Also created associated documentation.

* Newly exposed set of functions in :module:`tudatpy.io` regarding Missile DATCOM interfacing.

* New function ``write_force_and_moment_coefficients_to_files`` in :module:`tudatpy.kernel.io` to create aerodynamic force and moment coefficient files as a function of angle of attack and Mach number.

* New function ``compare_results`` in :module:`tudatpy.util` to compare the state or dependent variable history of two distinct propagations. Also created associated documentation.

* New function ``redirect_std`` in :module:`tudatpy.util` to redirect or mute all terminal messages when enclosed in a given with statement. Also created associated documentation.

* New function ``pareto_optimums`` in :module:`tudatpy.util` to to return a boolean mask stating whether each multi-dimensional point in a set is a Pareto optimum. Also created associated documentation.

**Changed:**

* Improved the documentation of the ``propagation`` module and all of the ``propagation_setup`` submodules, by adding code snippets, and correcting some errors.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* The :module:`tudatpy.plotting` API documentation page is now building.

* The process of building the API documentation page on readthedocs does not raise errors anymore, mainly caused by broken references or misnames functions.

**Security:**

* <news item>
