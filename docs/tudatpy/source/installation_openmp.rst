Optional parallel spherical-harmonic gravity
============================================

The spherical-harmonic acceleration factory supports an optional native OpenMP
backend. This parallelizes Legendre calculations and independent acceleration
terms within one propagation. It does not parallelize integration steps or
Hessian summation, and does not require free-threaded Python.

Building the backend
--------------------

Follow the normal source-build instructions and install the usual Tudat
prerequisites first. A C++ compiler and runtime supporting OpenMP are additionally
required. Configure a separate build directory with:

.. code-block:: bash

   cmake -S . -B build-openmp -DTUDAT_BUILD_WITH_OPENMP=ON
   cmake --build build-openmp --target kernel --parallel 4

Include the compiler, dependency-prefix and other options required by your
normal source build. CMake discovers the OpenMP C++ target and reports an error
if the requested backend is unavailable. The option defaults to ``OFF``.
Changing it requires rebuilding the native libraries and Python kernel
consistently; do not combine new headers with older library objects. Installed
CMake targets export their OpenMP dependency when enabled.

A prebuilt Python package only supports this feature if its maintainer built it
with this option. Setting an environment variable cannot enable a backend that
was not compiled into that package.

Using gravity threads
---------------------

.. code-block:: python

   from tudatpy.dynamics import propagation_setup

   acceleration = propagation_setup.acceleration
   print(acceleration.is_openmp_available())
   gravity = acceleration.spherical_harmonic_gravity(
       maximum_degree=100,
       maximum_order=100,
       number_of_threads=4,
   )
   print(gravity.number_of_threads)

``number_of_threads`` is keyword-only and defaults to one, even in an OpenMP
build. Existing calls, including the third positional ``remove_point_mass``
argument, retain their meaning. Zero or negative counts raise ``ValueError``.
Requesting multiple threads in a build without the backend raises
``RuntimeError``. The requested count is an upper bound; runtime limits may
produce a smaller team. Settings report the requested count, not the actual team.

Each native acceleration model owns its cache and reusable term buffer. A model
must not be evaluated concurrently by unrelated callers. Separate models can be
used independently. The environment is updated on the calling thread before
workers read it; inner workers do not invoke Python callbacks. This option does
not change the caller's GIL handling. Exceptions are delivered on the calling
thread after workers finish. Scientific archives and pickle retain their
existing schema: the thread policy is not stored, and loaded settings use one
thread. Construct new settings to select another policy.

Calls inside an active OpenMP region use sequential inner execution. This does
not detect outer workers implemented with ``std::thread`` or separate Python
processes. When combining independent arcs with inner gravity threads, budget
the product of outer workers and inner threads explicitly. More threads are not
necessarily faster; saved-per-term force output and near-pole force evaluation
retain their sequential path.

For controlled measurements, prevent unrelated BLAS/OpenMP work from adding
another thread layer, for example with ``OPENBLAS_NUM_THREADS=1`` and
``OMP_NUM_THREADS=1``. The gravity model's explicit setting still controls its
own team. Do not treat a particular machine's speedup as a portable guarantee.
