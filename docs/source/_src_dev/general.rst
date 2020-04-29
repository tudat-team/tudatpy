
General
=======

The primary aim of TudatPy is the expose the C++ libraries of `Tudat`_. As such, the current design pattern of TudatPy
for release 1.0 is to implement the namespaces within the libraries as self contained modules within the TudatPy package.
This means that the `Tudat Doxygen`_ page is *the current* primary resource in the design of the 1.0 release. The
`Tudat namespace list`_ is indicative of the final submodule list for 1.0, save for exampels such as the
``tudat::input_output`` namespace which is replaced by submodules of popular packages such as `pandas`_ and `NumPy`_
which are native to the Python ecosystem. For example uses of these popular Python packages in TudatPy applications,
see the :ref:`python_ecosystem_link`.

.. _Tudat: https://tudat.tudelft.nl/
.. _Tudat Doxygen: http://doxygen.tudat.tudelft.nl/
.. _Tudat namespace list: http://doxygen.tudat.tudelft.nl/namespaces.html
.. _pandas: https://pandas.pydata.org/
.. _NumPy: https://numpy.org/
