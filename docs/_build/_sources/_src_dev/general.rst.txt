
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

PEP 20 -- The Zen of Python
------------------------------------

    >>> import this

    | The Zen of Python, by Tim Peters
    |
    | Beautiful is better than ugly.
    | Explicit is better than implicit.
    | Simple is better than complex.
    | Complex is better than complicated.
    | Flat is better than nested.
    | Sparse is better than dense.
    | Readability counts.
    | Special cases aren't special enough to break the rules.
    | Although practicality beats purity.
    | Errors should never pass silently.
    | Unless explicitly silenced.
    | In the face of ambiguity, refuse the temptation to guess.
    | There should be one-- and preferably only one --obvious way to do it.
    | Although that way may not be obvious at first unless you're Dutch.
    | Now is better than never.
    | Although never is often better than *right* now.
    | If the implementation is hard to explain, it's a bad idea.
    | If the implementation is easy to explain, it may be a good idea.
    | Namespaces are one honking great idea -- let's do more of those!


PEP 8 -- Style Guide for Python Code
------------------------------------

Adapted from
`PEP 8`_


.. _`PEP 8`: https://www.python.org/dev/peps/pep-0008/