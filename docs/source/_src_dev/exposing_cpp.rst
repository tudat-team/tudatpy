
Exposing C++ in Python
======================

Modules
-------

core
~~~~

.. code-block:: cpp

    #include <pybind11/pybind11.h>

    int add(int i, int j) {
        return i + j;
    }

    PYBIND11_MODULE(core, m) {
        m.doc() = "Core exposed module from Tudat C++ libraries"; // optional module docstring

        m.def("add", &add, "A function which adds two numbers");
    }

Submodules
~~~~~~~~~~

Functions
-------------------

General Functions
~~~~~~~~~~~~~~~~~~

Let's say we want to expose :class:`~tudat::simulation_setup::setGlobalFrameBodyEphemerides` in the Python
environment. Since we wish our namespaces to be modules in Python, we would go to ``expose_simulation_setup.cpp`` under
the source ``tudatpy/tudatpy/``.


.. code-block:: cpp

     template< typename StateScalarType = double, typename TimeType = double >
     void setGlobalFrameBodyEphemerides( const NamedBodyMap& bodyMap,
                                         const std::string& globalFrameOrigin,
                                         const std::string& globalFrameOrientation )

Overloaded Functions
~~~~~~~~~~~~~~~~~~~~

Classes
-------

Overloaded Constructors
~~~~~~~~~~~~~~~~~~~~~~~

Inheritance
~~~~~~~~~~~

Virtual Members and Inheritance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Documented Errors
-----------------

TypeError: __init__()
~~~~~~~~~~~~~~~~~~~~~~

Not so evident case
'''''''''''''''''''

Attempting to instantiate an exposed object in Python:

.. code-block:: python

    instantiation = tudatpy.object_X(object_Y, 25.0E3, True, False, None)

Raises a TypeError:

.. code-block:: none
    :emphasize-lines: 6

    Traceback (most recent call last):

    ...

    TypeError: __init__(): incompatible constructor arguments. The following argument types are supported:
        1. tudatpy.object_X(arg1: tudatpy.object_Y, arg2: float, arg3: bool, arg4: bool = False, arg5: tudat::object_Z = None)

    Invoked with: kwargs: arg1=<tudatpy.object_Y object at 0x7f27040042f0>, arg2=25000.0, arg3=True, arg4=False, arg5=None

    Process finished with exit code 1

It may not be evident at first, but when going into the source code of the exposed class, something may look amiss:

.. code-block:: c++
    :emphasize-lines: 17

    py::class_<
            tudat::object_X,
            std::shared_ptr<tudat::object_X>,
            tudat::parent_of_object_X
    >(m, "object_X")
            .def(py::init<
                         const std::shared_ptr<tudat::object_Y>,
                         const double,
                         const bool,
                         const bool,
                         const std::shared_ptr<tudat::object_Z>
                 >(),
                 py::arg("arg1"),
                 py::arg("arg2"),
                 py::arg("arg3"),
                 py::arg("arg4") = false,
                 py::arg("arg5") = nullptr
            );


We can notice that the supported signature has an unexposed tudat object: ``arg5``. It's evident that it's unexposed, as
``arg1`` shows clearly that a ``tudatpy`` object is an accepted input. It's not expected that an error should be raised
as the default argument is set as a ``nullptr`` in the C++ module definition, which is itself a supported type conversion
in Pybind11. If an argument has a default associated to it, the type must be defined somewhere in the Pybind11 module,
prior to the definition of the relevant constructor.

