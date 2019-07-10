The following method is used within the `core.cpp` module of `pygmo` for the modules:

- `problems`
- `algorithms`
- `islands`

`[core.cpp, Line 148]`  External to `BOOST_PYTHON_MODULE` in `core.cpp`:

```cpp
// Exposed pagmo::problem.
std::unique_ptr<bp::class_<pagmo::problem>> problem_ptr;
```

`[core.cpp, Line 434]` 

```cpp
    // Create the problems submodule.
    std::string problems_module_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".problems");
    PyObject *problems_module_ptr = PyImport_AddModule(problems_module_name.c_str());
    if (!problems_module_ptr) {
        pygmo_throw(PyExc_RuntimeError, "error while creating the 'problems' submodule");
    }
    auto problems_module = bp::object(bp::handle<>(bp::borrowed(problems_module_ptr)));
    bp::scope().attr("problems") = problems_module;
```

`[core.cpp, Line 471]` 

```cpp
    // Store the pointers to the classes that can be extended by APs.
    bp::scope().attr("_problem_address") = reinterpret_cast<std::uintptr_t>(&pygmo::problem_ptr);
```