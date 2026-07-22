# Serialization dependency architecture

Tudat serialization is split into layers so that ordinary model headers do not parse archive,
stream, pybind11, or unrelated polymorphic-registration implementation code.

1. `core.h` contains archive-independent Cereal declarations and common standard-library adapters.
2. `eigen.h` contains the Eigen adapter and its archive-input validation.
3. `archives.h` selects the supported binary and JSON archive implementations.
4. `file_io_declarations.h` provides declaration-only macros for non-template public classes.
5. `file_io.h` contains archive/file implementations and definition macros. It belongs in the
   owning model library's serialization implementation files and in headers whose open templates
   cannot be defined out of line.
6. `pybind_helpers.h` contains only type-agnostic Python binding helpers.
7. `registrations_*.h` contain lightweight force-registration declarations. The actual type
   registrations are instantiated once per domain in `registrations*.cpp`.

Internal code must include the narrowest layer it needs. `base.h`, `serialization.h`, and
`registrations.h` are compatibility/convenience umbrella headers and must not be introduced into
new model or binding implementation headers.
