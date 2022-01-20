
class ClassNotDeclaredError(Exception):
    """Class not declared in API Declaration.

    """
    ...


class MethodNotDeclaredError(Exception):
    """Method not declared in API Declaration.

    Examples
    ========

    Example 1
    ---------

    .. code-block:: yaml
       :caption: ``example-1.yaml``

        classes:
          - name: ExampleClass
            methods:
              - name: constructor
              - name: foo_method
              # there is no bar_method !

    .. code-block:: C++
       :caption: ``ExampleClass.hpp``
       :emphasize-lines: 13

       #include <string>

       //! @get_docstring(ExampleClass.__docstring__)
       class ExampleClass {
       public:

          //! @get_docstring(ExampleClass.constructor)
          ExampleClass();

          //! @get_docstring(ExampleClass.foo_method)
          std::string foo_method();

          //! @get_docstring(ExampleClass.bar_method)
          std::string bar_method();

       }

    >>> from multidoc.parsing import parse_api_declaration
    >>> from multidoc.testing import TESTING_DIR
    >>> import os
    >>> s = parse_api_declaration(os.path.join(TESTING_DIR, "example-1.yaml"))

    """
    ...


class FunctionNotDeclaredError(Exception):
    """Function not declared in API Declaration.

    """
    ...


class OverloadNotFoundError(Exception):
    """Overload not declared in API Declaration.

    """
    ...
