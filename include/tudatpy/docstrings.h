
#include <string>
#include <exception>

namespace tudatpy {

    static inline std::string get_docstring(std::string name, int variant = 0) {

        if (name == "@docstring_key") {
            // for testing purposes
            return "@docstring_value";

        //! @generate_docstrings

        } else {

#ifdef DEBUG
            return "No documentation.";
#else
            throw std::runtime_error(
                    "Documentation was not found for:" + name + " (variant: " + std::to_string(variant) + ")");
#endif

        }

    }
}

