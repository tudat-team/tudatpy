#include <string>

namespace tudatpy {

static inline std::string get_docstring(std::string name, int variant = 0) {

  if (name == "test" && variant == 0) {
    return "test (variant: 0)";
  } else if (name == "test" && variant == 0) {
    return "test (variant: 1)";
  } else if (name == "test") {
    return "test (unknown variant)";
  } else {
    return "No documentation found.";
  }
}

}// namespace tudatpy
