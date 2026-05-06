// serialization/pybind11_helpers.h
#pragma once
#include <pybind11/pybind11.h>
#include "tudat/io/serialization/base.h"

namespace py = pybind11;

namespace tudat::serialization {

// For plain value types (non-polymorphic)
template<typename T>
auto make_pickle() {
    return py::pickle(
        [](const T& obj) {
            return py::bytes(serializeToBinaryString(obj));
        },
        [](py::bytes data) {
            return deserializeFromBinaryString<T>(data.cast<std::string>());
        }
    );
}

template<typename Base, typename Derived = Base>
auto make_pickle_polymorphic() {
    return py::pickle(
        [](const std::shared_ptr<Derived>& obj) {
            // upcast to Base so cereal sees the polymorphic pointer
            return py::bytes(serializeSharedPtrToBinaryString<Base>(
                std::static_pointer_cast<Base>(obj)));
        },
        [](py::bytes data) {
            // cereal reconstructs the real derived type inside the shared_ptr<Base>
            // pybind11 then resolves the correct Python class from the dynamic type
            return std::static_pointer_cast<Derived>(
                deserializeSharedPtrFromBinaryString<Base>(
                    data.cast<std::string>()));
        }
    );
}

} // namespace tudat::serialization