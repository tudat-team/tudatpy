// TODO: Disable warnings according to detected compiler. GCC/clang.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wpragmas"

#include <Python.h> // Must be first
#include <vector>
#include <boost/python.hpp>


#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "conversion.h"
#include "eigen_numpy.h"
#include "typedefs.h"

#include "boost/shared_ptr.hpp"
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"

namespace bp = boost::python;

// Test function for std::vector<int>.
std::vector<int> intsToTen() {
    std::vector<int> result;
    result.push_back(1);
    result.push_back(2);
    result.push_back(3);
    result.push_back(4);
    result.push_back(5);
    result.push_back(6);
    result.push_back(7);
    result.push_back(8);
    result.push_back(9);
    result.push_back(10);
    return result;
}

// Test function for std::vector<std::string> of Solar System.
std::vector<std::string> stringVectorSolarSystem() {
    std::vector<std::string> result;
    result.push_back("Mercury");
    result.push_back("Venus");
    result.push_back("Earth");
    result.push_back("Mars");
    result.push_back("Jupiter");
    result.push_back("Saturn");
    result.push_back("Uranus");
    result.push_back("Neptune");
    result.push_back("Pluto");
    return result;
};

// Test function for conversion of std::vector<std::string> to Python list of strings.
bp::list pythonWrapperForStringVectorSolarSystem() {
    return tudatpy::conversion::std_vector_to_py_list(stringVectorSolarSystem());
}

// Test function for conversion of std::vector<int> to Python list of ints.
bp::list pythonWrapperForIntsToTen() {
    return tudatpy::conversion::std_vector_to_py_list(intsToTen());
}

BOOST_PYTHON_MODULE (test_module) {
    bp::def("python_wrapper_for_string_vector_solar_system", &pythonWrapperForStringVectorSolarSystem);
    bp::def("python_wrapper_for_ints_to_ten", &pythonWrapperForIntsToTen);
}

#pragma GCC diagnostic pop