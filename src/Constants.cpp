//
// Created by ggarrett on 21-5-19.
//

#include <boost/python.hpp>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

using namespace boost::python;
using namespace tudat::physical_constants;

BOOST_PYTHON_MODULE(constants)
        {
                boost::python::scope().attr("JULIAN_DAY") = 1;
        }

