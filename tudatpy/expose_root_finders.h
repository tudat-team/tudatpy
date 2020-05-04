//
// Created by ggarrett on 04-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_ROOT_FINDERS_H
#define TUDATBUNDLE_EXPOSE_ROOT_FINDERS_H

#include <pybind11/pybind11.h>

#include "Tudat/Mathematics/RootFinders/halleyRootFinder.h"
#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/secantRootFinder.h"
#include "Tudat/Mathematics/RootFinders/bisection.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"
#include "Tudat/Mathematics/RootFinders/createRootFinder.h"

///RootFinders/UnitTests/testFunction.h"
///RootFinders/UnitTests/testFunction1.h"
///RootFinders/UnitTests/testFunction2.h"
///RootFinders/UnitTests/testFunction3.h"
//}/RootFinders/UnitTests/testFunctionWithLargeRootDifference.h"
///RootFinders/UnitTests/testFunctionWithZeroRoot.h"


namespace py=pybind11;
using namespace tudat::root_finders;

namespace tudatpy{
    void expose_root_finders(py::module &m);
}


#endif //TUDATBUNDLE_EXPOSE_ROOT_FINDERS_H
