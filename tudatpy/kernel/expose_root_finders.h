//
// Created by ggarrett on 04-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_ROOT_FINDERS_H
#define TUDATBUNDLE_EXPOSE_ROOT_FINDERS_H

#include <pybind11/pybind11.h>

#include "tudat/math/root_finders/halleyRootFinder.h"
#include "tudat/math/root_finders/newtonRaphson.h"
#include "tudat/math/root_finders/rootFinder.h"
#include "tudat/math/root_finders/secantRootFinder.h"
#include "tudat/math/root_finders/bisection.h"
#include "tudat/math/root_finders/terminationConditions.h"
#include "tudat/math/root_finders/createRootFinder.h"

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
