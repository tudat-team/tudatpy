/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Basics/timeType.h"

namespace tudat
{

namespace propagators
{



template class DynamicsSimulator< double, double >;
template class SingleArcDynamicsSimulator< double, double >;
template class MultiArcDynamicsSimulator< double, double >;
template class HybridArcDynamicsSimulator< double, double >;

#if( BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class MultiArcDynamicsSimulator< long double, double >;
template class MultiArcDynamicsSimulator< double, Time >;
template class MultiArcDynamicsSimulator< long double, Time >;

template class HybridArcDynamicsSimulator< long double, double >;
template class HybridArcDynamicsSimulator< double, Time >;
template class HybridArcDynamicsSimulator< long double, Time >;

template class SingleArcDynamicsSimulator< long double, double >;
template class SingleArcDynamicsSimulator< double, Time >;
template class SingleArcDynamicsSimulator< long double, Time >;

template class DynamicsSimulator< long double, double >;
template class DynamicsSimulator< double, Time >;
template class DynamicsSimulator< long double, Time >;
#endif

}

}

