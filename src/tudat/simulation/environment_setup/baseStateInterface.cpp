/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/baseStateInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameDoubleState( time );
}

#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameLongDoubleState( time );
}
#endif

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameDoubleState( time );
}

#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template<>
Eigen::Matrix< HighPrecisionStateScalar, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameLongDoubleState( time );
}
#endif

}  // namespace simulation_setup

}  // namespace tudat
