/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INITIALMASSSTATE_H
#define TUDAT_INITIALMASSSTATE_H

#include <boost/function.hpp>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of an initial rotational state.
template< typename InitialStateParameterType = double >
class InitialMassStateParameter: public EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:


    InitialMassStateParameter(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialMassState ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >(
            initial_mass_state, associatedBody ),
        initialMassState_( initialMassState ){ }


    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialMassState_;
    }

    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  parameterValue )
    {
        initialMassState_ = parameterValue;
    }


    int getParameterSize( )
    {
        return 1;
    }


private:

    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialMassState_;
};


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_INITIALMASSSTATE_H
