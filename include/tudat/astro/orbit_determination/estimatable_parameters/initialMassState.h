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
class InitialMassStateParameter : public EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:
    InitialMassStateParameter( const std::string& associatedBody,
                               const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialMassState ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >( initial_mass_state, associatedBody ),
        initialMassState_( initialMassState )
    { }

    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        // Retrieve state from propagator settings (to ensure consistency) if link is set
        if( initialStateGetFunction_ != nullptr )
        {
            initialMassState_ = initialStateGetFunction_( );
        }
        return initialMassState_;
    }

    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > parameterValue )
    {
        // Update state in propagator settings (to ensure consistency) if link is set
        if( initialStateSetFunction_ != nullptr )
        {
            initialStateSetFunction_( parameterValue );
        }
        initialMassState_ = parameterValue;
    }

    int getParameterSize( )
    {
        return 1;
    }

    // Add functions to get and set the state from the propagator settings, to ensure the states in propagator settings and parameters are
    // always identical
    void addStateClosureFunctions(
            const std::function< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >( ) > initialStateGetFunction,
            const std::function< void( const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& ) > initialStateSetFunction )
    {
        initialStateGetFunction_ = initialStateGetFunction;
        initialStateSetFunction_ = initialStateSetFunction;

        // Retrieve state from propagator settings (to ensure consistency) if link is set
        if( initialStateGetFunction_ != nullptr )
        {
            initialMassState_ = initialStateGetFunction_( );
        }
    }

private:
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialMassState_;



    std::function< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >( ) > initialStateGetFunction_;

    std::function< void( const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& ) > initialStateSetFunction_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_INITIALMASSSTATE_H
