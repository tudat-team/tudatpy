/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_HODOGRAPHIC_SHAPING_OPTIMISATION_SETUP_H
#define TUDAT_HODOGRAPHIC_SHAPING_OPTIMISATION_SETUP_H

#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <tudat/astro/basic_astro/physicalConstants.h>
#include <tudat/basics/testMacros.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include "tudat/astro/basic_astro/unitConversions.h"
#include <tudat/astro/basic_astro/orbitalElementConversions.h>
#include <tudat/io/basicInputOutput.h>

#include <random>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"
#include <pagmo/rng.hpp>
#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace pagmo;

namespace tudat
{

namespace shape_based_methods
{

//! Test function for a new low-thrust trajectory class in Tudat
struct FixedTimeHodographicShapingOptimisationProblem
{

    FixedTimeHodographicShapingOptimisationProblem( ){ }

    FixedTimeHodographicShapingOptimisationProblem(
            const Eigen::Vector6d& initialState,
            const Eigen::Vector6d& finalState,
            const double timeOfFlight,
            const double centralBodyGravitationalParameter,
            const int numberOfRevolutions,
            const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            const std::vector< std::vector< double > >& freeCoefficientsBounds ):
        initialState_( initialState ),
        finalState_( finalState ),
        timeOfFlight_( timeOfFlight ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        numberOfRevolutions_( numberOfRevolutions ),
        radialVelocityFunctionComponents_( radialVelocityFunctionComponents ),
        normalVelocityFunctionComponents_( normalVelocityFunctionComponents ),
        axialVelocityFunctionComponents_( axialVelocityFunctionComponents ),
        problemBounds_( freeCoefficientsBounds )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

protected:

private:

    Eigen::Vector6d initialState_;

    Eigen::Vector6d finalState_;

    double timeOfFlight_;

    double centralBodyGravitationalParameter_;

    int numberOfRevolutions_;

    const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents_;

    const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents_;

    const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents_;

    const std::vector< std::vector< double > > problemBounds_;

};


//! Test function for a new low-thrust trajectory class in Tudat
struct HodographicShapingOptimisationProblem
{
    typedef std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > BaseFunctionVector;

    HodographicShapingOptimisationProblem( ){ }

    HodographicShapingOptimisationProblem(
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double ) >& finalStateFunction,
            const double centralBodyGravitationalParameter,
            const int numberOfRevolutions,
            const std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            const bool minimizeMaximumThrust = false,
            const double initialMass = TUDAT_NAN ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        numberOfRevolutions_( numberOfRevolutions ),
        basisFunctionsFunction_( basisFunctionsFunction ),
        problemBounds_( freeCoefficientsBounds ),
        minimizeMaximumThrust_( minimizeMaximumThrust ),
        initialMass_( initialMass )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return minimizeMaximumThrust_ ? 2u : 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double ) > finalStateFunction_;

    double centralBodyGravitationalParameter_;

    int numberOfRevolutions_;

    std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction_;

    const std::vector< std::vector< double > > problemBounds_;

    bool minimizeMaximumThrust_;

    double initialMass_;

};

} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_HODOGRAPHIC_SHAPING_OPTIMISATION_SETUP_H
