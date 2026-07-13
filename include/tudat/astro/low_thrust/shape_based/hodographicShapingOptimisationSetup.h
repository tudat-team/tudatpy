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

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include "tudat/io/serialization/base.h"

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
struct FixedTimeHodographicShapingOptimisationProblem {
    FixedTimeHodographicShapingOptimisationProblem( ) {}

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
        initialState_( initialState ), finalState_( finalState ), timeOfFlight_( timeOfFlight ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ), numberOfRevolutions_( numberOfRevolutions ),
        radialVelocityFunctionComponents_( radialVelocityFunctionComponents ),
        normalVelocityFunctionComponents_( normalVelocityFunctionComponents ),
        axialVelocityFunctionComponents_( axialVelocityFunctionComponents ), problemBounds_( freeCoefficientsBounds )
    {}

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double >& x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    vector_double::size_type get_nobj( ) const
    {
        return 1u;
    }

    TUDAT_DEFINE_BINARY_IO( FixedTimeHodographicShapingOptimisationProblem )

protected:
private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( initialState_ ) );
        ar( CEREAL_NVP( finalState_ ) );
        ar( CEREAL_NVP( timeOfFlight_ ) );
        ar( CEREAL_NVP( centralBodyGravitationalParameter_ ) );
        ar( CEREAL_NVP( numberOfRevolutions_ ) );
        ar( CEREAL_NVP( radialVelocityFunctionComponents_ ) );
        ar( CEREAL_NVP( normalVelocityFunctionComponents_ ) );
        ar( CEREAL_NVP( axialVelocityFunctionComponents_ ) );
        ar( CEREAL_NVP( problemBounds_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( initialState_ ) );
        ar( CEREAL_NVP( finalState_ ) );
        ar( CEREAL_NVP( timeOfFlight_ ) );
        ar( CEREAL_NVP( centralBodyGravitationalParameter_ ) );
        ar( CEREAL_NVP( numberOfRevolutions_ ) );
        ar( CEREAL_NVP( radialVelocityFunctionComponents_ ) );
        ar( CEREAL_NVP( normalVelocityFunctionComponents_ ) );
        ar( CEREAL_NVP( axialVelocityFunctionComponents_ ) );
        ar( CEREAL_NVP( problemBounds_ ) );
    }

    Eigen::Vector6d initialState_;

    Eigen::Vector6d finalState_;

    double timeOfFlight_;

    double centralBodyGravitationalParameter_;

    int numberOfRevolutions_;

    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents_;

    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents_;

    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents_;

    std::vector< std::vector< double > > problemBounds_;

public:
    bool operator==( const FixedTimeHodographicShapingOptimisationProblem& rhs ) const
    {
        return equals( rhs );
    }

    bool operator!=( const FixedTimeHodographicShapingOptimisationProblem& rhs ) const
    {
        return !( *this == rhs );
    }

    //! Equality comparison via equals method
    bool equals( const FixedTimeHodographicShapingOptimisationProblem& rhs ) const
    {
        return initialState_ == rhs.initialState_ && finalState_ == rhs.finalState_ && timeOfFlight_ == rhs.timeOfFlight_ &&
                centralBodyGravitationalParameter_ == rhs.centralBodyGravitationalParameter_ &&
                numberOfRevolutions_ == rhs.numberOfRevolutions_ &&
                radialVelocityFunctionComponents_ == rhs.radialVelocityFunctionComponents_ &&
                normalVelocityFunctionComponents_ == rhs.normalVelocityFunctionComponents_ &&
                axialVelocityFunctionComponents_ == rhs.axialVelocityFunctionComponents_ && problemBounds_ == rhs.problemBounds_;
    }
};

//! Test function for a new low-thrust trajectory class in Tudat
struct HodographicShapingOptimisationProblem {
    typedef std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > BaseFunctionVector;

    HodographicShapingOptimisationProblem( ) {}

    HodographicShapingOptimisationProblem( const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
                                           const std::function< Eigen::Vector6d( const double ) >& finalStateFunction,
                                           const double centralBodyGravitationalParameter,
                                           const int numberOfRevolutions,
                                           const std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction,
                                           const std::vector< std::vector< double > >& freeCoefficientsBounds,
                                           const bool minimizeMaximumThrust = false,
                                           const double initialMass = TUDAT_NAN ):
        initialStateFunction_( initialStateFunction ), finalStateFunction_( finalStateFunction ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ), numberOfRevolutions_( numberOfRevolutions ),
        basisFunctionsFunction_( basisFunctionsFunction ), problemBounds_( freeCoefficientsBounds ),
        minimizeMaximumThrust_( minimizeMaximumThrust ), initialMass_( initialMass )
    {}

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double >& x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    vector_double::size_type get_nobj( ) const
    {
        return minimizeMaximumThrust_ ? 2u : 1u;
    }

    TUDAT_DEFINE_BINARY_IO( HodographicShapingOptimisationProblem )

protected:
private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( centralBodyGravitationalParameter_ ) );
        ar( CEREAL_NVP( numberOfRevolutions_ ) );
        ar( CEREAL_NVP( problemBounds_ ) );
        ar( CEREAL_NVP( minimizeMaximumThrust_ ) );
        ar( CEREAL_NVP( initialMass_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( centralBodyGravitationalParameter_ ) );
        ar( CEREAL_NVP( numberOfRevolutions_ ) );
        ar( CEREAL_NVP( problemBounds_ ) );
        ar( CEREAL_NVP( minimizeMaximumThrust_ ) );
        ar( CEREAL_NVP( initialMass_ ) );
    }

    double centralBodyGravitationalParameter_;

    int numberOfRevolutions_;

    std::vector< std::vector< double > > problemBounds_;

    bool minimizeMaximumThrust_;

    double initialMass_;

    // Non-serializable function members (ignored in save/load, default-constructed after deserialization)
    std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    std::function< Eigen::Vector6d( const double ) > finalStateFunction_;

    std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction_;

public:
    bool operator==( const HodographicShapingOptimisationProblem& rhs ) const
    {
        return equals( rhs );
    }

    bool operator!=( const HodographicShapingOptimisationProblem& rhs ) const
    {
        return !( *this == rhs );
    }

    //! Equality comparison via equals method
    bool equals( const HodographicShapingOptimisationProblem& rhs ) const
    {
        return centralBodyGravitationalParameter_ == rhs.centralBodyGravitationalParameter_ &&
                numberOfRevolutions_ == rhs.numberOfRevolutions_ && problemBounds_ == rhs.problemBounds_ &&
                minimizeMaximumThrust_ == rhs.minimizeMaximumThrust_ && initialMass_ == rhs.initialMass_;
    }
};

}  // namespace shape_based_methods
}  // namespace tudat

#endif  // TUDAT_HODOGRAPHIC_SHAPING_OPTIMISATION_SETUP_H
