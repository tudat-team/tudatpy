/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
//#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
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
//using namespace tudat::transfer_trajectories;
//using namespace tudat::shape_based_methods;
//using namespace tudat;
using namespace pagmo;

namespace tudat
{

namespace shape_based_methods
{

//! Test function for a new low-thrust trajectory class in Tudat
struct HodographicShapingOptimisationProblem
{

    HodographicShapingOptimisationProblem( ) /*: timeOfFlight_( physical_constants::JULIAN_DAY ), numberOfRevolutions_( 0 )*/{ }

    HodographicShapingOptimisationProblem(
            Eigen::Vector6d initialState,
            Eigen::Vector6d finalState,
            const double timeOfFlight,
            const int numberOfRevolutions,
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            std::vector< std::vector< double > >& freeCoefficientsBounds ):
        initialState_( initialState ),
        finalState_( finalState ),
        timeOfFlight_( timeOfFlight ),
        numberOfRevolutions_( numberOfRevolutions ),
        bodyMap_( bodyMap ),
        bodyToPropagate_( bodyToPropagate ),
        centralBody_( centralBody ),
        radialVelocityFunctionComponents_( radialVelocityFunctionComponents ),
        normalVelocityFunctionComponents_( normalVelocityFunctionComponents ),
        axialVelocityFunctionComponents_( axialVelocityFunctionComponents ),
        problemBounds_( freeCoefficientsBounds )
    {  }

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;
//    {

//        int numberFreeCoefficientsRadialFunction = radialVelocityFunctionComponents_.size( ) - 3;
//        int numberFreeCoefficientsNormalFunction = normalVelocityFunctionComponents_.size( ) - 3;
//        int numberFreeCoefficientsAxialFunction = axialVelocityFunctionComponents_.size( ) - 3;

//        if ( numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction != x.size( ) )
//        {
//            throw std::runtime_error( "Error, size of design variables vector unconsistent with number of base function components"
//                                      "when making a hodographic shaping optimisation problem." );
//        }

//        Eigen::VectorXd freeCoefficientsRadialVelocityFunction( numberFreeCoefficientsRadialFunction );
//        Eigen::VectorXd freeCoefficientsNormalVelocityFunction( numberFreeCoefficientsNormalFunction );
//        Eigen::VectorXd freeCoefficientsAxialVelocityFunction( numberFreeCoefficientsAxialFunction );

//        for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
//        {
//            freeCoefficientsRadialVelocityFunction[ i ] = x[ i ];
//        }
//        for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
//        {
//            freeCoefficientsNormalVelocityFunction[ i ] = x[ i + numberFreeCoefficientsRadialFunction ];
//        }
//        for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
//        {
//            freeCoefficientsAxialVelocityFunction[ i ] = x[ i + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
//        }

//        HodographicShaping hodographicShaping = HodographicShaping( initialState_, finalState_, timeOfFlight_, numberOfRevolutions_,
//                                                                    bodyMap_, bodyToPropagate_, centralBody_, radialVelocityFunctionComponents_,
//                                                                    normalVelocityFunctionComponents_, axialVelocityFunctionComponents_,
//                                                                    freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
//                                                                    freeCoefficientsAxialVelocityFunction );

//        std::vector< double > fitnessVector;
//        fitnessVector.push_back( hodographicShaping.computeDeltaV( ) );

//        return fitnessVector;
//    }

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

//    std::string get_name( ) const;

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

//    mutable simulation_setup::NamedBodyMap bodyMap_;

private:

    Eigen::Vector6d initialState_;
    Eigen::Vector6d finalState_;
    double timeOfFlight_;
    int numberOfRevolutions_;
    mutable simulation_setup::NamedBodyMap bodyMap_;
    std::string bodyToPropagate_;
    std::string centralBody_;
    mutable std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents_;
    mutable std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents_;
    mutable std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents_;

    const std::vector< std::vector< double > > problemBounds_;

};

} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_HODOGRAPHIC_SHAPING_OPTIMISATION_SETUP_H
