/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H
#define TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H


#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <tudat/astro/basic/physicalConstants.h>
#include <tudat/basics/testMacros.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include "tudat/astro/basic/unitConversions.h"
#include <tudat/astro/basic/orbitalElementConversions.h>

#include <tudat/io/basicInputOutput.h>

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/trajectory_design/trajectory.h"
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
using namespace tudat::transfer_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
using namespace tudat;
using namespace pagmo;

//! Test function for a new interplanetary trajectory class in Tudat
struct MultipleGravityAssist
{

    MultipleGravityAssist( const bool useTripTime = false ): useTripTime_( useTripTime ){ }

    MultipleGravityAssist( std::vector< std::vector< double > > &bounds,
                           std::vector< int > flybySequence,
                           const bool useTripTime = false );

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        if(useTripTime_ )
        {
            return 2u;
        }
        else
        {
            return 1u;
        }

    }

private:

    const std::vector< std::vector< double > > problemBounds_;

    bool useTripTime_;

    int numberOfLegs_;
    std::vector< TransferLegType > legTypeVector_;
    std::vector< std::string > bodyNamesVector_;
    std::vector< ephemerides::EphemerisPointer > ephemerisVector_;
    Eigen::VectorXd gravitationalParameterVector_;
    Eigen::VectorXd semiMajorAxes_;
    Eigen::VectorXd eccentricities_;
    Eigen::VectorXd minimumPericenterRadii_;
};

#endif // TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H
