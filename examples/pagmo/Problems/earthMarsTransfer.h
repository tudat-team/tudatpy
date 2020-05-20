/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
#define TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H

#include <vector>
#include <utility>
#include <limits>

#include <Eigen/Core>

#include <tudat/astro/basic/orbitalElementConversions.h>
#include <tudat/astro/basic/convertMeanToEccentricAnomalies.h>
#include <tudat/astro/mission_segments/multiRevolutionLambertTargeterIzzo.h>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"

/*!
 *  The class defined in this file is to be used in a Pagmo optimization. It defines the objective function for an Earth-Mars
 *  two-burn impulsive transfer, using a Lambert targeter. The independent variables are:
 *
 *  1) Departure Julian day
 *  2) Time-of-flight of Earth-Mars transfer
 *
 *  The problem minimized the Delta V, and if requested, minimizes the time of flight
 */

using namespace pagmo;

//! Test function for a new interplanetary trajectory class in Tudat
struct EarthMarsTransfer
{

    typedef Eigen::Matrix< double, 6, 1 > StateType;

    //! Default constructor, required for Pagmo compatibility
    EarthMarsTransfer( ): useTripTime_( false ){ }

    //! Constructor that sets boundaries of independent variables, and a boolean denoting whether the fitness is single-objective
    //! (Delta V), or dual objective (Delta V and time of flight).
    EarthMarsTransfer( std::vector< std::vector< double > > &bounds, const bool useTripTime = false );

    //! Calculate the fitness as a function of the parameter vector x
    std::vector< double > fitness( const std::vector< double > &x ) const;

    //! Retrieve the allowable limits of the parameter vector x: pair containing minima and maxima of parameter values
    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    //! Retrieve the name of the problem
    std::string get_name( ) const;

    //! Serialization function for Pagmo compatibility
    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    //! Retrieve the number of objectives in problem, e.g. the size of the vector returned by the fitness function
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

    StateType getPlanetPosition( const double date, const std::string planetName ) const;

    bool useTripTime_;
};

#endif // TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
