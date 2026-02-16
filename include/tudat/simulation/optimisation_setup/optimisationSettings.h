/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OPTIMISATION_SETTINGS_SETUP_H
#define TUDAT_OPTIMISATION_SETTINGS_SETUP_H

#include <utility>
#include <vector>

#include "pagmo/algorithm.hpp"

namespace tudat
{

namespace simulation_setup
{

//! Class defining settings for optimisation.
/*!
 *  Class defining settings for optimisation.
 */
class OptimisationSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param lowThrustLegType Type of low-thrust leg that is to be used.
     */
    OptimisationSettings(
            const pagmo::algorithm optimisationAlgorithm,
            const int numberOfGenerations,
            const int numberOfIndividualsPerPopulation,
            const double relativeToleranceConstraints = 1.0e-6,
            const std::pair< std::vector< double >, double > initialGuessThrustModel = std::make_pair( std::vector< double >( ), 0.0 ) ):
        optimisationAlgorithm_( optimisationAlgorithm ), numberOfGenerations_( numberOfGenerations ),
        numberOfIndividualsPerPopulation_( numberOfIndividualsPerPopulation ),
        relativeToleranceConstraints_( relativeToleranceConstraints ), initialGuessThrustModel_( initialGuessThrustModel )
    { }

    //! Destructor.
    virtual ~OptimisationSettings( ) { }

    //! Optimisation algorithm.
    pagmo::algorithm optimisationAlgorithm_;

    //! Number of generations.
    int numberOfGenerations_;

    //! Number of individuals per population.
    int numberOfIndividualsPerPopulation_;

    //! Relative tolerance for optimisation constraints.
    double relativeToleranceConstraints_;

    //! Initial guess for the optimisation.
    std::pair< std::vector< double >, double > initialGuessThrustModel_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OPTIMISATION_SETTINGS_SETUP_H
