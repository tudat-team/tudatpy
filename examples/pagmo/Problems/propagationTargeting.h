/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_EXAMPLE_PAGMO_PROBLEM_PROPAGATION_TARGETING_HPP
#define TUDAT_EXAMPLE_PAGMO_PROBLEM_PROPAGATION_TARGETING_HPP

#include <tudat/simulation/simulation.h>
#include <utility>
#include <vector>

using namespace tudat;

// Define the problem PaGMO-style
struct PropagationTargetingProblem {

    // Empty constructor
    PropagationTargetingProblem( ){ }

    PropagationTargetingProblem( const double altitudeOfPerigee, const double altitudeOfApogee,
                                 const double altitudeOfTarget, const double longitudeOfTarget,
                                 const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave,
                                 const bool useExtendedDynamics = false );

    // Fitness: takes the value of the RAAN and returns the value of the closest distance from target
    std::vector<double> fitness(const std::vector<double> &x) const;

    // Boundaries of the problem set between 0 and (360) degrees
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;

    Eigen::VectorXd getPreviousFinalState( )
    {
        return previousFinalState_;
    }

    std::map< double, Eigen::VectorXd > getPreviousStateHistory( )
    {
        return previousStateHistory_;
    }

    std::map< double, Eigen::VectorXd > getPreviousDependentVariablesHistory()
    {
        return previousDependentVariablesHistory_;
    }

    Eigen::VectorXd getPreviousDependentVariablesFinalValues()
    {
        return previousDependentVariablesFinalValues_;
    }



private:

    double altitudeOfPerigee_;
    double altitudeOfApogee_;
    double altitudeOfTarget_;
    double longitudeOfTarget_;
    double earthRadius_;
    double radiusOfPerigee_;
    double radiusOfApogee_;
    double earthGravitationalParameter_;
    double semiMajorAxis_;
    double simulationStartEpoch_;
    double simulationEndEpoch_;

    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave_;

    bool useExtendedDynamics_;

    mutable std::map< double, Eigen::VectorXd > previousStateHistory_;
    mutable Eigen::VectorXd previousFinalState_;
    mutable tudat::simulation_setup::NamedBodyMap bodyMap_;
    mutable std::map< double, Eigen::VectorXd > previousDependentVariablesHistory_;
    mutable Eigen::VectorXd previousDependentVariablesFinalValues_;




};

#endif // TUDAT_EXAMPLE_PAGMO_PROBLEM_PROPAGATION_TARGETING_HPP
