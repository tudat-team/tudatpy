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

#include <tudat/simulation/environment/body.h>
#include <tudat/astro/mission_segments/createTransferTrajectory.h>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::mission_segments;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;

//! Test function for a new interplanetary trajectory class in Tudat
struct MultipleGravityAssist
{

    MultipleGravityAssist( ){ }

    MultipleGravityAssist(
            const SystemOfBodies& bodyMap,
            const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
            const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
            const std::vector< std::string >& nodeIds,
            const std::string& centralBody,
            const std::vector< std::vector< double > > problemBounds_ ):
        bodyMap_( bodyMap ), legSettings_( legSettings ), nodeSettings_( nodeSettings ),
        nodeIds_( nodeIds ), centralBody_( centralBody ), problemBounds_( problemBounds_ )
    {
        numberOfNodes_ = nodeIds_.size( );

        getParameterVectorDecompositionIndices(
               legSettings_,  nodeSettings_, legParameterIndices_, nodeParameterIndices_ );
    }

    void getDecomposedDecisionVector(
            const Eigen::VectorXd rawDecisionVariables,
            std::vector< double >& currentNodeTimes,
            std::vector< Eigen::VectorXd >& currentLegFreeParameters,
            std::vector< Eigen::VectorXd >& currentNodeFreeParameters ) const;

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    std::vector< double >::size_type get_nobj() const
    {
        return 1u;
    }

private:

    SystemOfBodies bodyMap_;
    std::vector< std::shared_ptr< TransferLegSettings > > legSettings_;
    std::vector< std::shared_ptr< TransferNodeSettings > > nodeSettings_;
    std::vector< std::string > nodeIds_;
    std::string centralBody_;
    const std::vector< std::vector< double > > problemBounds_;

    std::vector< std::pair< int, int > > nodeParameterIndices_;
    std::vector< std::pair< int, int > > legParameterIndices_;

    unsigned int numberOfNodes_;

    mutable std::shared_ptr< TransferTrajectory > transferTrajectory_;
};

#endif // TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H
