/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_WEIGHT_SETTINGS_H
#define TUDAT_OBSERVATION_WEIGHT_SETTINGS_H

#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace observation_models
{

//! Settings describing observation weights before an ObservationDataset is created or extended.
struct ObservationWeightSettings {
    //! Supported ways to represent the observation weights.
    enum class WeightsBlockType {
        default_weights,           //!< Unit weight for every scalar observation component.
        constant_scalar,           //!< One scalar weight used for every observation and component.
        scalar_per_observation,    //!< One scalar weight per observation, repeated over its components.
        diagonal_per_observation,  //!< A separate component-level diagonal for every observation.
        constant_block,            //!< One component-level matrix used for every observation.
        block_per_observation,     //!< A separate component-level matrix for every observation.
        set_block                  //!< One matrix covering every component of every observation in the set.
    };

    //! Create settings that assign unit diagonal weights.
    static ObservationWeightSettings defaultWeights( )
    {
        return ObservationWeightSettings( );
    }

    //! Create settings for one scalar weight repeated over all observations and components.
    static ObservationWeightSettings constantScalar( const double weight )
    {
        ObservationWeightSettings settings;
        settings.type_ = WeightsBlockType::constant_scalar;
        settings.scalarWeight_ = weight;
        return settings;
    }

    //! Create settings with one scalar weight per observation event.
    static ObservationWeightSettings scalarPerObservation( const std::vector< double >& weights )
    {
        ObservationWeightSettings settings;
        settings.type_ = WeightsBlockType::scalar_per_observation;
        settings.scalarWeights_ = weights;
        return settings;
    }

    //! Create settings with one component-level diagonal per observation event.
    static ObservationWeightSettings diagonalPerObservation( const std::vector< Eigen::VectorXd >& weights )
    {
        ObservationWeightSettings settings;
        settings.type_ = WeightsBlockType::diagonal_per_observation;
        settings.diagonalWeights_ = weights;
        return settings;
    }

    //! Create settings with one component block repeated for every observation event.
    static ObservationWeightSettings constantBlock( const Eigen::MatrixXd& weightBlock )
    {
        ObservationWeightSettings settings;
        settings.type_ = WeightsBlockType::constant_block;
        settings.weightBlock_ = weightBlock;
        return settings;
    }

    //! Create settings with a separate component block for every observation event.
    static ObservationWeightSettings blockPerObservation( const std::vector< Eigen::MatrixXd >& weightBlocks )
    {
        ObservationWeightSettings settings;
        settings.type_ = WeightsBlockType::block_per_observation;
        settings.weightBlocks_ = weightBlocks;
        return settings;
    }

    //! Create settings with one full principal block for the complete set.
    static ObservationWeightSettings setBlock( const Eigen::MatrixXd& weightBlock )
    {
        ObservationWeightSettings settings;
        settings.type_ = WeightsBlockType::set_block;
        settings.weightBlock_ = weightBlock;
        return settings;
    }

    //! Selected weight representation.
    WeightsBlockType type_ = WeightsBlockType::default_weights;
    //! Scalar value used by the constant-scalar representation.
    double scalarWeight_ = 1.0;
    //! Per-observation values used by the scalar-per-observation representation.
    std::vector< double > scalarWeights_;
    //! Per-observation component diagonals used by the diagonal-per-observation representation.
    std::vector< Eigen::VectorXd > diagonalWeights_;
    //! Shared component block or complete set block, depending on type_.
    Eigen::MatrixXd weightBlock_;
    //! Per-observation component blocks used by the block-per-observation representation.
    std::vector< Eigen::MatrixXd > weightBlocks_;
};

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHT_SETTINGS_H
