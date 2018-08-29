/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_ORBITDETERMINATION_H
#define TUDAT_JSONINTERFACE_ORBITDETERMINATION_H

#include "Tudat/Astrodynamics/OrbitDetermination/podInputOutputTypes.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `EstimatableParameterSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< PodSettings >& parameterSettings );

//! Create a shared pointer to a `EstimatableParameterSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< PodSettings >& parameterSettings );

void updateInverseAPrioriCovarianceFromJSON(
        const nlohmann::json& jsonObject, const int numberOfParameters, Eigen::MatrixXd& inverseAprioriCovariance );

void updateObservationWeightsFromJSON(
        const nlohmann::json& jsonObject,
        const std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, int > > numberOfObservations,
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >& observableWeights );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ORBITDETERMINATION_H
