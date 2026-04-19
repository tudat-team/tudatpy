/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          T. Moyer (2000), Formulation for Observed and Computed Values of Deep Space Network Data
 *              Types for Navigation, DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 *          B. Nava, P. Coisson, and S.M. Radicella (2008), A new version of the NeQuick ionosphere
 *              electron density model, J. Atmos. Sol.-Terr. Phys., 70(15), doi:10.1016/j.jastp.2008.01.015
 *          M. Hernandez-Pajares et al. (2009), The IGS VTEC maps: a reliable source of ionospheric
 *              information since 1998, J. Geodesy, 83(3-4), doi:10.1007/s00190-008-0266-1
 */

#ifndef TUDAT_NEQUICK2IONOSPHERICCORRECTION_H
#define TUDAT_NEQUICK2IONOSPHERICCORRECTION_H

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/astro/basic_astro/neQuick2Model.h"

#include <functional>
#include <Eigen/Core>

namespace tudat
{
namespace observation_models
{

//! Path-integrated ionospheric correction using NeQuick-2 electron density model.
/*!
 * Computes the ionospheric light-time correction by numerically integrating the
 * NeQuick-2 electron density along the transmitter-receiver ray path using
 * Gauss-Legendre quadrature. Optionally constrained by IONEX VTEC maps.
 *
 * This correction works for arbitrary geometries: ground-to-space, space-to-space,
 * and embedded-in-ionosphere configurations (e.g., ISS at ~420 km).
 *
 * delay = sign * 40.3 * STEC / f^2 / c
 * where STEC is the path-integrated electron density (el/m^2).
 */
class NeQuick2IonosphericCorrection : public LightTimeCorrection
{
public:

    NeQuick2IonosphericCorrection(
        std::shared_ptr< environment::NeQuick2Model > neQuick2Model,
        std::shared_ptr< environment::IonexConstrainedNeQuick2Model > rescaledModel,
        const ObservableType observableType,
        std::function< Eigen::Vector6d( double ) > earthStateFunction,
        std::function< Eigen::Matrix3d( double ) > earthRotationToBodyFixedFunction,
        double earthEquatorialRadius,
        double firstOrderDelayCoefficient = 40.3,
        int quadratureOrder = 50 );

    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings ) override;

private:

    std::shared_ptr< environment::NeQuick2Model > neQuick2Model_;
    std::shared_ptr< environment::IonexConstrainedNeQuick2Model > rescaledModel_;

    std::function< Eigen::Vector6d( double ) > earthStateFunction_;
    std::function< Eigen::Matrix3d( double ) > earthRotationToBodyFixedFunction_;

    double earthEquatorialRadius_;
    double firstOrderDelayCoefficient_;
    int quadratureOrder_;
    int sign_;
};

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_NEQUICK2IONOSPHERICCORRECTION_H
