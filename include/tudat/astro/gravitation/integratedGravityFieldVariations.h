/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_INTEGRATEDGRAVITYFIELDVARIATIONS_H
#define TUDAT_INTEGRATEDGRAVITYFIELDVARIATIONS_H

#include <map>
#include <memory>

#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

namespace gravitation
{

//! Gravity-field variation defined by numerically integrated degree-two coefficients.
/*!
 * The integrated state contains the unnormalised coefficient corrections in the order
 * [ C20, C21, C22, S21, S22 ]. During propagation the current state supplied by the
 * environment updater is used. Outside propagation, a state-history interpolator is used
 * once propagation results have been installed in the environment.
 */
class IntegratedGravityFieldVariations : public GravityFieldVariations
{
public:
    IntegratedGravityFieldVariations( );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( const double time ) override;

    void setCurrentCoefficientCorrections( const Eigen::VectorXd& coefficientCorrections );

    void setCurrentCoefficientCorrectionDerivative( const Eigen::VectorXd& coefficientCorrectionDerivative );

    void setCoefficientCorrectionHistory( const std::map< double, Eigen::Vector5d >& coefficientCorrectionHistory,
                                          const std::shared_ptr< interpolators::InterpolatorSettings >& interpolatorSettings = nullptr );

    Eigen::Vector5d getCoefficientCorrections( const double time ) const;

    Eigen::Vector5d getCoefficientCorrectionDerivative( const double time ) const;

    void setIsBodyInPropagation( const bool isBodyInPropagation );

    bool getIsBodyInPropagation( ) const;

private:
    Eigen::Vector5d currentCoefficientCorrections_;

    Eigen::Vector5d currentCoefficientCorrectionDerivative_;

    std::map< double, Eigen::Vector5d > coefficientCorrectionHistory_;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector5d > > coefficientCorrectionInterpolator_;

    bool isBodyInPropagation_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_INTEGRATEDGRAVITYFIELDVARIATIONS_H
