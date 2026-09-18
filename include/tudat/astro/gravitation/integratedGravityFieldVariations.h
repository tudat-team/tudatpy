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
    //! Create an initially zero degree-two variation that is not in propagation.
    IntegratedGravityFieldVariations( );

    //! Return the current corrections as normalized cosine and sine coefficient blocks.
    /*!
     * \param time Epoch at which the corrections are requested. The epoch is used by the
     * installed history interpolator outside propagation and ignored for the live state during
     * propagation.
     * \return Pair containing the normalized degree-two cosine and sine correction blocks.
     */
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( const double time ) override;

    //! Set the live unnormalised coefficient corrections read from the propagated state.
    /*!
     * \param coefficientCorrections Five corrections ordered as
     * [ C20, C21, C22, S21, S22 ].
     */
    void setCurrentCoefficientCorrections( const Eigen::VectorXd& coefficientCorrections );

    //! Set the live unnormalised coefficient rates computed by the state derivative model.
    /*!
     * \param coefficientCorrectionDerivative Five rates ordered as
     * [ dC20/dt, dC21/dt, dC22/dt, dS21/dt, dS22/dt ].
     */
    void setCurrentCoefficientCorrectionDerivative( const Eigen::VectorXd& coefficientCorrectionDerivative );

    //! Install a completed propagated history for use outside propagation.
    /*!
     * A multi-point history is interpolated using the supplied settings, or linearly when no
     * settings are supplied. A single-point history is treated as a constant correction.
     * \param coefficientCorrectionHistory Map of epochs to unnormalised five-element corrections.
     * \param interpolatorSettings Settings used to interpolate a multi-point history.
     */
    void setCoefficientCorrectionHistory( const std::map< double, Eigen::Vector5d >& coefficientCorrectionHistory,
                                          const std::shared_ptr< interpolators::InterpolatorSettings >& interpolatorSettings = nullptr );

    //! Retrieve the live or installed unnormalised coefficient corrections at an epoch.
    Eigen::Vector5d getCoefficientCorrections( const double time ) const;

    //! Retrieve the live coefficient rates or a secant rate from the installed history.
    Eigen::Vector5d getCoefficientCorrectionDerivative( const double time ) const;

    //! Select live propagated state or post-propagation history as the correction source.
    void setIsBodyInPropagation( const bool isBodyInPropagation );

    //! Return whether coefficient corrections are currently read from the propagated state.
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
