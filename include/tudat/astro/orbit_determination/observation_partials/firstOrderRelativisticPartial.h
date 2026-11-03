/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FIRSTORDERRELATIVISTICLIGHTTIMECORRECTIONPARTIAL_H
#define TUDAT_FIRSTORDERRELATIVISTICLIGHTTIMECORRECTIONPARTIAL_H

#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
/*!
 * Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
 * \param singleBodyLightTimeCorrection Value of light-time correction for which partial is to be evaluated.
 * \param bodyGravitationalParameter Value of gravitational parameter
 * \return Partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
 */
double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter( const double singleBodyLightTimeCorrection,
                                                                                             const double bodyGravitationalParameter );

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
/*!
 * Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
 * \param totalLightTimeCorrection Value of light-time correction for which partial is to be evaluated.
 * \param ppnParameterGamma Value of PPN parameter gamma.
 * \return Partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
 */
double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtPpnParameterGamma( const double totalLightTimeCorrection,
                                                                                  const double ppnParameterGamma );

//! Class for computing the partial derivatives of 1st-order relativistic light-time correction.
class FirstOrderRelativisticLightTimeCorrectionPartial : public LightTimeCorrectionPartial
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param correctionCalculator Object used to compute light-time correction.
     */
    FirstOrderRelativisticLightTimeCorrectionPartial(
            const std::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > correctionCalculator ):
        LightTimeCorrectionPartial( observation_models::first_order_relativistic ), correctionCalculator_( correctionCalculator )
    {
        perturbingBodies_ = correctionCalculator_->getPerturbingBodyNames( );
        perturbingBodyGravitationalParameterFunctions_ = correctionCalculator_->getPerturbingBodyGravitationalParameterFunctions( );
    }

    //! Destructor.
    ~FirstOrderRelativisticLightTimeCorrectionPartial( ) {}

    //! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
    /*!
     * Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
     * \param states States of link ends for observation for which partial is to be computed
     * \param times Times of link ends for observation for which partial is to be computed
     * \param bodyIndex Index of body in correctionCalculator_->currentLighTimeCorrectionComponents_ for which partial is to
     * be computed. NOTE: this index must be smaller than the number of bodies perturbing the light-time.
     * \return Partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
     */
    SingleOneWayRangePartialReturnType wrtBodyGravitationalParameter( const std::vector< Eigen::Vector6d >& states,
                                                                      const std::vector< double >& times,
                                                                      const int bodyIndex )
    {
        if( !( static_cast< int >( perturbingBodies_.size( ) ) > bodyIndex ) )
        {
            throw std::runtime_error( "Error, bodyIndex in FirstOrderRelativisticLightTimeCorrectionPartial not consistent with contents" );
        }

        double partialValue = getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter(
                correctionCalculator_->getCurrentLightTimeCorrectionComponent( bodyIndex ),
                perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( ) );
        return std::make_pair( ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partialValue ).finished( ),
                               ( times[ 0 ] + times[ 1 ] ) / 2.0 );
    }

    //! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
    /*!
     * Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
     * \param states States of link ends for observation for which partial is to be computed
     * \param times Times of link ends for observation for which partial is to be computed
     * \return Partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
     */
    SingleOneWayRangePartialReturnType wrtPpnParameterGamma( const std::vector< Eigen::Vector6d >& states,
                                                             const std::vector< double >& times )
    {
        double partialValue = getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtPpnParameterGamma(
                correctionCalculator_->getCurrentTotalLightTimeCorrection( ), correctionCalculator_->getPpnParameterGammaFunction_( )( ) );
        return std::make_pair( ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partialValue ).finished( ),
                               ( times[ 0 ] + times[ 1 ] ) / 2.0 );
    }

    SingleOneWayRangeGradientPartialReturnType gradientWrtPpnParameterGamma( const std::vector< Eigen::Vector6d >& states,
                                                                             const std::vector< double >& times,
                                                                             const observation_models::LinkEndType gradientLinkEnd )
    {
        return std::make_pair( correctionCalculator_->calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                                       states.at( 0 ), states.at( 1 ), times.at( 0 ), times.at( 1 ), gradientLinkEnd, nullptr ) /
                                       ( correctionCalculator_->getPpnParameterGammaFunction_( )( ) + 1.0 ),
                               gradientLinkEnd == observation_models::receiver ? times.at( 1 ) : times.at( 0 ) );
    }

    //! Gamma partial of the explicit link-end time derivative, at fixed link-end states.
    SingleOneWayRangePartialReturnType timeDerivativeWrtPpnParameterGamma( const std::vector< Eigen::Vector6d >& states,
                                                                           const std::vector< double >& times,
                                                                           const observation_models::LinkEndType linkEnd )
    {
        // Use the same first-order gamma dependence as the delay and position-gradient partials.
        const double partial = correctionCalculator_->calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
                                       states.at( 0 ), states.at( 1 ), times.at( 0 ), times.at( 1 ), linkEnd, nullptr ) /
                ( correctionCalculator_->getPpnParameterGammaFunction_( )( ) + 1.0 );
        return std::make_pair( ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partial ).finished( ),
                               linkEnd == observation_models::receiver ? times.at( 1 ) : times.at( 0 ) );
    }

    //! Gravitational-parameter partial of a single body's light-time correction gradient w.r.t. link-end position.
    /*!
     * Uses the first-order mu dependence of the delay (as wrtBodyGravitationalParameter), which is not exact when bending
     * is included.
     * \param states States of link ends for observation for which partial is to be computed
     * \param times Times of link ends for observation for which partial is to be computed
     * \param gradientLinkEnd Link end at which the position gradient is evaluated
     * \param bodyIndex Index of the perturbing body for which partial is to be computed
     * \return Partial derivative of the position gradient w.r.t. the body's gravitational parameter, and associated time
     */
    SingleOneWayRangeGradientPartialReturnType gradientWrtBodyGravitationalParameter( const std::vector< Eigen::Vector6d >& states,
                                                                                      const std::vector< double >& times,
                                                                                      const observation_models::LinkEndType gradientLinkEnd,
                                                                                      const int bodyIndex )
    {
        return std::make_pair( correctionCalculator_->calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                                       states.at( 0 ), states.at( 1 ), times.at( 0 ), times.at( 1 ), gradientLinkEnd, bodyIndex ) /
                                       perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( ),
                               gradientLinkEnd == observation_models::receiver ? times.at( 1 ) : times.at( 0 ) );
    }

    //! Gravitational-parameter partial of a single body's explicit link-end time derivative, at fixed link-end states.
    /*!
     * Uses the first-order mu dependence of the delay (as wrtBodyGravitationalParameter), which is not exact when bending
     * is included.
     * \param states States of link ends for observation for which partial is to be computed
     * \param times Times of link ends for observation for which partial is to be computed
     * \param linkEnd Link end whose time is differentiated
     * \param bodyIndex Index of the perturbing body for which partial is to be computed
     * \return Partial derivative of the time derivative w.r.t. the body's gravitational parameter, and associated time
     */
    SingleOneWayRangePartialReturnType timeDerivativeWrtBodyGravitationalParameter( const std::vector< Eigen::Vector6d >& states,
                                                                                    const std::vector< double >& times,
                                                                                    const observation_models::LinkEndType linkEnd,
                                                                                    const int bodyIndex )
    {
        const double partial = correctionCalculator_->calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
                                       states.at( 0 ), states.at( 1 ), times.at( 0 ), times.at( 1 ), linkEnd, bodyIndex ) /
                perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( );
        return std::make_pair( ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partial ).finished( ),
                               linkEnd == observation_models::receiver ? times.at( 1 ) : times.at( 0 ) );
    }

    //! Function to get the names of bodies causing light-time correction.
    /*!
     * Function to get the names of bodies causing light-time correction.
     * \return Names of bodies causing light-time correction.
     */
    std::vector< std::string > getPerturbingBodies( )
    {
        return perturbingBodies_;
    }

protected:
    //! Object used to compute light-time correction.
    std::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > correctionCalculator_;

    //! Names of bodies causing light-time correction.
    std::vector< std::string > perturbingBodies_;

    //! Set of functions returning the gravitational parameters of the gravitating bodies.
    std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // FIRSTORDERRELATIVISTICLIGHTTIMECORRECTIONPARTIAL_H
