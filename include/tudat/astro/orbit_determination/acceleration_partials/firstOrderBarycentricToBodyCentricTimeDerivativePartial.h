/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_FIRSTORDERBARYCENTRICTOBODYCENTICTIMEDERIVATIVEPARTIAL_H
#define TUDAT_FIRSTORDERBARYCENTRICTOBODYCENTICTIMEDERIVATIVEPARTIAL_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/orbit_determination/acceleration_partials/relativisticTimeDerivativePartial.h"
#include "tudat/astro/propagators/relativisticTimeStateDerivative.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Compute partial of first-order barycentric->bodycentric time derivative w.r.t. body position.
/*!
 *  \param currentBodyDistance Distance between central and external body.
 *  \param externalBodyGravitationalParameter Gravitational parameter of external body.
 *  \param externalBodyRelativePosition Relative position of external body w.r.t. central body.
 *  \return Row vector partial w.r.t. body position.
 */
Eigen::Matrix< double, 1, 3 > getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
        double currentBodyDistance,
        double externalBodyGravitationalParameter,
        const Eigen::Vector3d& externalBodyRelativePosition );

//! Compute partial of first-order barycentric->bodycentric time derivative w.r.t. body position.
/*!
 *  \param externalBodyGravitationalParameter Gravitational parameter of external body.
 *  \param externalBodyRelativePosition Relative position of external body w.r.t. central body.
 *  \return Row vector partial w.r.t. body position.
 */
Eigen::Matrix< double, 1, 3 > getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
        double externalBodyGravitationalParameter,
        const Eigen::Vector3d& externalBodyRelativePosition );

//! Partial model for first-order barycentric->bodycentric relativistic time derivative.
class FirstOrderBarycentricToBodyCentricTimeDerivativePartial : public RelativisticTimeDerivativePartial
{
public:
    //! Constructor.
    /*!
     *  \param stateDerivativeModel First-order barycentric->bodycentric state-derivative model.
     */
    explicit FirstOrderBarycentricToBodyCentricTimeDerivativePartial(
            const std::shared_ptr< FirstOrderBarycentricToBodyCentricTimeStateDerivative< double, double > >
            stateDerivativeModel );

    ~FirstOrderBarycentricToBodyCentricTimeDerivativePartial( ) override = default;

    //! Get current partial w.r.t. translational state of one external body.
    /*!
     *  \param bodyIndex Index of external body in model list.
     *  \return 1x6 state partial.
     */
    Eigen::Matrix< double, 1, 6 > wrtTranslationalStateOfBody( const int bodyIndex )
    {
        return ( Eigen::Matrix< double, 1, 6 >( )
                 << currentPartialsWrtExternalBodyPositions_.at( bodyIndex ),
                 Eigen::Matrix< double, 1, 3 >::Zero( ) ).finished( );
    }

    //! Retrieve function for partial w.r.t. translational state of a body.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return Function returning a 1x6 partial matrix.
     */
    std::function< Eigen::Matrix< double, 1, 6 >( ) > getDerivativeFunctionWrtTranslationalStateOfBody(
            const std::string& bodyName ) override;

    //! Check whether translational-state partial for a body is non-zero.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return True if partial is non-zero.
     */
    bool isStateDerivativePartialWrtTranslationalStateNonNull(
            const std::string& bodyName ) override;

    //! Retrieve function for partial w.r.t. scalar parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Retrieve function for partial w.r.t. vector parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Pair of (partial function, parameter size). Always empty for this model.
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Get current partial w.r.t. one body's gravitational parameter.
    /*!
     *  \param bodyIndex Index of external body in model list.
     *  \return 1x1 parameter partial.
     */
    Eigen::Matrix< double, 1, 1 > getPartialWrtBodyGravitationalParameter( int bodyIndex )
    {
        return ( Eigen::Matrix< double, 1, 1 >( )
                 << -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT /
                    stateDerivativeModel_->getCurrentExternalBodyDistance( bodyIndex ) ).finished( );
    }

    //! Get current partial w.r.t. one external-body position.
    /*!
     *  \param bodyIndex Index of external body in model list.
     *  \return 1x3 position partial.
     */
    Eigen::Matrix< double, 1, 3 > getPartialWrtBodyPosition( int bodyIndex )
    {
        return currentPartialsWrtExternalBodyPositions_.at( bodyIndex );
    }

    //! Get current partial w.r.t. central-body position.
    /*!
     *  \return 1x3 position partial.
     */
    Eigen::Matrix< double, 1, 3 > getPartialWrtCentralBodyPosition( )
    {
        return currentPartialWrtCentralBodyPosition_;
    }

    //! Get current partial w.r.t. central-body velocity.
    /*!
     *  \return 1x3 velocity partial.
     */
    Eigen::Matrix< double, 1, 3 > getPartialWrtCentralBodyVelocity( )
    {
        return -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                stateDerivativeModel_->getCurrentCentralBodyState( ).segment( 3, 3 ).transpose( );
    }

    //! Get current partial w.r.t. full central-body state.
    /*!
     *  \return 1x6 state partial.
     */
    Eigen::Matrix< double, 1, 6 > getPartialWrtCentralBodyState( )
    {
        return ( Eigen::Matrix< double, 1, 6 >( )
                 << getPartialWrtCentralBodyPosition( ), getPartialWrtCentralBodyVelocity( ) ).finished( );
    }

    //! Update cached partial quantities.
    /*!
     *  \param currentTime Current evaluation time.
     */
    void update( const double currentTime = TUDAT_NAN ) override;

protected:
    std::vector< Eigen::Matrix< double, 1, 3 > > currentPartialsWrtExternalBodyPositions_;

    Eigen::Matrix< double, 1, 3 > currentPartialWrtCentralBodyPosition_;

    std::vector< std::string > externalBodies_;

    std::map< std::string, int > externalBodyIndexMap_;

    std::shared_ptr< FirstOrderBarycentricToBodyCentricTimeStateDerivative< double, double > >
        stateDerivativeModel_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_FIRSTORDERBARYCENTRICTOBODYCENTICTIMEDERIVATIVEPARTIAL_H
