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

Eigen::Matrix< double, 1, 3 > getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
        double currentBodyDistance,
        double externalBodyGravitationalParameter,
        const Eigen::Vector3d& externalBodyRelativePosition );

Eigen::Matrix< double, 1, 3 > getFirstOrderBarycentricToBodyCentricTimeDerivativeWrtBodyPosition(
        double externalBodyGravitationalParameter,
        const Eigen::Vector3d& externalBodyRelativePosition );

class FirstOrderBarycentricToBodyCentricTimeDerivativePartial : public RelativisticTimeDerivativePartial
{
public:
    explicit FirstOrderBarycentricToBodyCentricTimeDerivativePartial(
            const std::shared_ptr< FirstOrderBarycentricToBodyCentricTimeStateDerivative< double, double > >
            stateDerivativeModel );

    ~FirstOrderBarycentricToBodyCentricTimeDerivativePartial( ) override = default;

    Eigen::Matrix< double, 1, 6 > wrtTranslationalStateOfBody( const int bodyIndex )
    {
        return ( Eigen::Matrix< double, 1, 6 >( )
                 << currentPartialsWrtExternalBodyPositions_.at( bodyIndex ),
                 Eigen::Matrix< double, 1, 3 >::Zero( ) ).finished( );
    }

    std::function< Eigen::Matrix< double, 1, 6 >( ) > getDerivativeFunctionWrtTranslationalStateOfBody(
            const std::string& bodyName ) override;

    bool isStateDerivativePartialWrtTranslationalStateNonNull(
            const std::string& bodyName ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    Eigen::Matrix< double, 1, 1 > getPartialWrtBodyGravitationalParameter( int bodyIndex )
    {
        return ( Eigen::Matrix< double, 1, 1 >( )
                 << -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT /
                    stateDerivativeModel_->getCurrentExternalBodyDistance( bodyIndex ) ).finished( );
    }

    Eigen::Matrix< double, 1, 3 > getPartialWrtBodyPosition( int bodyIndex )
    {
        return currentPartialsWrtExternalBodyPositions_.at( bodyIndex );
    }

    Eigen::Matrix< double, 1, 3 > getPartialWrtCentralBodyPosition( )
    {
        return currentPartialWrtCentralBodyPosition_;
    }

    Eigen::Matrix< double, 1, 3 > getPartialWrtCentralBodyVelocity( )
    {
        return -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                stateDerivativeModel_->getCurrentCentralBodyState( ).segment( 3, 3 ).transpose( );
    }

    Eigen::Matrix< double, 1, 6 > getPartialWrtCentralBodyState( )
    {
        return ( Eigen::Matrix< double, 1, 6 >( )
                 << getPartialWrtCentralBodyPosition( ), getPartialWrtCentralBodyVelocity( ) ).finished( );
    }

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
