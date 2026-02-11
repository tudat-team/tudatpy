/*    Copyright (c) 2010-2019,
 *    Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_CHRISTOFFELSYMBOLPARTIAL_H
#define TUDAT_CHRISTOFFELSYMBOLPARTIAL_H

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/propagators/propagationSettings.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

class ChristoffelSymbolPartial
{
public:
    explicit ChristoffelSymbolPartial( const std::string& acceleratedBody ):
        acceleratedBody_( acceleratedBody ) { }

    virtual ~ChristoffelSymbolPartial( ) = default;

    std::vector< Eigen::Matrix< double, 4, 4 > > wrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return getDerivativeFunctionWrtStateOfIntegratedBody( stateReferencePoint, integratedStateType ).first( );
    }

    virtual std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    std::vector< std::vector< Eigen::Matrix4d > > wrtParameter(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >& parameter )
    {
        std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int > partialFunction =
                getParameterPartialFunction( parameter );

        if( partialFunction.second > 0 )
        {
            return partialFunction.first( );
        }

        return std::vector< std::vector< Eigen::Matrix4d > >( );
    }

    virtual std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    std::vector< std::vector< Eigen::Matrix4d > > wrtParameter(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >& parameter )
    {
        std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int > partialFunction =
                getParameterPartialFunction( parameter );

        if( partialFunction.second > 0 )
        {
            return partialFunction.first( );
        }

        return std::vector< std::vector< Eigen::Matrix4d > >( );
    }

    virtual std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    virtual void update( const double currentTime ) = 0;

protected:
    std::string acceleratedBody_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_CHRISTOFFELSYMBOLPARTIAL_H
