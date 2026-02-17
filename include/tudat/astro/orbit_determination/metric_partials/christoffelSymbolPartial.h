/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
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

//! Base interface for Christoffel-symbol partial derivatives.
class ChristoffelSymbolPartial
{
public:
    //! Constructor.
    /*!
     *  \param acceleratedBody Name of accelerated body for which Christoffel partials are evaluated.
     */
    explicit ChristoffelSymbolPartial( const std::string& acceleratedBody ):
        acceleratedBody_( acceleratedBody ) { }

    virtual ~ChristoffelSymbolPartial( ) = default;

    //! Evaluate partials w.r.t. state of an integrated body.
    /*!
     *  \param stateReferencePoint Pair containing body and optional reference point identifier.
     *  \param integratedStateType Integrated state type for which partial is requested.
     *  \return Vector of Christoffel-symbol partial matrices.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > wrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return getDerivativeFunctionWrtStateOfIntegratedBody( stateReferencePoint, integratedStateType ).first( );
    }

    //! Retrieve function that computes partials w.r.t. state of an integrated body.
    /*!
     *  \param stateReferencePoint Pair containing body and optional reference point identifier.
     *  \param integratedStateType Integrated state type for which partial is requested.
     *  \return Pair of (partial function, number of output columns).
     */
    virtual std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    //! Evaluate partials w.r.t. scalar estimatable parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Christoffel-symbol partial matrices; empty when unavailable.
     */
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

    //! Retrieve function for partials w.r.t. scalar parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size). Default implementation returns an empty function.
     */
    virtual std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Evaluate partials w.r.t. vector estimatable parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Christoffel-symbol partial matrices; empty when unavailable.
     */
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

    //! Retrieve function for partials w.r.t. vector parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Pair of (partial function, parameter size). Default implementation returns an empty function.
     */
    virtual std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Update cached quantities used in partial evaluation.
    /*!
     *  \param currentTime Current evaluation time.
     */
    virtual void update( const double currentTime ) = 0;

protected:
    std::string acceleratedBody_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_CHRISTOFFELSYMBOLPARTIAL_H
