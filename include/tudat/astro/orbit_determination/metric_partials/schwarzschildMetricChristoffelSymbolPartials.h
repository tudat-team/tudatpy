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


#ifndef TUDAT_SCHWARZSCHILDMETRICCHRISTOFFELSYMBOLPARTIALS_H
#define TUDAT_SCHWARZSCHILDMETRICCHRISTOFFELSYMBOLPARTIALS_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/metric_partials/christoffelSymbolPartial.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Christoffel-symbol partial model for the harmonic Schwarzschild metric.
class SchwarzschildMetricChristoffelSymbolPartial : public ChristoffelSymbolPartial
{
public:
    //! Constructor.
    /*!
     *  \param metricObject Harmonic Schwarzschild metric object.
     *  \param centralBody Name of central gravitating body in the metric model.
     */
    SchwarzschildMetricChristoffelSymbolPartial(
            const std::shared_ptr< relativity::HarmonicSchwarzschildMetric > metricObject,
            const std::string& centralBody );

    ~SchwarzschildMetricChristoffelSymbolPartial( ) override = default;

    //! Retrieve function for partials w.r.t. integrated-body state.
    /*!
     *  \param stateReferencePoint Body/reference-point identifier.
     *  \param integratedStateType Integrated state type.
     *  \return Pair of (partial function, number of output columns).
     */
    std::pair< std::function< std::vector< std::vector< Eigen::Matrix< double, 4, 4 > > >( ) >, int >
        getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override;

    //! Retrieve function for partials w.r.t. scalar parameter.
    /*!
     *  \param parameter Scalar estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    //! Retrieve function for partials w.r.t. vector parameter.
    /*!
     *  \param parameter Vector estimatable parameter.
     *  \return Pair of (partial function, parameter size).
     */
    std::pair< std::function< std::vector< std::vector< Eigen::Matrix4d > >( ) >, int >
        getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    //! Update internal cached partial terms.
    /*!
     *  \param currentTime Current evaluation time.
     */
    void update( const double currentTime ) override;

private:
    std::vector< std::vector< Eigen::Matrix4d > > currentChristoffelSymbolsPositionPartials_;

    std::shared_ptr< relativity::HarmonicSchwarzschildMetric > metricObject_;

    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

    std::string centralBody_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_SCHWARZSCHILDMETRICCHRISTOFFELSYMBOLPARTIALS_H
