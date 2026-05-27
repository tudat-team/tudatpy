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

#ifndef TUDAT_SCHWARZSCHILDMETRIC_H
#define TUDAT_SCHWARZSCHILDMETRIC_H

#include <memory>
#include <string>
#include <functional>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/relativity/relativisticPotentials.h"
#include "tudat/astro/relativity/metric.h"

namespace tudat
{

namespace relativity
{

//! Harmonic-gauge Schwarzschild metric with optional second post-Newtonian terms.
class HarmonicSchwarzschildMetric : public Metric
{
public:
    //! Constructor.
    /*!
     *  \param centralGravitationalParameterFunction Function returning central-body gravitational parameter.
     *  \param ppnParameterSet PPN parameter set used in metric model.
     *  \param centralBodyName Name of central gravitating body.
     *  \param includeSecondPostNewtonianOrder If true, include second post-Newtonian terms.
     */
    HarmonicSchwarzschildMetric( const std::function< double( ) >& centralGravitationalParameterFunction,
                                 const std::shared_ptr< PPNParameterSet >& ppnParameterSet,
                                 const std::string& centralBodyName,
                                 const bool includeSecondPostNewtonianOrder = false ):
        centralGravitationalParameterFunction_( centralGravitationalParameterFunction ), ppnParameterSet_( ppnParameterSet ),
        centralBodyName_( centralBodyName ), includeSecondPostNewtonianOrder_( includeSecondPostNewtonianOrder )
    {}

    //! Copy constructor.
    /*!
     *  \param originalMetric Original metric object.
     */
    HarmonicSchwarzschildMetric( const HarmonicSchwarzschildMetric& originalMetric ):
        centralGravitationalParameterFunction_( originalMetric.centralGravitationalParameterFunction_ ),
        ppnParameterSet_( originalMetric.ppnParameterSet_ ), centralBodyName_( originalMetric.centralBodyName_ ),
        includeSecondPostNewtonianOrder_( originalMetric.includeSecondPostNewtonianOrder_ )
    {}

    //! Clone metric object.
    /*!
     *  \return Shared pointer to copied metric object.
     */
    std::shared_ptr< Metric > Clone( ) override
    {
        return std::make_shared< HarmonicSchwarzschildMetric >( *this );
    }

    //! Update metric and/or Christoffel symbols at requested state and time.
    /*!
     *  \param state Current translational state at evaluation point.
     *  \param time Current evaluation time.
     *  \param updateCurrentMetric If true, update covariant metric.
     *  \param updateCurrentChristoffelSymbols If true, update Christoffel symbols.
     */
    void update( const Eigen::Matrix< double, 6, 1 >& state,
                 const double time,
                 const bool updateCurrentMetric,
                 const bool updateCurrentChristoffelSymbols ) override;

    //! Get current central-body gravitational parameter.
    /*!
     *  \return Current gravitational parameter.
     */
    double getCurrentCentralGravitationalParameter( ) const
    {
        return currentCentralBodyGravitationalParameter_;
    }

    //! Get PPN parameter set.
    /*!
     *  \return PPN parameter set.
     */
    std::shared_ptr< PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    //! Check if second post-Newtonian terms are included.
    /*!
     *  \return True if second-order terms are included.
     */
    bool getIncludeSecondPostNewtonianOrder( ) const
    {
        return includeSecondPostNewtonianOrder_;
    }

    //! Get current first-order covariant metric contribution.
    /*!
     *  \return 4x4 first-order metric contribution.
     */
    Eigen::Matrix< double, 4, 4 > getCurrentFirstOrderCovariantMetricContributions( ) const
    {
        return currentFirstOrderCovariantMetricContributions_;
    }

    //! Get current second-order covariant metric contribution.
    /*!
     *  \return 4x4 second-order metric contribution.
     */
    Eigen::Matrix< double, 4, 4 > getCurrentSecondOrderCovariantMetricContributions( ) const
    {
        return currentSecondOrderCovariantMetricContributions_;
    }

    //! Get central body name.
    /*!
     *  \return Central body name.
     */
    std::string getCentralBodyName( ) const
    {
        return centralBodyName_;
    }

protected:
    void updateMetric( );

    void updateChristoffelSymbols( );

    std::function< double( ) > centralGravitationalParameterFunction_;
    std::shared_ptr< PPNParameterSet > ppnParameterSet_;
    std::string centralBodyName_;
    bool includeSecondPostNewtonianOrder_;

    Eigen::Matrix< double, 4, 4 > currentFirstOrderCovariantMetricContributions_{};
    Eigen::Matrix< double, 4, 4 > currentSecondOrderCovariantMetricContributions_{};
    double currentCentralBodyGravitationalParameter_{};
};

}  // namespace relativity

}  // namespace tudat

#endif  // TUDAT_SCHWARZSCHILDMETRIC_H
