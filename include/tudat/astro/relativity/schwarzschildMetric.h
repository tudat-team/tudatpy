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

class HarmonicSchwarzschildMetric : public Metric
{
public:
    HarmonicSchwarzschildMetric(
        const std::function< double( ) >& centralGravitationalParameterFunction,
        const std::shared_ptr< PPNParameterSet >& ppnParameterSet,
        const std::string& centralBodyName,
        const bool includeSecondPostNewtonianOrder = false )
        : centralGravitationalParameterFunction_( centralGravitationalParameterFunction ),
          ppnParameterSet_( ppnParameterSet ),
          centralBodyName_( centralBodyName ),
          includeSecondPostNewtonianOrder_( includeSecondPostNewtonianOrder ) { }

    HarmonicSchwarzschildMetric( const HarmonicSchwarzschildMetric& originalMetric )
        : centralGravitationalParameterFunction_( originalMetric.centralGravitationalParameterFunction_ ),
          ppnParameterSet_( originalMetric.ppnParameterSet_ ),
          centralBodyName_( originalMetric.centralBodyName_ ),
          includeSecondPostNewtonianOrder_( originalMetric.includeSecondPostNewtonianOrder_ ) { }

    std::shared_ptr< Metric > Clone( ) override
    {
        return std::make_shared< HarmonicSchwarzschildMetric >( *this );
    }

    void update( const Eigen::Matrix< double, 6, 1 >& state,
                 const double time,
                 const bool updateCurrentMetric,
                 const bool updateCurrentChristoffelSymbols ) override;

    double getCurrentCentralGravitationalParameter( ) const
    {
        return currentCentralBodyGravitationalParameter_;
    }

    std::shared_ptr< PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    bool getIncludeSecondPostNewtonianOrder( ) const
    {
        return includeSecondPostNewtonianOrder_;
    }

    Eigen::Matrix< double, 4, 4 > getCurrentFirstOrderCovariantMetricContributions( ) const
    {
        return currentFirstOrderCovariantMetricContributions_;
    }

    Eigen::Matrix< double, 4, 4 > getCurrentSecondOrderCovariantMetricContributions( ) const
    {
        return currentSecondOrderCovariantMetricContributions_;
    }

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

    Eigen::Matrix< double, 4, 4 > currentFirstOrderCovariantMetricContributions_{ };
    Eigen::Matrix< double, 4, 4 > currentSecondOrderCovariantMetricContributions_{ };
    double currentCentralBodyGravitationalParameter_{ };
};

} // namespace relativity

} // namespace tudat

#endif // TUDAT_SCHWARZSCHILDMETRIC_H
