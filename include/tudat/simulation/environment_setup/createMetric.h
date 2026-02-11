#ifndef CREATEMETRIC_H
#define CREATEMETRIC_H

#include <memory>
#include <map>
#include <string>
#include <vector>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/relativity/metric.h"

namespace tudat
{

extern std::map< std::pair< std::string, std::string >, std::shared_ptr< relativity::Metric > > evaluatedMetricObjects;

namespace simulation_setup
{

enum SpaceTimeMetricTypes
{
    schwardschild_metric,
    solar_system_metric
};

// Base settings class
class SpaceTimeMetricSettings
{
public:
    explicit SpaceTimeMetricSettings( const SpaceTimeMetricTypes metricType )
        : metricType_( metricType ) { }

    virtual ~SpaceTimeMetricSettings( ) = default;

    SpaceTimeMetricTypes getMetricType( ) const
    {
        return metricType_;
    }

protected:
    SpaceTimeMetricTypes metricType_;
};

// Schwardschild metric settings
class SchwardschildSpaceTimeMetricSettings : public SpaceTimeMetricSettings
{
public:
    SchwardschildSpaceTimeMetricSettings(
            const std::string& bodyName,
            const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSetToUse,
            const bool includeSecondPostNewtonianOrder = false )
        : SpaceTimeMetricSettings( schwardschild_metric ),
          bodyName_( bodyName ),
          ppnParameterSet_( ppnParameterSetToUse ),
          includeSecondPostNewtonianOrder_( includeSecondPostNewtonianOrder ) { }

    std::string getBodyName( ) const
    {
        return bodyName_;
    }

    std::shared_ptr< relativity::PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    bool getIncludeSecondPostNewtonianOrder( ) const
    {
        return includeSecondPostNewtonianOrder_;
    }

private:
    std::string bodyName_;
    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;
    bool includeSecondPostNewtonianOrder_;
};

// Solar System metric settings
class SolarSystemSpaceTimeMetricSettings : public SpaceTimeMetricSettings
{
public:
    SolarSystemSpaceTimeMetricSettings(
            const std::vector< std::string >& bodiesWithFirstOrderExpansion,
            const std::vector< std::string >& bodiesWithSecondOrderExpansion = { },
            const std::map< std::string, std::pair< int, int > >& bodySphericalHarmonicExpansions = { },
            const std::vector< std::string >& angularMomentumBodies = { },
            const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSetToUse = nullptr,
            const bool useBodyAccelerations = true )
        : SpaceTimeMetricSettings( solar_system_metric ),
          bodiesWithFirstOrderExpansion_( bodiesWithFirstOrderExpansion ),
          bodiesWithSecondOrderExpansion_( bodiesWithSecondOrderExpansion ),
          bodySphericalHarmonicExpansions_( bodySphericalHarmonicExpansions ),
          angularMomentumBodies_( angularMomentumBodies ),
          ppnParameterSet_( ppnParameterSetToUse ),
          useBodyAccelerations_( useBodyAccelerations ) { }

    std::vector< std::string > getBodiesWithFirstOrderExpansion( ) const
    {
        return bodiesWithFirstOrderExpansion_;
    }

    std::vector< std::string > getBodiesWithSecondOrderExpansion( ) const
    {
        return bodiesWithSecondOrderExpansion_;
    }

    std::map< std::string, std::pair< int, int > > getBodySphericalHarmonicExpansions( ) const
    {
        return bodySphericalHarmonicExpansions_;
    }

    std::vector< std::string > getAngularMomentumBodies( ) const
    {
        return angularMomentumBodies_;
    }

    std::shared_ptr< relativity::PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    bool getUseBodyAccelerations( ) const
    {
        return useBodyAccelerations_;
    }

private:
    std::vector< std::string > bodiesWithFirstOrderExpansion_;
    std::vector< std::string > bodiesWithSecondOrderExpansion_;
    std::map< std::string, std::pair< int, int > > bodySphericalHarmonicExpansions_;
    std::vector< std::string > angularMomentumBodies_;
    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;
    bool useBodyAccelerations_;
};

std::shared_ptr< relativity::Metric > createSpaceTimeMetric(
        const std::shared_ptr< SpaceTimeMetricSettings >& spaceTimeMetricSettings,
        const simulation_setup::SystemOfBodies& bodyMap );

} // namespace simulation_setup
} // namespace tudat

#endif // CREATEMETRIC_H
