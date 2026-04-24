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

#ifndef TUDAT_CREATE_METRIC_H
#define TUDAT_CREATE_METRIC_H

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

//! Base settings class for space-time metric construction.
class SpaceTimeMetricSettings
{
public:
    //! Constructor.
    /*!
     *  \param metricType Identifier of metric model to construct.
     */
    explicit SpaceTimeMetricSettings( const SpaceTimeMetricTypes metricType )
        : metricType_( metricType ) { }

    virtual ~SpaceTimeMetricSettings( ) = default;

    //! Get configured metric model type.
    /*!
     *  \return Metric model type.
     */
    SpaceTimeMetricTypes getMetricType( ) const
    {
        return metricType_;
    }

protected:
    SpaceTimeMetricTypes metricType_;
};

//! Settings for the harmonic Schwarzschild metric.
class SchwardschildSpaceTimeMetricSettings : public SpaceTimeMetricSettings
{
public:
    //! Constructor.
    /*!
     *  \param bodyName Central body used in the Schwarzschild metric.
     *  \param ppnParameterSetToUse PPN parameter set used in metric evaluation.
     *  \param includeSecondPostNewtonianOrder If true, include second post-Newtonian terms.
     */
    SchwardschildSpaceTimeMetricSettings(
            const std::string& bodyName,
            const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSetToUse,
            const bool includeSecondPostNewtonianOrder = false )
        : SpaceTimeMetricSettings( schwardschild_metric ),
          bodyName_( bodyName ),
          ppnParameterSet_( ppnParameterSetToUse ),
          includeSecondPostNewtonianOrder_( includeSecondPostNewtonianOrder ) { }

    //! Get central body name.
    std::string getBodyName( ) const
    {
        return bodyName_;
    }

    //! Get PPN parameter set.
    std::shared_ptr< relativity::PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    //! Get second post-Newtonian inclusion flag.
    bool getIncludeSecondPostNewtonianOrder( ) const
    {
        return includeSecondPostNewtonianOrder_;
    }

private:
    std::string bodyName_;
    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;
    bool includeSecondPostNewtonianOrder_;
};

//! Settings for the solar-system post-Newtonian metric.
class SolarSystemSpaceTimeMetricSettings : public SpaceTimeMetricSettings
{
public:
    //! Constructor.
    /*!
     *  \param bodiesWithFirstOrderExpansion Bodies included at first post-Newtonian order.
     *  \param bodiesWithSecondOrderExpansion Bodies included in second-order terms.
     *  \param bodySphericalHarmonicExpansions Optional map of body name to maximum spherical harmonic degree/order.
     *  \param angularMomentumBodies Bodies for which angular momentum terms are included.
     *  \param ppnParameterSetToUse Optional PPN parameter set (default parameter set used when null).
     *  \param useBodyAccelerations If true, include explicit body acceleration terms where required.
     */
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

    //! Get list of first-order bodies.
    std::vector< std::string > getBodiesWithFirstOrderExpansion( ) const
    {
        return bodiesWithFirstOrderExpansion_;
    }

    //! Get list of second-order bodies.
    std::vector< std::string > getBodiesWithSecondOrderExpansion( ) const
    {
        return bodiesWithSecondOrderExpansion_;
    }

    //! Get spherical-harmonic expansion settings per body.
    std::map< std::string, std::pair< int, int > > getBodySphericalHarmonicExpansions( ) const
    {
        return bodySphericalHarmonicExpansions_;
    }

    //! Get list of bodies with angular-momentum terms.
    std::vector< std::string > getAngularMomentumBodies( ) const
    {
        return angularMomentumBodies_;
    }

    //! Get PPN parameter set.
    std::shared_ptr< relativity::PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    //! Get body-acceleration inclusion flag.
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

//! Factory function for creating a configured space-time metric object.
/*!
 *  \param spaceTimeMetricSettings Metric settings defining model type and configuration.
 *  \param bodyMap System of bodies providing required environment models/functions.
 *  \return Shared pointer to metric instance.
 */
std::shared_ptr< relativity::Metric > createSpaceTimeMetric(
        const std::shared_ptr< SpaceTimeMetricSettings >& spaceTimeMetricSettings,
        const simulation_setup::SystemOfBodies& bodyMap );

//! Factory function for Schwarzschild metric settings.
/*!
 *  \param bodyName Central body used in the Schwarzschild metric.
 *  \param includeSecondPostNewtonianOrder If true, include second post-Newtonian terms.
 *  \param ppnParameterSetToUse PPN parameter set used in metric evaluation.
 *  \return Shared pointer to Schwarzschild metric settings.
 */
inline std::shared_ptr< SchwardschildSpaceTimeMetricSettings > schwardschildSpaceTimeMetricSettings(
        const std::string& bodyName,
        const bool includeSecondPostNewtonianOrder = false,
        const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSetToUse = relativity::ppnParameterSet )
{
    return std::make_shared< SchwardschildSpaceTimeMetricSettings >(
                bodyName, ppnParameterSetToUse, includeSecondPostNewtonianOrder );
}

//! Factory function for solar-system metric settings.
/*!
 *  \param bodiesWithFirstOrderExpansion Bodies included at first post-Newtonian order.
 *  \param bodiesWithSecondOrderExpansion Bodies included in second-order terms.
 *  \param bodySphericalHarmonicExpansions Optional map of body name to maximum spherical harmonic degree/order.
 *  \param angularMomentumBodies Bodies for which angular momentum terms are included.
 *  \param useBodyAccelerations If true, include explicit body acceleration terms where required.
 *  \param ppnParameterSetToUse PPN parameter set used in metric evaluation.
 *  \return Shared pointer to solar-system metric settings.
 */
inline std::shared_ptr< SolarSystemSpaceTimeMetricSettings > solarSystemSpaceTimeMetricSettings(
        const std::vector< std::string >& bodiesWithFirstOrderExpansion,
        const std::vector< std::string >& bodiesWithSecondOrderExpansion = { },
        const std::map< std::string, std::pair< int, int > >& bodySphericalHarmonicExpansions = { },
        const std::vector< std::string >& angularMomentumBodies = { },
        const bool useBodyAccelerations = true,
        const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSetToUse = relativity::ppnParameterSet )
{
    return std::make_shared< SolarSystemSpaceTimeMetricSettings >(
                bodiesWithFirstOrderExpansion,
                bodiesWithSecondOrderExpansion,
                bodySphericalHarmonicExpansions,
                angularMomentumBodies,
                ppnParameterSetToUse,
                useBodyAccelerations );
}

//! Create and assign global base metric from metric settings.
/*!
 *  This function sets the global ``baseMetric`` used by direct-from-metric
 *  relativistic time propagation and clears cached cloned metric objects.
 *
 *  \param spaceTimeMetricSettings Metric settings defining model type and configuration.
 *  \param bodyMap System of bodies providing required environment models/functions.
 */
inline void createBaseMetric(
        const std::shared_ptr< SpaceTimeMetricSettings >& spaceTimeMetricSettings,
        const simulation_setup::SystemOfBodies& bodyMap )
{
    baseMetric = createSpaceTimeMetric( spaceTimeMetricSettings, bodyMap );
    evaluatedMetricObjects.clear( );
}

} // namespace simulation_setup
} // namespace tudat

#endif // TUDAT_CREATE_METRIC_H
