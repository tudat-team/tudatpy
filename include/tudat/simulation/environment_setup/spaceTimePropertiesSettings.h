/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPACE_TIME_PROPERTIES_SETTINGS_H
#define TUDAT_SPACE_TIME_PROPERTIES_SETTINGS_H

#include <memory>
#include <stdexcept>

#include "tudat/astro/relativity/metric.h"

namespace tudat
{

namespace simulation_setup
{

class SpaceTimeMetricSettings;

class SpaceTimePropertiesSettings
{
public:
    SpaceTimePropertiesSettings( const std::shared_ptr< SpaceTimeMetricSettings >& metricSettings = nullptr,
                                 const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSet =
                                         std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 ),
                                 const double equivalencePrincipleLpiViolationParameter = 0.0 ):
        metricSettings_( metricSettings ), ppnParameterSet_( ppnParameterSet ),
        equivalencePrincipleLpiViolationParameter_( equivalencePrincipleLpiViolationParameter )
    {
        if( ppnParameterSet_ == nullptr )
        {
            ppnParameterSet_ = std::make_shared< relativity::PPNParameterSet >( 1.0, 1.0 );
        }
    }

    std::shared_ptr< SpaceTimeMetricSettings > getMetricSettings( ) const
    {
        return metricSettings_;
    }

    void setMetricSettings( const std::shared_ptr< SpaceTimeMetricSettings >& metricSettings )
    {
        metricSettings_ = metricSettings;
    }

    std::shared_ptr< relativity::PPNParameterSet > getPpnParameterSet( ) const
    {
        return ppnParameterSet_;
    }

    void setPpnParameterSet( const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSet )
    {
        if( ppnParameterSet == nullptr )
        {
            throw std::runtime_error( "Error when setting PPN parameter set, input is nullptr." );
        }
        ppnParameterSet_ = ppnParameterSet;
    }

    double getEquivalencePrincipleLpiViolationParameter( ) const
    {
        return equivalencePrincipleLpiViolationParameter_;
    }

    void setEquivalencePrincipleLpiViolationParameter( const double equivalencePrincipleLpiViolationParameter )
    {
        equivalencePrincipleLpiViolationParameter_ = equivalencePrincipleLpiViolationParameter;
    }

private:
    std::shared_ptr< SpaceTimeMetricSettings > metricSettings_;

    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

    double equivalencePrincipleLpiViolationParameter_;
};

//! Factory function for PPN parameter settings.
/*!
 *  \param parameterGamma First-order PPN parameter gamma.
 *  \param parameterBeta First-order PPN parameter beta.
 *  \param parameterEpsilon Second-order post-Newtonian parameter epsilon.
 *  \return Shared pointer to PPN parameter settings.
 */
inline std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet( const double parameterGamma = 1.0,
                                                                       const double parameterBeta = 1.0,
                                                                       const double parameterEpsilon = 0.0 )
{
    return std::make_shared< relativity::PPNParameterSet >( parameterGamma, parameterBeta, 0.0, parameterEpsilon );
}

//! Factory function for system-level space-time properties settings.
/*!
 *  \param metricSettings Optional settings defining which metric model to build.
 *  \param ppnParameterSet Optional PPN parameter settings. If null, GR defaults are used.
 *  \param equivalencePrincipleLpiViolationParameter Equivalence-principle local-position-invariance violation parameter.
 *  \return Shared pointer to space-time properties settings.
 */
inline std::shared_ptr< SpaceTimePropertiesSettings > spaceTimePropertiesSettings(
        const std::shared_ptr< SpaceTimeMetricSettings >& metricSettings = nullptr,
        const std::shared_ptr< relativity::PPNParameterSet >& ppnParameterSet = nullptr,
        const double equivalencePrincipleLpiViolationParameter = 0.0 )
{
    return std::make_shared< SpaceTimePropertiesSettings >( metricSettings, ppnParameterSet, equivalencePrincipleLpiViolationParameter );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_SPACE_TIME_PROPERTIES_SETTINGS_H
