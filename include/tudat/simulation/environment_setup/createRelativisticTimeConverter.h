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

#ifndef TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H
#define TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H

#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <iostream>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/propagators/relativisticTimeStateDerivative.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"

namespace tudat
{
namespace simulation_setup
{

template< typename StateScalarType = double, typename TimeType = double >
class DirectRelativisticTimeConverterSettings
{
public:
    //! Constructor.
    /*!
     *  \param baryCentricToBodyCentricConversionSettings Settings defining the barycentric-to-bodycentric
     *  conversion model for the associated body.
     *  \param numericalIntegrationSettings Integration settings used by the internal relativistic-time propagations.
     *  \param bodyCentricToTopocentricConversionSettings Optional set of bodycentric-to-topocentric settings
     *  (typically one per ground station/reference point on the associated body).
     */
    DirectRelativisticTimeConverterSettings(
            const std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > >&
                    baryCentricToBodyCentricConversionSettings,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& numericalIntegrationSettings,
            const std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > >&
                    bodyCentricToTopocentricConversionSettings = std::vector<
                            std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > >( ) ):
        baryCentricToBodyCentricConversionSettings_( baryCentricToBodyCentricConversionSettings ),
        bodyCentricToTopocentricConversionSettings_( bodyCentricToTopocentricConversionSettings ),
        numericalIntegrationSettings_( numericalIntegrationSettings )
    {
        using namespace propagators;

        if( !( baryCentricToBodyCentricConversionSettings_->getRelativisticStateDerivativeType( ) ==
                       first_order_barycentric_to_bodycentric ||
               baryCentricToBodyCentricConversionSettings_->getRelativisticStateDerivativeType( ) ==
                       second_order_barycentric_to_bodycentric ) )
        {
            throw std::runtime_error( "[ERROR] Invalid type for barycentric-to-bodycentric conversion: " +
                                      std::to_string( static_cast< int >(
                                              baryCentricToBodyCentricConversionSettings_->getRelativisticStateDerivativeType( ) ) ) );
        }

        for( unsigned int i = 0; i < bodyCentricToTopocentricConversionSettings_.size( ); ++i )
        {
            auto topocentricSetting = bodyCentricToTopocentricConversionSettings_.at( i );
            if( topocentricSetting->getRelativisticStateDerivativeType( ) != first_order_bodycentric_to_topocentric )
            {
                throw std::runtime_error(
                        "[ERROR] Invalid type for topocentric setting for station " + topocentricSetting->getReferencePointId( ).second +
                        ": " + std::to_string( static_cast< int >( topocentricSetting->getRelativisticStateDerivativeType( ) ) ) );
            }

            if( topocentricSetting->getReferencePointId( ).first !=
                baryCentricToBodyCentricConversionSettings_->getReferencePointId( ).first )
            {
                throw std::runtime_error( "[ERROR] Station " + topocentricSetting->getReferencePointId( ).second + " is on body " +
                                          topocentricSetting->getReferencePointId( ).first + ", but barycentric setting uses body " +
                                          baryCentricToBodyCentricConversionSettings_->getReferencePointId( ).first );
            }
        }

        associatedBody_ = baryCentricToBodyCentricConversionSettings_->getReferencePointId( ).first;
    }

    //! Get barycentric-to-bodycentric conversion settings.
    /*!
     *  \return Settings object for the barycentric/bodycentric conversion leg.
     */
    std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > >
    getBaryCentricToBodyCentricConversionSettings( ) const
    {
        return baryCentricToBodyCentricConversionSettings_;
    }

    //! Get bodycentric-to-topocentric conversion settings.
    /*!
     *  \return Vector of settings for topocentric conversion legs.
     */
    std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > >
    getBodyCentricToTopocentricConversionSettings( ) const
    {
        return bodyCentricToTopocentricConversionSettings_;
    }

    //! Get numerical integration settings.
    /*!
     *  \return Integrator settings used in converter construction.
     */
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getNumericalIntegrationSettings( ) const
    {
        return numericalIntegrationSettings_;
    }

    //! Get name of body associated with this converter.
    /*!
     *  \return Associated body name.
     */
    std::string getAssociatedBody( ) const
    {
        return associatedBody_;
    }

private:
    std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > >
            baryCentricToBodyCentricConversionSettings_;

    std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > >
            bodyCentricToTopocentricConversionSettings_;

    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > numericalIntegrationSettings_;

    std::string associatedBody_;
};

//! Setup a single relativistic time converter (barycentric -> bodycentric -> topocentric).
/*!
 *  \param conversionSettings Settings defining the converter and all conversion legs.
 *  \param bodies System of bodies in which the converter is created and stored.
 */
template< typename StateScalarType = double, typename TimeType = double >
void setRelativisticTimeConverter(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > >& conversionSettings,
        const SystemOfBodies& bodies );

//! Setup multiple relativistic time converters (one per associated body).
/*!
 *  \param bodies System of bodies in which converters are created and stored.
 *  \param converterSettings Map from body name to converter settings.
 */
template< typename StateScalarType = double, typename TimeType = double >
void setRelativisticTimeConverters(
        const SystemOfBodies& bodies,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > > >&
                converterSettings );

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H
