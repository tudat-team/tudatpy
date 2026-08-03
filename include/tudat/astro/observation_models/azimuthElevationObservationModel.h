/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AZIMUTHELEVATIONOBSERVATIONMODEL_H
#define TUDAT_AZIMUTHELEVATIONOBSERVATIONMODEL_H

#include <map>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating azimuth/elevation observables at a ground station.
template< typename ObservationScalarType = double, typename TimeType = double >
class AzimuthElevationObservationModel : public ObservationModel< 2, ObservationScalarType, TimeType >
{
public:
    AzimuthElevationObservationModel(
            const LinkEnds linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator,
            const std::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = nullptr,
            const bool normalizeAzimuth = false ):
        ObservationModel< 2, ObservationScalarType, TimeType >(
                azimuth_elevation_angle,
                linkEnds,
                observationBiasCalculator,
                std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{ std::make_shared<
                        FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                        std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                                lightTimeCalculator },
                        std::make_shared< LightTimeConvergenceCriteria >( ),
                        false ) } ),
        pointingAnglesCalculator_( pointingAnglesCalculator ), normalizeAzimuth_( normalizeAzimuth )
    {
        if( linkEnds.count( receiver ) == 0 || linkEnds.at( receiver ).getReferencePointName( ) == "" )
        {
            throw std::runtime_error( "Error when creating azimuth/elevation model: receiver link end must be a ground station." );
        }
        if( pointingAnglesCalculator_ == nullptr )
        {
            throw std::runtime_error( "Error when creating azimuth/elevation model: pointing angle calculator is null." );
        }
    }

    ~AzimuthElevationObservationModel( ) {}

    Eigen::Matrix< ObservationScalarType, 2, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetingsInput = nullptr ) override
    {
        std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetings;
        std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator =
                getLightTimeCalculator( );
        this->setFrequencyProperties( time, linkEndAssociatedWithTime, lightTimeCalculator, ancillarySetingsInput, ancillarySetings );

        if( ancillarySetings != nullptr )
        {
            throw std::runtime_error( "Error, calling azimuth/elevation observable with ancillary settings, but none are supported." );
        }

        this->getFullLinkLightTimeCalculatorFromBase( )->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySetings );

        Eigen::Vector3d relativePosition = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment( 0, 3 );

        std::pair< double, double > elevationAzimuth =
                pointingAnglesCalculator_->calculatePointingAngles( relativePosition, linkEndTimes.at( 1 ) );

        double azimuth = elevationAzimuth.second;
        if( normalizeAzimuth_ )
        {
            azimuth = basic_mathematics::computeModulo( azimuth, 2.0 * mathematical_constants::PI );
        }

        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << azimuth, elevationAzimuth.first ).finished( );
    }

    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return this->getSingleLegLightTimeCalculator( );
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( ) } } };
    }

    std::shared_ptr< ground_stations::PointingAnglesCalculator > getPointingAnglesCalculator( )
    {
        return pointingAnglesCalculator_;
    }

    bool getNormalizeAzimuth( )
    {
        return normalizeAzimuth_;
    }

private:
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator_;

    bool normalizeAzimuth_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_AZIMUTHELEVATIONOBSERVATIONMODEL_H
