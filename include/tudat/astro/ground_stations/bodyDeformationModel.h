/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#ifndef TUDAT_BODYDEFORMATIONMODEL_H
#define TUDAT_BODYDEFORMATIONMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <type_traits>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"

namespace tudat
{

namespace ground_stations
{
struct BodyDeformationStationMotionModel;
}

namespace basic_astrodynamics
{

class BodyDeformationModel
{
public:
    virtual ~BodyDeformationModel( ) {}

    //! Body-fixed Cartesian displacement in metres at the supplied precise TDB epoch.
    Eigen::Vector3d calculateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station );
    //! Displacement and its d/dTDB velocity, in metres and metres per second.
    //! Position evaluation remains independent of rate availability.
    Eigen::Vector6d calculateDisplacementState( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station );

    void setDisplacementDerivativeTimeStep( double timeStep );
    double getDisplacementDerivativeTimeStep( ) const
    {
        return displacementDerivativeTimeStep_;
    }

protected:
    explicit BodyDeformationModel( const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter = nullptr );
    virtual Eigen::Vector3d evaluateDisplacement( const Time& tdb,
                                                  const std::shared_ptr< ground_stations::GroundStationState > station ) = 0;
    //! Guarded five-point precise-Time fallback; analytic models override this hook.
    virtual Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                                       const std::shared_ptr< ground_stations::GroundStationState > station );
    Eigen::Vector3d differentiateDisplacement( const Time& time,
                                               const std::function< Eigen::Vector3d( const Time& ) >& displacement,
                                               const std::function< void( const Time&, const Time& ) >& validateInterval ) const;
    Eigen::Vector3d differentiateDisplacement( double time,
                                               const std::function< Eigen::Vector3d( const double& ) >& displacement,
                                               const std::function< void( const double&, const double& ) >& validateInterval ) const;
    virtual void validateDisplacementDerivativeInterval( const Time&, const Time& ) const {}
    //! Only the original-double station boundary invokes this adapter.
    virtual Time normalizeScalarTdb( double tdb ) const;

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeConverter_;

private:
    friend struct ground_stations::BodyDeformationStationMotionModel;
    double displacementDerivativeTimeStep_ = 1.0;
};

void validateStationLoadingComposition( const std::vector< std::shared_ptr< BodyDeformationModel > >& models,
                                        const std::vector< bool >& selected,
                                        const std::shared_ptr< ground_stations::GroundStationState >& station );

}  // namespace basic_astrodynamics

namespace ground_stations
{

struct BodyDeformationStationMotionModel : public StationMotionModel {
public:
    BodyDeformationStationMotionModel(
            std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >&( ) > modelList,
            const bool throwExceptionWhenNotAvailable ):
        modelList_( modelList ), throwExceptionWhenNotAvailable_( throwExceptionWhenNotAvailable )
    {}

    Eigen::Vector6d getBodyFixedStationMotion( const Time& tdb,
                                               const std::shared_ptr< GroundStationState > station,
                                               const std::string& origin = "" ) override
    {
        return evaluateMotion< true >( tdb, station );
    }
    Eigen::Vector6d getBodyFixedStationMotion( double tdb,
                                               const std::shared_ptr< GroundStationState > station,
                                               const std::string& origin = "" ) override
    {
        return evaluateMotion< true >( tdb, station );
    }
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& tdb,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin = "" ) override
    {
        return evaluateMotion< false >( tdb, station );
    }
    Eigen::Vector3d getBodyFixedStationPositionCorrection( double tdb,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin = "" ) override
    {
        return evaluateMotion< false >( tdb, station );
    }

private:
    template< bool State, typename Epoch >
    Eigen::Matrix< double, State ? 6 : 3, 1 > evaluateMotion( const Epoch& input, const std::shared_ptr< GroundStationState > station )
    {
        Eigen::Matrix< double, State ? 6 : 3, 1 > result = Eigen::Matrix < double, State ? 6 : 3, 1 > ::Zero( );
        if( !std::isfinite( static_cast< double >( input ) ) )
            throw std::invalid_argument( "Body deformation requires a finite TDB epoch." );
        if( !station ) throw std::invalid_argument( "Body deformation requires a station." );
        const auto& models = modelList_( );
        std::vector< bool > selected( models.size( ), false );
        for( std::size_t i = 0; i < models.size( ); ++i )
        {
            const auto& model = models[ i ];
            Time tdb;
            Eigen::Vector3d displacement;
            try
            {
                if constexpr( std::is_same< Epoch, double >::value )
                    tdb = model->normalizeScalarTdb( input );
                else
                    tdb = input;
                if( !State || !throwExceptionWhenNotAvailable_ ) displacement = model->evaluateDisplacement( tdb, station );
            }
            catch( const std::invalid_argument& )
            {
                throw;
            }
            catch( const std::exception& )
            {
                if( throwExceptionWhenNotAvailable_ ) throw;
                continue;
            }
            if( ( !State || !throwExceptionWhenNotAvailable_ ) && !displacement.allFinite( ) )
                throw std::runtime_error( "Non-finite body deformation displacement." );
            if constexpr( State )
                result += model->calculateDisplacementState( tdb, station );
            else
                result += displacement;
            selected[ i ] = true;
        }
        basic_astrodynamics::validateStationLoadingComposition( models, selected, station );
        return result;
    }

    std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >&( ) > modelList_;
    bool throwExceptionWhenNotAvailable_;
};

}  // namespace ground_stations

}  // namespace tudat

#endif  // TUDAT_BODYDEFORMATIONMODEL_H
