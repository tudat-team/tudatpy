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

    //! Displacement of the station (body-fixed Cartesian, in metres) at the given TDB epoch (seconds since J2000).
    Eigen::Vector3d calculateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station );

    //! Displacement in metres of a site given only by its body-fixed Cartesian position (existing interface).
    //! Models that need a station identifier, such as ocean tide loading, cannot evaluate this.
    Eigen::Vector3d calculateDisplacement( const Time& tdb, const Eigen::Vector3d& bodyFixedPosition );

    //! Displacement of the station (in metres) and its time derivative (in metres per TDB second), at the given TDB epoch.
    //! A model that cannot compute the velocity at an epoch (for instance because its data has a gap there) throws an
    //! exception here, but calculateDisplacement still works at that epoch.
    Eigen::Vector6d calculateDisplacementState( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station );

    //! Set the time step (in seconds) of the numerical derivative used for models without an analytical velocity.
    void setDisplacementDerivativeTimeStep( double timeStep );

    //! Time step (in seconds) of the numerical derivative used for models without an analytical velocity.
    double getDisplacementDerivativeTimeStep( ) const
    {
        return displacementDerivativeTimeStep_;
    }

protected:
    //! Constructor; the time scale converter is used to convert TDB to the time scale used by the model. If it is not
    //! provided, getGeocentricTimeConverter( ) is used.
    explicit BodyDeformationModel( const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter = nullptr );

    //! Displacement computed by the derived class (called by calculateDisplacement, which checks its input and output).
    virtual Eigen::Vector3d evaluateDisplacement( const Time& tdb,
                                                  const std::shared_ptr< ground_stations::GroundStationState > station ) = 0;

    //! Displacement and velocity computed by the derived class. By default, the velocity is computed numerically from the
    //! displacement at four nearby epochs (fourth-order central difference); derived classes with an analytical velocity
    //! override this function.
    virtual Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                                       const std::shared_ptr< ground_stations::GroundStationState > station );

    //! Numerical derivative (fourth-order central difference) of the given displacement function at the given epoch. If
    //! validateInterval is provided, it is first called with the first and last epoch at which the displacement is evaluated.
    Eigen::Vector3d differentiateDisplacement( const Time& time,
                                               const std::function< Eigen::Vector3d( const Time& ) >& displacement,
                                               const std::function< void( const Time&, const Time& ) >& validateInterval ) const;
    Eigen::Vector3d differentiateDisplacement( double time,
                                               const std::function< Eigen::Vector3d( const double& ) >& displacement,
                                               const std::function< void( const double&, const double& ) >& validateInterval ) const;

    //! Check that the numerical derivative of the displacement can be computed over the given TDB interval; throws an
    //! exception if it cannot. By default, any interval is accepted.
    virtual void validateDisplacementDerivativeInterval( const Time&, const Time& ) const {}

    //! Object to convert between time scales.
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeConverter_;

private:
    friend struct ground_stations::BodyDeformationStationMotionModel;

    //! Time step (seconds) of the numerical derivative of the displacement.
    double displacementDerivativeTimeStep_ = 1.0;
};

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
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& tdb,
                                                           const std::shared_ptr< GroundStationState > station,
                                                           const std::string& origin = "" ) override
    {
        return evaluateMotion< false >( tdb, station );
    }

private:
    //! Sum of the displacements of all deformation models (and of their velocities, if State is true). If
    //! throwExceptionWhenNotAvailable_ is false, models that have no data at this epoch are skipped.
    template< bool State >
    Eigen::Matrix< double, State ? 6 : 3, 1 > evaluateMotion( const Time& tdb, const std::shared_ptr< GroundStationState > station )
    {
        Eigen::Matrix< double, State ? 6 : 3, 1 > result = Eigen::Matrix < double, State ? 6 : 3, 1 > ::Zero( );
        if( !std::isfinite( static_cast< double >( tdb ) ) )
        {
            throw std::invalid_argument( "Body deformation requires a finite TDB epoch." );
        }
        if( !station )
        {
            throw std::invalid_argument( "Body deformation requires a station." );
        }
        const auto& models = modelList_( );
        for( std::size_t i = 0; i < models.size( ); ++i )
        {
            const auto& model = models[ i ];
            Eigen::Vector3d displacement;
            try
            {
                if( !State || !throwExceptionWhenNotAvailable_ )
                {
                    displacement = model->evaluateDisplacement( tdb, station );
                }
            }
            catch( const std::out_of_range& )
            {
                if( throwExceptionWhenNotAvailable_ )
                {
                    throw;
                }
                continue;
            }
            if( ( !State || !throwExceptionWhenNotAvailable_ ) && !displacement.allFinite( ) )
            {
                throw std::runtime_error( "Non-finite body deformation displacement." );
            }

            // A model that has data at this epoch can still fail to compute its velocity; that exception is not caught.
            if constexpr( State )
            {
                result += model->calculateDisplacementState( tdb, station );
            }
            else
            {
                result += displacement;
            }
        }
        return result;
    }

    std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >&( ) > modelList_;
    bool throwExceptionWhenNotAvailable_;
};

}  // namespace ground_stations

}  // namespace tudat

#endif  // TUDAT_BODYDEFORMATIONMODEL_H
