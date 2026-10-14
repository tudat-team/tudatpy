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

#include <functional>
#include <memory>
#include <stdexcept>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"

namespace tudat
{

namespace basic_astrodynamics
{

class BodyDeformationModel
{
public:
    virtual ~BodyDeformationModel( ) {}

    virtual Eigen::Vector3d calculateDisplacement( const double time, const Eigen::Vector3d& bodyFixedPosition ) = 0;

    virtual Eigen::Vector3d calculateDisplacement( const double time,
                                                   const std::shared_ptr< ground_stations::GroundStationState > stationState ) = 0;

    //! Environment context. Generic/custom models receive the original TDB scalar unchanged.
    virtual Eigen::Vector3d calculateDisplacement( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& bodyFixedPosition );
    virtual Eigen::Vector3d calculateDisplacement( const EarthDeformationEpoch& epoch,
                                                   const std::shared_ptr< ground_stations::GroundStationState > stationState );

    //! Body-fixed displacement velocity in metres per second at fixed nominal station coordinates.
    /*!
     * Velocity is the derivative with respect to the concrete model's displacement time argument,
     * in seconds, using that same epoch and time scale. No implicit time-scale conversion occurs.
     * For example, a displacement evaluated against TT produces a rate with respect to TT.
     * Existing displacement-only models inherit a five-point numerical derivative. The stencil
     * must remain in a smooth interval: jumps, interpolation knots and leap-second boundaries
     * require a model-specific override or a step that avoids the boundary. Callback coverage
     * errors propagate. Mutable model/callback caches require serial evaluation; the central
     * displacement is evaluated last to restore their epoch, also best-effort after an error.
     */
    virtual Eigen::Vector3d calculateDisplacementVelocity( const double time, const Eigen::Vector3d& bodyFixedPosition );

    virtual Eigen::Vector3d calculateDisplacementVelocity( const double time,
                                                           const std::shared_ptr< ground_stations::GroundStationState > stationState );

    //! Body-fixed displacement and velocity; models may override this combined evaluation.
    virtual Eigen::Vector6d calculateDisplacementState( const double time, const Eigen::Vector3d& bodyFixedPosition );

    virtual Eigen::Vector6d calculateDisplacementState( const double time,
                                                        const std::shared_ptr< ground_stations::GroundStationState > stationState );

    //! Context velocities/states are derivatives with respect to environment TDB.
    //! The generic default forwards the legacy scalar rate/state without conversion.
    virtual Eigen::Vector3d calculateDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& bodyFixedPosition );
    virtual Eigen::Vector3d calculateDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                           const std::shared_ptr< ground_stations::GroundStationState > stationState );
    virtual Eigen::Vector6d calculateDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& bodyFixedPosition );
    virtual Eigen::Vector6d calculateDisplacementState( const EarthDeformationEpoch& epoch,
                                                        const std::shared_ptr< ground_stations::GroundStationState > stationState );

    void setDisplacementDerivativeTimeStep( const double timeStep );

    double getDisplacementDerivativeTimeStep( ) const
    {
        return displacementDerivativeTimeStep_;
    }

protected:
    //! Assemble native state using context displacement/rate overrides, with the centre evaluated last.
    Eigen::Vector6d calculateNativeDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position );
    Eigen::Vector6d calculateNativeDisplacementState( const EarthDeformationEpoch& epoch,
                                                      const std::shared_ptr< ground_stations::GroundStationState > station );
    //! TT-only models: obtain their TT rate and apply the shared dTT/dTDB factor once.
    Eigen::Vector3d calculateTerrestrialDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position );
    Eigen::Vector3d calculateTerrestrialDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                              const std::shared_ptr< ground_stations::GroundStationState > station );
    //! Mixed-clock models: differentiate the complete context displacement with respect to TDB.
    Eigen::Vector3d calculateRoutedDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position );
    Eigen::Vector3d calculateRoutedDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                         const std::shared_ptr< ground_stations::GroundStationState > station );
    virtual void validateContextDerivativeInterval( const EarthDeformationEpoch& epoch,
                                                    const EarthDeformationEpoch& lower,
                                                    const EarthDeformationEpoch& upper ) const
    {
        validateDisplacementDerivativeInterval( epoch.getTdbSeconds( ), lower.getTdbSeconds( ), upper.getTdbSeconds( ) );
    }
    //! Model-specific interval policy; the generic default imposes no time-scale rules.
    virtual void validateDisplacementDerivativeInterval( const double, const double, const double ) const {}

private:
    double displacementDerivativeTimeStep_ = 1.0;
};

}  // namespace basic_astrodynamics

namespace ground_stations
{

struct BodyDeformationStationMotionModel : public StationMotionModel {
public:
    BodyDeformationStationMotionModel(
            std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >&( ) > modelList,
            const bool throwExceptionWhenNotAvailable,
            const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& timeConverter = nullptr ):
        modelList_( modelList ), throwExceptionWhenNotAvailable_( throwExceptionWhenNotAvailable ), timeConverter_( timeConverter )
    {}

    Eigen::Vector6d getBodyFixedStationMotion( const double time,
                                               const std::shared_ptr< ground_stations::GroundStationState > groundStationState,
                                               const std::string& targetFrameOrigin = "" )
    {
        Eigen::Vector6d motion = Eigen::Vector6d::Zero( );
        const basic_astrodynamics::EarthDeformationEpoch epoch( time, timeConverter_ );
        std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >& currentModels = modelList_( );
        for( unsigned int i = 0; i < currentModels.size( ); i++ )
        {
            if( !throwExceptionWhenNotAvailable_ )
            {
                // Preserve the optional central-displacement availability policy. Once the
                // central displacement exists, derivative/state failures must propagate.
                try
                {
                    currentModels.at( i )->calculateDisplacement( epoch, groundStationState );
                }
                catch( const std::exception& )
                {
                    continue;
                }
            }
            const Eigen::Vector6d modelState = currentModels.at( i )->calculateDisplacementState( epoch, groundStationState );
            if( !modelState.allFinite( ) )
            {
                throw std::runtime_error( "Non-finite body deformation state." );
            }
            motion += modelState;
        }
        return motion;
    }

    Eigen::Vector3d getBodyFixedStationPositionCorrection( const double time,
                                                           const std::shared_ptr< ground_stations::GroundStationState > groundStationState,
                                                           const std::string& targetFrameOrigin = "" ) override
    {
        Eigen::Vector3d correction = Eigen::Vector3d::Zero( );
        const basic_astrodynamics::EarthDeformationEpoch epoch( time, timeConverter_ );
        for( const auto& model : modelList_( ) )
        {
            Eigen::Vector3d displacement;
            try
            {
                displacement = model->calculateDisplacement( epoch, groundStationState );
            }
            catch( const std::exception& )
            {
                if( throwExceptionWhenNotAvailable_ )
                {
                    throw;
                }
                continue;
            }
            if( !displacement.allFinite( ) )
            {
                throw std::runtime_error( "Non-finite body deformation displacement." );
            }
            correction += displacement;
        }
        return correction;
    }

protected:
    std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >&( ) > modelList_;

    bool throwExceptionWhenNotAvailable_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeConverter_;
};

}  // namespace ground_stations

}  // namespace tudat

#endif  // TUDAT_BODYDEFORMATIONMODEL_H
