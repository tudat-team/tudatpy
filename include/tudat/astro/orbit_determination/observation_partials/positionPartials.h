/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONPARTIALS_H
#define TUDAT_POSITIONPARTIALS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate the partial of position of a point on body wrt position of that body.
/*!
 *  Function to calculate the partial of position of a point on body wrt position of that body,
 *  with both positions expressed in the same frame.
 *  \return Requested partial (3x3 identity matrix)
 */
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyPosition( );

//! Function to calculate the partial of position of a point on a body wrt its body-fixed position
/*!
 *  Function to calculate the partial of position of a point on a body (expressed in a non-corotating, non-body fixed frame)
 *  wrt its body-fixed position, e.g. the partial of a ground station's inertial
 *  position wrt its body-fixed position.
 *  \param rotationMatrixToInertialFrame Rotation matrix from body-fixed to inertial frame.
 *  \return Partial of position of a point on a body wrt its body-fixed position
 */
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyFixedPointPosition( const Eigen::Matrix3d& rotationMatrixToInertialFrame );

//! Base class for calculating the partial of an inertial position wrt a parameter.
/*!
 *  Base class for calculating the partial of an inertial position wrt a parameter. A derived class is implemented for
 *  each (type of) estimatable parameter. A separate instance of the class must be made for each distinct state.
 *  Note that partials wrt parameters describing a property of a rotation matrix from a local to the inertial frame is
 *  implemented in the RotationMatrixPartial class, which is
 *  then used by the CartesianStatePartialWrtRotationMatrixParameter derived class of this class
 */
class CartesianStatePartial
{
public:
    CartesianStatePartial( ) {}

    //! Destructor.
    virtual ~CartesianStatePartial( ) {}

    //! Pure virtual base class function for determining partial at current time and body state.
    /*!
     *  Pure virtual base class function for determining partial at current time and body state wrt parameter of specific
     *  implemented derived. class.
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point position wrt parameter (with specific parameter determined by derived class implementation).
     */
    virtual Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition( const Eigen::Vector6d& state, const double time ) = 0;

    virtual Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity( const Eigen::Vector6d& state, const double time ) = 0;
};

//! Class to compute the partial derivative of the Cartesian state of a body w.r.t. to inertial three-dimensional
//! position of this body
class CartesianStatePartialWrtCartesianState : public CartesianStatePartial
{
public:
    //! Constructor
    CartesianStatePartialWrtCartesianState( )
    {
        positionPartial_.setZero( );
        positionPartial_.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Identity( );
        positionPartial_.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Zero( );

        velocityPartial_.setZero( );
        velocityPartial_.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Zero( );
        velocityPartial_.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Identity( );
    }

    //! Destructor
    ~CartesianStatePartialWrtCartesianState( ) {}

    //! Function for determining partial of position at current time and body state.
    /*!
     *  Function for determining partial of position  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt position
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition( const Eigen::Vector6d& state, const double time )
    {
        return positionPartial_;
    }
    //! Function for determining partial of velocity at current time and body state.
    /*!
     *  Function for determining partial of velocity  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt velocity
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity( const Eigen::Vector6d& state, const double time )
    {
        return velocityPartial_;
    }

private:
    //! Partial of current state w.r.t. position
    Eigen::Matrix< double, 3, 6 > positionPartial_;

    //! Partial of current state w.r.t. velocity
    Eigen::Matrix< double, 3, 6 > velocityPartial_;
};

//! Class to compute the partial derivative of the Cartesian state of a body w.r.t. to a property of a rotation
//! matrix to/from a body-fixed frame.
class CartesianStatePartialWrtRotationMatrixParameter : public CartesianStatePartial
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param rotationMatrixPartialObject Object to compute the associated partial of a rotation matrix
     * \param positionFunctionInLocalFrame Function returning the body-fixed position of the point at which the partial
     * is to be computed.
     */
    CartesianStatePartialWrtRotationMatrixParameter( const std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject,
                                                     const std::function< Eigen::Vector3d( const double ) > positionFunctionInLocalFrame ):
        rotationMatrixPartialObject_( rotationMatrixPartialObject ), positionFunctionInLocalFrame_( positionFunctionInLocalFrame )
    {}

    //! Destructor
    ~CartesianStatePartialWrtRotationMatrixParameter( ) {}

    //! Function for determining partial of position at current time and body state.
    /*!
     *  Function for determining partial of position  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt rotation property
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition( const Eigen::Vector6d& state, const double time )
    {
        return rotationMatrixPartialObject_->calculatePartialOfInertialPositionWrtParameter( time, positionFunctionInLocalFrame_( time ) );
    }

    //! Function for determining partial of velocity at current time and body state.
    /*!
     *  Function for determining partial of velocity  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt rotation property
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity( const Eigen::Vector6d& state, const double time )
    {
        return rotationMatrixPartialObject_->calculatePartialOfInertialVelocityWrtParameter( time, positionFunctionInLocalFrame_( time ) );
    }

private:
    //! Object to compute the associated partial of a rotation matrix
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject_;

    //! Function returning the body-fixed position of the point at which the partial is to be computed.
    std::function< Eigen::Vector3d( const double ) > positionFunctionInLocalFrame_;
};

//! Class to compute the partial derivative of the inertial Cartesian state of a point on a body w.r.t. the constant
//! body-fixed position of that point.
class CartesianPartialWrtBodyFixedPosition : public CartesianStatePartial
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param bodyRotationModel Rotation model for body.
     */
    CartesianPartialWrtBodyFixedPosition( const std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel ):
        bodyRotationModel_( bodyRotationModel )
    {}

    //! Destructor
    ~CartesianPartialWrtBodyFixedPosition( ) {}

    //! Function for determining partial of position at current time and body state.
    /*!
     *  Function for determining partial of position  at current time and body state wrt three-dimensional body-fixed state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt body-fixed position.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition( const Eigen::Vector6d& state, const double time )
    {
        return Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) );
    }

    //! Function for determining partial of velocity at current time and body state.
    /*!
     *  Function for determining partial of velocity  at current time and body state wrt three-dimensional state
     *  \param state Current inertial state of point of which partial is to be calculated.
     *  \param time Current time
     *  \return Partial of point state wrt body-fixed position.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity( const Eigen::Vector6d& state, const double time )
    {
        return bodyRotationModel_->getDerivativeOfRotationToBaseFrame( time );
    }

private:
    //! Rotation model for body.
    std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel_;
};

//! Partial of the state of a ground station w.r.t. the body-fixed velocity of its linear motion model.
class CartesianPartialWrtStationLinearVelocity : public CartesianStatePartial
{
public:
    CartesianPartialWrtStationLinearVelocity( const std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel,
                                              const double referenceEpoch ):
        bodyRotationModel_( bodyRotationModel ), referenceEpoch_( referenceEpoch )
    {}

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition( const Eigen::Vector6d& state, const double time )
    {
        return Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) ) * ( time - referenceEpoch_ );
    }

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity( const Eigen::Vector6d& state, const double time )
    {
        return bodyRotationModel_->getDerivativeOfRotationToBaseFrame( time ) * ( time - referenceEpoch_ ) +
                Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) );
    }

private:
    std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel_;

    double referenceEpoch_;
};

//! Partial of the state of a point on the surface of a body w.r.t. the displacement Love number h_n or Shida number l_n of the basic
//! tidal deformation models of the body. The displacement is linear in this number, so that its partial is the displacement computed
//! with a unit Love or Shida number. The velocity partial is that of the rotation of the body-fixed displacement partial, without the
//! (much smaller) time variation of the displacement partial itself.
class CartesianPartialWrtDisplacementLoveNumber : public CartesianStatePartial
{
public:
    /*!
     *  \param bodyRotationModel Rotation model of the body
     *  \param deformationModels Basic tidal deformation models of the body
     *  \param degrees Degree of each Love or Shida number
     *  \param useShidaNumber Whether the partial is w.r.t. the Shida number l_n (true) or the Love number h_n (false)
     *  \param nominalPosition Nominal body-fixed position of the point
     */
    CartesianPartialWrtDisplacementLoveNumber(
            const std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel,
            const std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > >& deformationModels,
            const std::vector< int >& degrees,
            const bool useShidaNumber,
            const Eigen::Vector3d& nominalPosition ):
        bodyRotationModel_( bodyRotationModel ), deformationModels_( deformationModels ), degrees_( degrees ),
        unitLoveNumber_( useShidaNumber ? 0.0 : 1.0 ), unitShidaNumber_( useShidaNumber ? 1.0 : 0.0 ), nominalPosition_( nominalPosition )
    {}

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfPosition( const Eigen::Vector6d& state, const double time )
    {
        return Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) ) * getBodyFixedPartial( time );
    }

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfVelocity( const Eigen::Vector6d& state, const double time )
    {
        return bodyRotationModel_->getDerivativeOfRotationToBaseFrame( time ) * getBodyFixedPartial( time );
    }

private:
    //! Partial of the body-fixed displacement, with one column per degree, in the order of degrees_
    Eigen::Matrix< double, 3, Eigen::Dynamic > getBodyFixedPartial( const double time )
    {
        Eigen::Matrix< double, 3, Eigen::Dynamic > partial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, degrees_.size( ) );
        const Eigen::Vector3d siteUnitVector = nominalPosition_.normalized( );
        for( const auto& model : deformationModels_ )
        {
            const Eigen::Vector3d deformedBodyPosition = model->getDeformedBodyStateFunction( )( time ).segment( 0, 3 );
            const Eigen::Quaterniond deformedBodyRotation = model->getDeformedBodyRotationFunction( )( time );
            const double deformedBodyGravitationalParameter = model->getGravitionalParameterOfDeformedBody( )( );
            const double equatorialRadius = model->getDeformedBodyEquatorialRadius( );
            const std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyStateFunctions =
                    model->getDeformingBodyStateFunctions( );
            const std::vector< std::function< double( ) > > deformingBodyGravitationalParameters =
                    model->getGravitionalParametersOfDeformingBodies( );
            for( unsigned int j = 0; j < deformingBodyStateFunctions.size( ); j++ )
            {
                const Eigen::Vector3d relativeBodyPosition =
                        deformedBodyRotation * ( deformingBodyStateFunctions.at( j )( time ).segment( 0, 3 ) - deformedBodyPosition );
                const double gravitationalParameterRatio =
                        deformingBodyGravitationalParameters.at( j )( ) / deformedBodyGravitationalParameter;
                for( unsigned int i = 0; i < degrees_.size( ); i++ )
                {
                    switch( degrees_.at( i ) )
                    {
                        case 2:
                            partial.col( i ) += basic_astrodynamics::calculateDegreeTwoBasicTidalDisplacement( gravitationalParameterRatio,
                                                                                                               siteUnitVector,
                                                                                                               relativeBodyPosition,
                                                                                                               equatorialRadius,
                                                                                                               unitLoveNumber_,
                                                                                                               unitShidaNumber_ );
                            break;
                        case 3:
                            partial.col( i ) +=
                                    basic_astrodynamics::calculateDegreeThreeBasicTidalDisplacement( gravitationalParameterRatio,
                                                                                                     siteUnitVector,
                                                                                                     relativeBodyPosition,
                                                                                                     equatorialRadius,
                                                                                                     unitLoveNumber_,
                                                                                                     unitShidaNumber_ );
                            break;
                        default:
                            throw std::runtime_error(
                                    "Error, partials of basic tidal displacements of degrees other than 2 or 3 "
                                    "not implemented" );
                    }
                }
            }
        }
        return partial;
    }

    std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel_;

    std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > > deformationModels_;

    std::vector< int > degrees_;

    //! Love and Shida number for which the displacement is the partial: (1, 0) for h_n, (0, 1) for l_n
    double unitLoveNumber_;

    double unitShidaNumber_;

    Eigen::Vector3d nominalPosition_;
};

//! Derived class for scaling three-dimensional position partial to position observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to position observable partial. Although the implementation
 *  is trivial for non-relativistic reference frames, it is included in the architecture pending future implementation
 *  of more rigorous reference frames.
 */
class PositionObservationScaling : public DirectPositionPartialScaling< 3 >
{
public:
    PositionObservationScaling( ): DirectPositionPartialScaling< 3 >( observation_models::position_observable ) {}

    //! Destructor
    ~PositionObservationScaling( ) {}

    //! Update the scaling object to the current times and states (no functionality needed).
    /*!
     *  Update the scaling object to the current times and states (no functionality needed).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation,
                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
    {}

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 3, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return Eigen::Matrix3d::Identity( );
    }

    Eigen::Matrix< double, 3, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return getPositionScalingFactor( linkEndType );
    }

    virtual Eigen::Matrix< double, 3, 1 > getLightTimePartialScalingFactor( )
    {
        throw std::runtime_error( "Error when calculating position partial scaling factor; term non-existent for position observables" );
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return observation_models::observed_body;
    }

private:
};

//! Derived class for scaling three-dimensional position partial to relative position observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to relative position observable partial. Although the implementation
 *  is trivial for non-relativistic reference frames, it is included in the architecture pending future implementation
 *  of more rigorous reference frames.
 */
class RelativePositionObservationScaling : public DirectPositionPartialScaling< 3 >
{
public:
    RelativePositionObservationScaling( ): DirectPositionPartialScaling< 3 >( observation_models::relative_position_observable ) {}

    //! Destructor
    ~RelativePositionObservationScaling( ) {}

    //! Update the scaling object to the current times and states (no functionality needed).
    /*!
     *  Update the scaling object to the current times and states (no functionality needed).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation,
                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
    {
        currentLinkEndType_ = fixedLinkEnd;
    }

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 3, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        if( linkEndType == observation_models::observed_body )
        {
            return Eigen::Matrix3d::Identity( );
        }
        else if( linkEndType == observation_models::observer )
        {
            return -Eigen::Matrix3d::Identity( );
        }
        else
        {
            throw std::runtime_error( "Error when getting relative position scaling factor, did not recognize reference link end " +
                                      observation_models::getLinkEndTypeString( linkEndType ) );
        }
    }

    Eigen::Matrix< double, 3, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return getPositionScalingFactor( linkEndType );
    }

    virtual Eigen::Matrix< double, 3, 1 > getLightTimePartialScalingFactor( )
    {
        throw std::runtime_error( "Error when calculating position partial scaling factor; term non-existent for position observables" );
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:
    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;
};

//! Derived class for scaling three-dimensional velocity partial to velocity observable partial
/*!
 *  Derived class for scaling three-dimensional velocity partial to velocity observable partial. Although the implementation
 *  is trivial for non-relativistic reference frames, it is included in the architecture pending future implementation
 *  of more rigorous reference frames.
 */
class VelocityObservationScaling : public DirectPositionPartialScaling< 3 >
{
public:
    VelocityObservationScaling( ): DirectPositionPartialScaling< 3 >( observation_models::velocity_observable ) {}

    //! Destructor
    ~VelocityObservationScaling( ) {}

    //! Update the scaling object to the current times and states (no functionality needed).
    /*!
     *  Update the scaling object to the current times and states (no functionality needed).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation,
                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
    {}

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Velocity partial scaling factor at current link end
     */
    Eigen::Matrix< double, 3, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return Eigen::Matrix3d::Zero( );
    }

    Eigen::Matrix< double, 3, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return getPositionScalingFactor( linkEndType );
    }

    Eigen::Matrix< double, 3, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return Eigen::Matrix3d::Identity( );
    }

    Eigen::Matrix< double, 3, 3 > getFixedTimeVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return getVelocityScalingFactor( linkEndType );
    }

    virtual Eigen::Matrix< double, 3, 1 > getLightTimePartialScalingFactor( )
    {
        throw std::runtime_error( "Error when calculating velocity partial scaling factor; term non-existent for position observables" );
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return observation_models::observed_body;
    }

    virtual bool isVelocityScalingNonZero( const observation_models::LinkEndType linkEndType )
    {
        return true;
    }

private:
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_POSITIONPARTIALS_H
