#ifndef TUDAT_VELOCITYPARTIALS_H
#define TUDAT_VELOCITYPARTIALS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/function.hpp>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Base class for calculating the partial of an inertial velocity wrt a parameter.
/*!
 *  Base class for calculating the partial of an inertial velocity wrt a parameter. A derived class is implemented for
 *  each (type of) estimatable parameter. A separate instance of the class must be made for each distinct state.
 *  Note that partials wrt parameters describing a property of a rotation matrix from a local to the inertial frame is
 *  implemented in the RotationMatrixPartial class, which is
 *  then used by the VelocityPartialWrtRotationMatrixParameter derived class of this class
 */
class VelocityPartial
{
public:

    //! Destructor.
    virtual ~VelocityPartial( ){ }

    //! Pure virtual base class function for determining partial at current time and body state.
    /*!
     *  Pure virtual base class function for determining partial at current time and body state wrt parameter of specific
     *  implemented derived. class.
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point velocity wrt parameter (with specific parameter determined by derived class implementation).
     */
    virtual Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const Eigen::Vector6d& state, const double time ) = 0;
};

//! Class to compute the partial derivative of the three-dimensional velocity of a body w.r.t. to inertial three-dimensional
//! velocity of this body
class VelocityPartialWrtVelocity: public VelocityPartial
{
public:

    //! Constructor
    VelocityPartialWrtVelocity( ){ }

    //! Destructor
    ~VelocityPartialWrtVelocity( ){ }

    //! Function for determining partial at current time and body state.
    /*!
     *  Function for determining partial at current time and body state wrt inertial three-dimensional velocity
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point velocity wrt velocity
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const Eigen::Vector6d& state,
            const double time )
    {
        return Eigen::Matrix3d::Identity( );
    }
};

//! Class to compute the partial derivative of the three-dimensional velocity of a body w.r.t. to a property of a rotation
//! matrix to/from a body-fixed frame.
class VelocityPartialWrtRotationMatrixParameter: public VelocityPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param rotationMatrixPartialObject Object to compute the associated partial of a rotation matrix
     * \param velocityFunctionInLocalFrame Function returning the body-fixed velocity of the point at which the partial
     * is to be computed.
     */
    VelocityPartialWrtRotationMatrixParameter(
            const boost::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject,
            const boost::function< Eigen::Vector3d( const double ) > velocityFunctionInLocalFrame ):
        rotationMatrixPartialObject_( rotationMatrixPartialObject ),
        velocityFunctionInLocalFrame_( velocityFunctionInLocalFrame ){ }

    //! Destructor
    ~VelocityPartialWrtRotationMatrixParameter( ){ }

    //! Function for determining partial at current time and body state.
    /*!
     *  Function for determining partial at current time and body state wrt rotation property
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point velocity wrt rotation propert.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const Eigen::Vector6d& state,
            const double time )
    {
        throw std::runtime_error( "Error, velocity w.r.t. rotation matrix parameter not yet implemented" );
        //return rotationMatrixPartialObject_->calculatePartialOfRotatedVector( time, positionFunctionInLocalFrame_( time ) );
    }
private:

    //! Object to compute the associated partial of a rotation matrix
    boost::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject_;

    //! Function returning the body-fixed velocity of the point at which the partial is to be computed.
    boost::function< Eigen::Vector3d( const double ) > velocityFunctionInLocalFrame_;
};

//! Class to compute the partial derivative of the three-dimensional velocity of a body w.r.t. to body-fixed velocity of
//! some reference point (e.g. ground station)
class VelocityPartialWrtBodyFixedState: public VelocityPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyRotationModel Object to compute the rotation to/from the body-fixed frame.
     */
    VelocityPartialWrtBodyFixedState( const boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel ):
        bodyRotationModel_( bodyRotationModel ){ }

    //! Destructor
    ~VelocityPartialWrtBodyFixedState( ){ }

    //! Function for determining partial at current time and body state.
    /*!
     *  Function for determining partial at current time and body state wrt body-fixed point velocity
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point velocity wrt body-fixed point velocity
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const Eigen::Vector6d& state,
            const double time )
    {
        Eigen::Matrix3d currentRotationToBaseFrame = Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) );
        statePartial_.block( 0, 0, 3, 3 ) = currentRotationToBaseFrame;
        statePartial_.block( 0, 3, 3, 3 ) =
                -currentRotationToBaseFrame * bodyRotationModel_->getRotationalVelocityVectorInTargetFrame( time ) *
                currentRotationToBaseFrame;

        return statePartial_;
    }

private:

    Eigen::Matrix< double, 3, 6 > statePartial_;

    //! Object to compute the rotation to/from the body-fixed frame.
    boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel_;
};

}

}

#endif // TUDAT_VELOCITYPARTIALS_H
