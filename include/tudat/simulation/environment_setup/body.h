/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_BODY_H
#define TUDAT_BODY_H

#include <map>
#include <vector>

#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include "tudat/astro/basic_astro/polyhedronFuntions.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/propulsion/thrustGuidance.h"
//#include "tudat/astro/reference_frames/dependentOrientationCalculator.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/numericalDerivative.h"

namespace tudat {

namespace simulation_setup {

//! Base class used for the determination of the inertial state of a Body's ephemeris origin
/*!
 *  Base class used for the determination of the inertial state of a Body's ephemeris origin. This base class is used
 *  to provide an untemplated interface class through which to call the base frame state. The state may be defined
 *  in a templated manner in the derived class.
 */
class BaseStateInterface {
public:
    //! Constructor
    /*!
     * Constructor
     * \param baseFrameId Name of frame origin for which inertial state is computed by this class
     */
    BaseStateInterface(
            const std::string baseFrameId) : baseFrameId_(baseFrameId) {}

    //! Destructor
    virtual ~BaseStateInterface() {}

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    template<typename OutputTimeType, typename OutputStateScalarType>
    Eigen::Matrix<OutputStateScalarType, 6, 1> getBaseFrameState(
            const OutputTimeType time);

protected:
    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix<double, 6, 1> getBaseFrameDoubleState(const double time) = 0;

    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double long state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix<long double, 6, 1> getBaseFrameLongDoubleState(const double time) = 0;

    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix<double, 6, 1> getBaseFrameDoubleState(const Time &time) = 0;

    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and long double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix<long double, 6, 1> getBaseFrameLongDoubleState(const Time &time) = 0;

    //! Name of frame origin for which inertial state is computed by this class
    std::string baseFrameId_;
};

//! Class used for the determination of the inertial state of a Body's ephemeris origin
template<typename TimeType, typename StateScalarType>
class BaseStateInterfaceImplementation : public BaseStateInterface {
public:
    //! Constructor
    /*!
     * Constructor
     * \param baseFrameId Name of frame origin for which inertial state is computed by this class
     * \param stateFunction Function returning frame's inertial state as a function of time.
     * \param subtractStateFunction Boolean denoting whether to subtract or add the state function (i.e. whether to multiply
     * result of stateFunction by -1).
     */
    BaseStateInterfaceImplementation(
            const std::string baseFrameId,
            const std::function<Eigen::Matrix<StateScalarType, 6, 1>(const TimeType)> stateFunction,
            const bool subtractStateFunction = 0) : BaseStateInterface(baseFrameId),
        stateFunction_(stateFunction), stateMultiplier_((subtractStateFunction == 0) ? 1.0 : -1.0) {}

    //! Destructor
    ~BaseStateInterfaceImplementation() {}

protected:
    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix<double, 6, 1> getBaseFrameDoubleState(const double time) {
        return static_cast<double>(stateMultiplier_) * std::move(stateFunction_(time)).template cast<double>();
    }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double long state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix<long double, 6, 1> getBaseFrameLongDoubleState(const double time)
    {
        return static_cast<long double>(stateMultiplier_) * stateFunction_(time).template cast<long double>();
    }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix<double, 6, 1> getBaseFrameDoubleState(const Time &time)
    {
        return static_cast<double>(stateMultiplier_) * stateFunction_(time).template cast<double>();
    }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and long double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix<long double, 6, 1> getBaseFrameLongDoubleState(const Time &time)
    {
        return static_cast<long double>(stateMultiplier_) * std::move(stateFunction_(time)).template cast<long double>();
    }

private:
    //! Function returning frame's inertial state as a function of time.
    std::function<Eigen::Matrix<StateScalarType, 6, 1>(const TimeType)> stateFunction_;

    //! Value (1 or -1) by which to multiply the state returned by stateFunction_.
    int stateMultiplier_;
};

class RigidBodyProperties {
  public:
    RigidBodyProperties()
        : currentMass_(TUDAT_NAN),
          currentCenterOfMass_(Eigen::Vector3d::Constant(TUDAT_NAN)),
          currentInertiaTensor_(Eigen::Matrix3d::Constant(TUDAT_NAN)),
          isBodyInPropagation_(false),
          isMassComputed_(false),
          isComComputed_(false),
          isInertiaTensorComputed_(false){}

    virtual ~RigidBodyProperties() {}

    void update(const double currentTime)
    {
        updateMass( currentTime );
        updateMassDistribution( currentTime );
    }

    virtual void updateMass(const double currentTime) = 0;

    virtual void updateMassDistribution(const double currentTime) = 0;

    virtual void resetCurrentTime( )
    {
        isMassComputed_ = false;
        isComComputed_ = false;
        isInertiaTensorComputed_ = false;
    }

    double getCurrentMass( )
    {
        if( !isMassComputed_ )
        {
            throw std::runtime_error( "Error when retrieveing mass, mass is not computed/defined." );
        }
        return currentMass_;
    }

    Eigen::Vector3d getCurrentCenterOfMass( )
    {
        if( !isComComputed_ )
        {
            throw std::runtime_error( "Error when retrieving center of mass, center of mass is not computed/defined." );
        }
        return currentCenterOfMass_;
    }

    Eigen::Matrix3d getCurrentInertiaTensor( )
    {
        if( !isInertiaTensorComputed_ )
        {
            throw std::runtime_error( "Error when retrieving inertia tensor, inertia tensor is not computed/defined." );
        }
        return currentInertiaTensor_;
    }

    virtual void setCurrentMass( const double currentMass ) = 0;

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation )
    {
        isBodyInPropagation_ = isBodyInPropagation;
    }


protected:
    
    double currentMass_;

    Eigen::Vector3d currentCenterOfMass_;

    Eigen::Matrix3d currentInertiaTensor_;

    bool isBodyInPropagation_;

    bool isMassComputed_;

    bool isComComputed_;

    bool isInertiaTensorComputed_;
};

class TimeDependentRigidBodyProperties: public RigidBodyProperties
{
public:
    TimeDependentRigidBodyProperties(
            const std::function< double( const double ) > massFunction,
            const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction = nullptr,
            const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction = nullptr ):
            massFunction_( massFunction ),
            centerOfMassFunction_( centerOfMassFunction ),
            inertiaTensorFunction_( inertiaTensorFunction )
    {
    }

    TimeDependentRigidBodyProperties(
            const double mass,
            const Eigen::Vector3d& centerOfMass = Eigen::Vector3d::Constant( TUDAT_NAN ),
            const Eigen::Matrix3d& inertiaTensor = Eigen::Matrix3d::Constant( TUDAT_NAN ) ):
            massFunction_( nullptr ), centerOfMassFunction_( nullptr ), inertiaTensorFunction_( nullptr )
    {
        currentMass_ = mass;
        isMassComputed_ = true;

        if( !centerOfMass.hasNaN( ) )
        {
            currentCenterOfMass_ = centerOfMass;
            isComComputed_ = true;
        }

        if( !inertiaTensor.hasNaN( ) )
        {
            currentInertiaTensor_ = inertiaTensor;
            isInertiaTensorComputed_ = true;
        }
    }


    virtual ~TimeDependentRigidBodyProperties( ){ }

    virtual void resetCurrentTime( )
    {
        if( massFunction_ != nullptr )
        {
            isMassComputed_ = false;
        }

        if( centerOfMassFunction_ != nullptr )
        {
            isComComputed_ = false;
        }

        if( inertiaTensorFunction_ != nullptr )
        {
            isInertiaTensorComputed_ = false;
        }
    }

    std::function< double( const double ) > getMassFunction( )
    {
        return massFunction_;
    }

    void setMassFunction( const std::function< double( const double ) > massFunction )
    {
        massFunction_ = massFunction;
    }

    virtual void setCurrentMass( const double currentMass )
    {
        currentMass_ = currentMass;
        isMassComputed_ = true;
    }

    void setInertiaTensor( const Eigen::Matrix3d& inertiaTensor )
    {
        currentInertiaTensor_ = inertiaTensor;
        isInertiaTensorComputed_ = true;
    }

    virtual void updateMass( const double currentTime )
    {
        if( massFunction_ != nullptr && ( !isMassComputed_ || !isBodyInPropagation_ ) )
        {
            this->currentMass_ = massFunction_( currentTime );
            isMassComputed_ = true;
        }
    }

    virtual void updateMassDistribution( const double currentTime )
    {
        if( centerOfMassFunction_ != nullptr && ( !isComComputed_ || !isBodyInPropagation_ )  )
        {
            this->currentCenterOfMass_ = centerOfMassFunction_( currentTime );
            isComComputed_ = true;
        }

        if( inertiaTensorFunction_ != nullptr && ( !isInertiaTensorComputed_ || !isBodyInPropagation_ ) )
        {
            this->currentInertiaTensor_ = inertiaTensorFunction_( currentTime );
            isInertiaTensorComputed_ = true;
        }
    }

protected:

    std::function< double( const double ) > massFunction_;

    std::function< Eigen::Vector3d( const double ) > centerOfMassFunction_;

    std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction_;

};

class MassDependentRigidBodyProperties: public RigidBodyProperties
{
public:
    MassDependentRigidBodyProperties(
            const double currentMass,
            const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
            const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction ):
            centerOfMassFunction_( centerOfMassFunction ),
            inertiaTensorFunction_( inertiaTensorFunction )
    {
        setCurrentMass( currentMass );
        updateMassDistribution( TUDAT_NAN );
    }

    virtual ~MassDependentRigidBodyProperties( ){ }

    virtual void updateMass( const double currentTime )
    {
    }

    virtual void updateMassDistribution( const double currentTime )
    {
        if( centerOfMassFunction_ != nullptr && ( !isComComputed_ || !isBodyInPropagation_ ) )
        {
            currentCenterOfMass_ = centerOfMassFunction_( currentMass_ );
            isComComputed_ = true;
        }

        if( inertiaTensorFunction_ != nullptr && ( !isInertiaTensorComputed_ || !isBodyInPropagation_ ) )
        {
            currentInertiaTensor_ = inertiaTensorFunction_( currentMass_ );
            isInertiaTensorComputed_ = true;
        }
    }

    virtual void setCurrentMass( const double currentMass )
    {
        currentMass_ = currentMass;
        isMassComputed_ = true;
    }

protected:

    std::function< Eigen::Vector3d( const double ) > centerOfMassFunction_;

    std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction_;

};


class FromGravityFieldRigidBodyProperties: public RigidBodyProperties
{
public:
    FromGravityFieldRigidBodyProperties(
            const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel ):
            gravityFieldModel_( gravityFieldModel )
    {
        currentMass_ = gravityFieldModel_->getGravitationalParameter() / physical_constants::GRAVITATIONAL_CONSTANT;
        currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
        currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );

        isMassComputed_ = true;
        isComComputed_ = true;
        isInertiaTensorComputed_ = true;

        if( std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel ) != nullptr )
        {
            modelIsTimeDependent_  = true;
        }
        else
        {
            modelIsTimeDependent_ = false;
        }
    }

    virtual ~FromGravityFieldRigidBodyProperties( ){ }

    virtual void resetCurrentTime( )
    {
        if( modelIsTimeDependent_ )
        {
            isMassComputed_ = false;
            isComComputed_ = false;
            isInertiaTensorComputed_ = false;
        }
    }

    virtual void updateMass( const double currentTime )
    {
        if( ( modelIsTimeDependent_ && !isMassComputed_ ) || !isBodyInPropagation_ )
        {
            currentMass_ = gravityFieldModel_->getGravitationalParameter() /
                           physical_constants::GRAVITATIONAL_CONSTANT;
            isMassComputed_ = true;
        }
    }

    virtual void updateMassDistribution( const double currentTime )
    {
        if( ( modelIsTimeDependent_ && ( !isComComputed_ || !isInertiaTensorComputed_ ) ) || !isBodyInPropagation_ )
        {
            currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass();
            currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor();
            isComComputed_ = true;
            isInertiaTensorComputed_ = true;
        }
    }

    virtual void setCurrentMass( const double currentMass )
    {
        throw std::runtime_error( "Error when resetting body mass; mass cannot be reset for bodies with a gravity field. Reset the gravity field's gravitational parameter instead." );
    }

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation )
    {
        currentMass_ = gravityFieldModel_->getGravitationalParameter() / physical_constants::GRAVITATIONAL_CONSTANT;
        currentCenterOfMass_ = gravityFieldModel_->getCenterOfMass( );
        currentInertiaTensor_ = gravityFieldModel_->getInertiaTensor( );

        isBodyInPropagation_ = isBodyInPropagation;
    }

protected:

    const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

    bool modelIsTimeDependent_;


};


//! Body class representing the properties of a celestial body (natural or artificial).
/*!
 *  Body class representing the properties of a celestial body (natural or artificial). By storing
 *  all properties of bodies (ephemeris, rotation, gravity, etc.) in a set of body objects,
 *  the simulation environment can be defined in a clear and modular way. To create body
 *  objects, the createBodies.h function provides a range of functionality. The
 *  createAccelerationModels.h file provides functions to use body objects to create acceleration
 *  objects.
 */
class Body {
public:
    //! Constructor for a body
    /*!
     * Constructor for a body, sets current state (with zero default value).
     * \param state Current state of body at initialization (default = zeroes).
     */
    Body( const Eigen::Vector6d& state =
            Eigen::Vector6d::Zero( ) )
        : bodyIsGlobalFrameOrigin_( -1 ), currentState_( state ), timeOfCurrentState_( TUDAT_NAN ),
          ephemerisFrameToBaseFrame_( std::make_shared< BaseStateInterfaceImplementation< double, double > >(
                                          "", [ = ]( const double ){ return Eigen::Vector6d::Zero( ); } ) ),
          currentRotationToLocalFrame_( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
          currentRotationToGlobalFrame_( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
          currentRotationToLocalFrameDerivative_( Eigen::Matrix3d::Zero( ) ),
          currentAngularVelocityVectorInGlobalFrame_( Eigen::Vector3d::Zero( ) ),
          currentAngularVelocityVectorInLocalFrame_( Eigen::Vector3d::Zero( ) ),
          bodyName_( "unnamed_body" )
    {
        currentLongState_ = currentState_.cast< long double >( );
        isStateSet_ = false;
        isRotationSet_ = false;
    }

    //! Function to retrieve the class returning the state of this body's ephemeris origin w.r.t. the global origin
    /*!
     * Function to retrieve the class returning the state of this body's ephemeris origin w.r.t. the global origin
     * \return Class returning the state of this body's ephemeris origin w.r.t. the global origin
     */
    std::shared_ptr<BaseStateInterface> getEphemerisFrameToBaseFrame() {
        return ephemerisFrameToBaseFrame_;
    }

    //! Function to set the class returning the state of this body's ephemeris origin w.r.t. the global origin
    /*!
     * Function to set the class returning the state of this body's ephemeris origin w.r.t. the global origin
     * \param ephemerisFrameToBaseFrame Class returning the state of this body's ephemeris origin w.r.t. the global origin
     */
    void setEphemerisFrameToBaseFrame(const std::shared_ptr<BaseStateInterface> ephemerisFrameToBaseFrame) {
        ephemerisFrameToBaseFrame_ = ephemerisFrameToBaseFrame;
    }

    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    Eigen::Vector6d getState( )
    {
        if( !isStateSet_ )
        {
            std::vector< double > test;
            test.at( 0 );
            throw std::runtime_error( "Error when retrieving state from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentState_;
        }
    }

    //! Set current state of body manually
    /*!
     * Set current state of body manually, which must be in the global frame. Note that this
     * function does not set the currentLongState_, use the setLongState when needing the use of the
     * long precision current state.
     * \param state Current state of the body that is set.
     */
    void setState(const Eigen::Vector6d &state)
    {
        currentState_ = state;
        isStateSet_ = true;
    }

    //! Set current state of body manually in long double precision.
    /*!
     * Set current state of body manually in long double precision. State must be in the global
     * frame.  Note that this function sets both the currentState_ and currentLongState_ variables
     * (currentLongState_ directly and currentState_ by casting the input to double entries).
     * \param longState Current state of the body that is set, in long double precision.
     */
    void setLongState(const Eigen::Matrix<long double, 6, 1> &longState) {
        currentLongState_ = longState;
        currentState_ = longState.cast<double>();
        isStateSet_ = true;

    }

    //! Templated function to set the state manually.
    /*!
     * Templated function to set the state manually, calls either setState or setLongState function.
     * \param state Current state of the body that is set, with StateScalarType precision.
     */
    template<typename StateScalarType>
    void setTemplatedState(const Eigen::Matrix<StateScalarType, 6, 1> &state);

    //! Templated function to set the current state of the body from its ephemeris and
    //! global-to-ephemeris-frame function.
    /*!
     * Templated function to set the current state of the body from its ephemeris and
     * global-to-ephemeris-frame function. It sets both the currentState_ and currentLongState_ variables. F
     * FUndamental coputation is done on state with StateScalarType precision as a function of TimeType time
     * \param time Time at which the global state is to be set.
     */
    template<typename StateScalarType = double, typename TimeType = double>
    void setStateFromEphemeris(const TimeType &time)
    {
        if (!(static_cast<Time>(time) == timeOfCurrentState_))
        {
            if( bodyEphemeris_ == nullptr )
            {
                throw std::runtime_error( "Error when requesting state from ephemeris of body " + bodyName_ + ", body has no ephemeris" );
            }
            // If body is not global frame origin, set state.
            if (bodyIsGlobalFrameOrigin_ == 0)
            {
                if (sizeof(StateScalarType) == 8)
                {
                    currentState_ =
                            (bodyEphemeris_->getTemplatedStateFromEphemeris<StateScalarType, TimeType>(time) + ephemerisFrameToBaseFrame_->getBaseFrameState<TimeType, StateScalarType>(time)).template cast<double>();
                    currentLongState_ = currentState_.template cast<long double>();
                }
                else
                {
                    currentLongState_ =
                            (bodyEphemeris_->getTemplatedStateFromEphemeris<StateScalarType, TimeType>(time) + ephemerisFrameToBaseFrame_->getBaseFrameState<TimeType, StateScalarType>(time) ).template cast<long double>();
                    currentState_ = currentLongState_.template cast<double>();
                }
            }
            // If body is global frame origin, set state to zeroes, and barycentric state value.
            else if (bodyIsGlobalFrameOrigin_ == 1)
            {
                currentState_.setZero();
                currentLongState_.setZero();

                if (sizeof(StateScalarType) == 8)
                {
                    currentBarycentricState_ =
                            ephemerisFrameToBaseFrame_->getBaseFrameState<TimeType, StateScalarType>(time).template cast<double>();
                    currentBarycentricLongState_ = currentBarycentricState_.template cast<long double>();
                }
                else
                {
                    currentBarycentricLongState_ =
                            ephemerisFrameToBaseFrame_->getBaseFrameState<TimeType, StateScalarType>(time).template cast<long double>();
                    currentBarycentricState_ = currentBarycentricLongState_.template cast<double>();
                }
            }
            else
            {
                throw std::runtime_error("Error when setting body state, global origin not yet defined.");
            }

            timeOfCurrentState_ = static_cast<TimeType>(time);
        }
        isStateSet_ = true;
    }

    //    extern template void setStateFromEphemeris< double, double >( const double& time );

    //! Templated function to get the current state of the body from its ephemeris and
    //! global-to-ephemeris-frame function.
    /*!
     * Templated function to get the current state of the body from its ephemeris and
     * global-to-ephemeris-frame function.  It calls the setStateFromEphemeris state, resetting the currentState_ /
     * currentLongState_ variables, and returning the state with the requested precision
     * \param time Time at which to evaluate states.
     * \return State at requested time
     */
    template<typename StateScalarType = double, typename TimeType = double>
    Eigen::Matrix<StateScalarType, 6, 1> getStateInBaseFrameFromEphemeris(const TimeType time)
    {
        setStateFromEphemeris<StateScalarType, TimeType>(time);
        if (sizeof(StateScalarType) == 8)
        {
            return currentState_.template cast<StateScalarType>();
        }
        else
        {
            return currentLongState_.template cast<StateScalarType>();
        }
    }

    //! Templated function to get the current berycentric state of the body from its ephemeris andcglobal-to-ephemeris-frame
    //! function.
    /*!
     * Templated function to get the current berycentric state of the body from its ephemeris andcglobal-to-ephemeris-frame
     * function. It calls the setStateFromEphemeris state, resetting the currentBarycentricState_ /
     * currentBarycentricLongState_ variables, and returning the state with the requested precision. This function can ONLY be
     * called if this body is the global frame origin, otherwise an exception is thrown
     * \param time Time at which to evaluate states.
     * \return Barycentric State at requested time
     */
    template<typename StateScalarType = double, typename TimeType = double>
    Eigen::Matrix<StateScalarType, 6, 1> getGlobalFrameOriginBarycentricStateFromEphemeris(const TimeType time) {
        if (bodyIsGlobalFrameOrigin_ != 1) {
            throw std::runtime_error("Error, calling global frame origin barycentric state on body that is not global frame origin");
        }

        setStateFromEphemeris<StateScalarType, TimeType>(time);

        if (sizeof(StateScalarType) == 8) {
            return currentBarycentricState_.template cast<StateScalarType>();
        } else {
            return currentBarycentricLongState_.template cast<StateScalarType>();
        }
    }

    //! Get current rotational state.
    /*!
     * Returns the internally stored current rotational state vector.
     * \return Current rotational state.
     */
    Eigen::Vector7d getRotationalStateVector() {
        Eigen::Vector7d rotationalStateVector;

        rotationalStateVector.segment( 0, 4 ) =
                linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( currentRotationToGlobalFrame_ ) );
        rotationalStateVector.segment( 4, 3 ) = currentAngularVelocityVectorInLocalFrame_;
        return rotationalStateVector;
    }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getPosition()
    {
        if( !isStateSet_ )
        {
            throw std::runtime_error( "Error when retrieving position from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentState_.segment(0, 3);
        }
    }

    void getPositionByReference( Eigen::Vector3d& position );

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getVelocity()
    {
        if( !isStateSet_ )
        {
            throw std::runtime_error( "Error when retrieving velociy from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentState_.segment(3, 3);
        }
    }

    //! Get current state, in long double precision
    /*!
     * Returns the internally stored current state vector, in long double precision
     * \return Current state, in long double precisio
     */
    Eigen::Matrix<long double, 6, 1> getLongState()
    {
        if( !isStateSet_ )
        {
            throw std::runtime_error( "Error when retrieving long state from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentLongState_;
        }
    }

    //! Get current position, in long double precision
    /*!
     * Returns the internally stored current position vector, in long double precision
     * \return Current position, in long double precision
     */
    Eigen::Matrix<long double, 3, 1> getLongPosition( )
    {
        if( !isStateSet_ )
        {
            throw std::runtime_error( "Error when retrieving long position from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentLongState_.segment(0, 3);
        }
    }

    //! Get current velocity, in long double precision.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity, in long double precision
     */
    Eigen::Matrix<long double, 3, 1> getLongVelocity( )
    {
        if( !isStateSet_ )
        {
            throw std::runtime_error( "Error when retrieving long velocity from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentLongState_.segment(3, 3);
        }
    }

    //! Templated function to retrieve the state.
    /*!
     * Templated function to retrieve the state, calls either getState or getLongState function.
     * \return  Current state of the body, with StateScalarType precision.
     */
    template<typename ScalarStateType>
    Eigen::Matrix<ScalarStateType, 6, 1> getTemplatedState();

    //! Function to set the rotation from global to body-fixed frame at given time
    /*!
     * Function to set the rotation from global to body-fixed frame at given time, using the
     * rotationalEphemeris_ member object
     * \param time Time at which the rotation is to be retrieved.
     */
    void setCurrentRotationToLocalFrameFromEphemeris( const double time )
    {
        if( rotationalEphemeris_!= nullptr )
        {
            currentRotationToLocalFrame_ = rotationalEphemeris_->getRotationToTargetFrame( time );
        }
//        else if( dependentOrientationCalculator_ != nullptr )
//        {
//            currentRotationToLocalFrame_ = dependentOrientationCalculator_->computeAndGetRotationToLocalFrame( time );
//        }
        else
        {
            throw std::runtime_error(
                        "Error, no rotation model found in Body::setCurrentRotationToLocalFrameFromEphemeris" );
        }
        currentRotationToGlobalFrame_ = currentRotationToLocalFrame_.inverse( );
        isRotationSet_ = true;
    }

    //! Function to set the full rotational state at given time
    /*!
     * Function to set the full rotational state at (rotation from global to body-fixed frame
     * rotation matrix derivative from global to body-fixed frame and angular velocity vector in the
     * global frame) at given time, using the rotationalEphemeris_ member object.
     * \param time Time at which the angular velocity vector in the global frame is to be retrieved.
     */
    template< typename TimeType >
    void setCurrentRotationalStateToLocalFrameFromEphemeris( const TimeType time )
    {
        if( rotationalEphemeris_ != nullptr )
        {
            rotationalEphemeris_->getFullRotationalQuantitiesToTargetFrameTemplated< TimeType >(
                        currentRotationToLocalFrame_, currentRotationToLocalFrameDerivative_,
                        currentAngularVelocityVectorInGlobalFrame_, time );
            currentAngularVelocityVectorInLocalFrame_ = currentRotationToLocalFrame_ * currentAngularVelocityVectorInGlobalFrame_;
        }
//        else if( dependentOrientationCalculator_ != nullptr )
//        {
//            currentRotationToLocalFrame_ = dependentOrientationCalculator_->computeAndGetRotationToLocalFrame( time );
//            currentRotationToLocalFrameDerivative_.setZero( );
//            currentAngularVelocityVectorInGlobalFrame_.setZero( );
//            currentAngularVelocityVectorInLocalFrame_.setZero( );
//        }
        else
        {
            throw std::runtime_error(
                        "Error, no rotationalEphemeris_ found in Body::setCurrentRotationalStateToLocalFrameFromEphemeris" );
        }
        currentRotationToGlobalFrame_ = currentRotationToLocalFrame_.inverse( );
        isRotationSet_ = true;

    }

    //! Function to set the full rotational state directly
    /*!
     * Function to set the full rotational state  directly (rotation from global to body-fixed frame
     * rotation matrix derivative from global to body-fixed frame and angular velocity vector in the
     * global frame) directly, by providing the current rotational state as input.
     * \param currentRotationalStateFromLocalToGlobalFrame Quaternion from body-fixed to propagation frame
     * (in vector form) and the body's angular velocity vector in body-fixed frame.
     */
  void setCurrentRotationalStateToLocalFrame(const Eigen::Vector7d currentRotationalStateFromLocalToGlobalFrame) {
    currentRotationToGlobalFrame_ =
        Eigen::Quaterniond(currentRotationalStateFromLocalToGlobalFrame(0),
                           currentRotationalStateFromLocalToGlobalFrame(1),
                           currentRotationalStateFromLocalToGlobalFrame(2),
                           currentRotationalStateFromLocalToGlobalFrame(3));

    currentRotationToGlobalFrame_.normalize();
    currentRotationToLocalFrame_ = currentRotationToGlobalFrame_.inverse();
    currentAngularVelocityVectorInGlobalFrame_ =
        currentRotationToGlobalFrame_ * currentRotationalStateFromLocalToGlobalFrame.block< 3, 1 >(4, 0);

    currentAngularVelocityVectorInLocalFrame_ = currentRotationalStateFromLocalToGlobalFrame.block< 3, 1 >(4, 0);

    Eigen::Matrix3d currentRotationMatrixToLocalFrame = (currentRotationToLocalFrame_).toRotationMatrix();
    currentRotationToLocalFrameDerivative_ = linear_algebra::getCrossProductMatrix(
                                                 currentRotationalStateFromLocalToGlobalFrame.block< 3, 1 >(4, 0 ))
        * currentRotationMatrixToLocalFrame;
    isRotationSet_ = true;

  }

    //! Get current rotation from body-fixed to inertial frame.
    /*!
     *  Get current rotation from body-fixed to inertial frame, as set from the rotationalEphemeris_
     *  by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameFromEphemeris function.  If body has no rotational ephemeris,
     *  an identity quaternion (no rotation) is returned.
     *  \return Current rotation from body-fixed to inertial frame.
     */
    Eigen::Quaterniond getCurrentRotationToGlobalFrame( )
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving rotation to global frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentRotationToGlobalFrame_;
        }
    }

    Eigen::Quaterniond& getCurrentRotationToGlobalFrameReference( )
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving rotation to global frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentRotationToGlobalFrame_;
        }
    }

    //! Get current rotation from inertial to body-fixed frame.
    /*!
     *  Get current rotation from inertial to body-fixed frame, as set from the rotationalEphemeris_
     *  by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameFromEphemeris function.  If body has no rotational ephemeris,
     *  an identity quaternion (no rotation) is returned.
     *  \return Current rotation from inertial to body-fixed frame.
     */
    Eigen::Quaterniond getCurrentRotationToLocalFrame()
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving rotation to local frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentRotationToLocalFrame_;
        }
    }

    Eigen::Matrix3d getCurrentRotationMatrixToGlobalFrame()
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving rotation to global frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return Eigen::Matrix3d( currentRotationToLocalFrame_.inverse() );
        }
    }

    Eigen::Matrix3d getCurrentRotationMatrixToLocalFrame()
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving rotation to local frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return Eigen::Matrix3d( currentRotationToLocalFrame_ );
        }
    }

    //! Get current rotational state.
    /*!
     *  Get current rotational state, expressed as a quaternion from global to body-fixed frame
     *  (in vector form) and the body's angular velocity vector in inertial frame.
     *  \return Current rotational state in quaternions and rotational velocity.
     */
    Eigen::Vector7d getCurrentRotationalState()
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving rotation from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return (Eigen::VectorXd(7) << linear_algebra::convertQuaternionToVectorFormat(getCurrentRotationToGlobalFrame()),
                    getCurrentAngularVelocityVectorInGlobalFrame())
                    .finished();
        }
    }

    //! Get current rotation matrix derivative from body-fixed to global frame.
    /*!
     *  Get current rotation matrix derivative from body-fixed frame to global, as set from the
     *  rotationalEphemeris_ by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameDerivativeFromEphemeris function. If body has no rotational
     *  ephemeris, an zero matrix (no rotation) is returned.
     *  \return Current otation matrix derivative from global to body-fixed frame.
     */
    Eigen::Matrix3d getCurrentRotationMatrixDerivativeToGlobalFrame( )
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving derivative of rotation to global frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else if( currentRotationToLocalFrameDerivative_.hasNaN( ) )
        {
            throw std::runtime_error( "Error when retrieving derivative of rotation to global frame from body " + bodyName_ + ", matrix is undefined" );
        }
        else
        {
            return currentRotationToLocalFrameDerivative_.transpose();
        }
    }

    //! Get current rotation matrix derivative from global to body-fixed frame.
    /*!
     *  Get current rotation matrix derivative from global to body-fixed frame, as set from the
     *  rotationalEphemeris_ by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameDerivativeFromEphemeris function. If body has no rotational
     *  ephemeris, an zero matrix (no rotation) is returned.
     *  \return Current otation matrix derivative from global to body-fixed frame.
     */
    Eigen::Matrix3d getCurrentRotationMatrixDerivativeToLocalFrame()
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving derivative of rotation to local frame from body " + bodyName_ + ", state of body is not yet defined" );
        }
        else if( currentRotationToLocalFrameDerivative_.hasNaN( ) )
        {
            throw std::runtime_error( "Error when retrieving derivative of rotation to local frame from body " + bodyName_ + ", matrix is undefined" );
        }
        else
        {
            return currentRotationToLocalFrameDerivative_;
        }
    }

    //! Get current angular velocity vector for body's rotation, expressed in the global frame.
    /*!
     *  Get current angular velocity vector for body's rotation, expressed in the global frame.
     *  \return Current angular velocity vector for body's rotation, expressed in the global frame.
     */
    Eigen::Vector3d getCurrentAngularVelocityVectorInGlobalFrame( )
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving angular velocioty of body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentAngularVelocityVectorInGlobalFrame_;
        }
    }

    //! Get current angular velocity vector for body's rotation, expressed in the local frame.
    /*!
     *  Get current angular velocity vector for body's rotation, expressed in the local frame.
     *  Transformation from the global to the local frame is done by rotating the vector with the
     *  current quaternion to local frame.
     *  \return Current angular velocity vector for body's rotation, expressed in the local frame.
     */
    Eigen::Vector3d getCurrentAngularVelocityVectorInLocalFrame( )
    {
        if( !isRotationSet_ )
        {
            throw std::runtime_error( "Error when retrieving angular velocioty of body " + bodyName_ + ", state of body is not yet defined" );
        }
        else
        {
            return currentAngularVelocityVectorInLocalFrame_;
        }
    }

    //! Function to set the ephemeris of the body.
    /*!
     *  Function to set the ephemeris of the body, which is used to represent the (a priori)
     *  state history of the body.
     *  \param bodyEphemeris New ephemeris of the body.
     */
    void setEphemeris( const std::shared_ptr< ephemerides::Ephemeris > bodyEphemeris )
    {
        bodyEphemeris_ = bodyEphemeris;
    }

    //! Function to set the gravity field of the body.
    /*!
     *  Function to set the gravity field of the body; input is also used to (re)set the mass of the
     *  body.
     *  \param gravityFieldModel New gravity field of the body.
     */
    void setGravityFieldModel(
            const std::shared_ptr<gravitation::GravityFieldModel> gravityFieldModel)
    {
        gravityFieldModel_ = gravityFieldModel;

        // Update current mass of body, provide warning
        if( massProperties_ != nullptr )
        {
            std::cerr << "Warning when settings gravity field model for body, mass interface already found: overrriding existing mass interface" << std::endl;
        }
        else
        {
            massProperties_ = std::make_shared< FromGravityFieldRigidBodyProperties >( gravityFieldModel );
        }

    }

    //! Function to set the atmosphere model of the body.
    /*!
     *  Function to set the atmosphere model of the body.
     *  \param atmosphereModel Atmosphere model of the body.
     */
    void setAtmosphereModel(
            const std::shared_ptr<aerodynamics::AtmosphereModel> atmosphereModel) {
        atmosphereModel_ = atmosphereModel;
    }

    //! Function to set the rotation model of the body.
    /*!
     *  Function to set the rotation model of the body.
     *  \param rotationalEphemeris Rotation model of the body.
     */
    void setRotationalEphemeris(
            const std::shared_ptr<ephemerides::RotationalEphemeris> rotationalEphemeris) {
//        if (dependentOrientationCalculator_ != nullptr) {
//            std::cerr << "Warning when setting rotational ephemeris, dependentOrientationCalculator_ already found, NOT setting closure" << std::endl;
//        }
        rotationalEphemeris_ = rotationalEphemeris;
    }

    //! Function to set the shape model of the body.
    /*!
     *  Function to set the shape model of the body.
     *  \param shapeModel Shape model of the body.
     */
    void setShapeModel(const std::shared_ptr<basic_astrodynamics::BodyShapeModel> shapeModel) {
        shapeModel_ = shapeModel;
    }

    //! Function to set the aerodynamic coefficient interface of the body.
    /*!
     *  Function to set the aerodynamic coefficient interface of the body.
     *  \param aerodynamicCoefficientInterface Aerodynamic coefficient interface of the body.
     */
    void setAerodynamicCoefficientInterface(
            const std::shared_ptr<aerodynamics::AerodynamicCoefficientInterface>
            aerodynamicCoefficientInterface)
    {
        aerodynamicCoefficientInterface_ = aerodynamicCoefficientInterface;
        std::vector< std::string > controlSurfaceList = aerodynamicCoefficientInterface->getControlSurfaceNames( );
        for( unsigned int i = 0; i < controlSurfaceList.size( ); i++ )
        {
            if ( getVehicleSystems( )->doesControlSurfaceExist( controlSurfaceList.at( i )) == 0 )
            {
                vehicleSystems_->setCurrentControlSurfaceDeflection( controlSurfaceList.at( i ), 0.0 );
            }
        }

        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( aerodynamicFlightConditions_ ) != nullptr )
        {
            std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( aerodynamicFlightConditions_ )->resetAerodynamicCoefficientInterface( aerodynamicCoefficientInterface_ );
        }
    }

    //! Function to set the body flight conditions
    /*!
     * Function to set the body flight conditions, which calculates current aerodynamic angles,
     * altitude, etc.
     * \param aerodynamicFlightConditions Body flight conditions
     */
    void setFlightConditions(
            const std::shared_ptr<aerodynamics::FlightConditions> aerodynamicFlightConditions) {
        aerodynamicFlightConditions_ = aerodynamicFlightConditions;

//        // If dependentOrientationCalculator_ object already exists, provide a warning and create closure between the two
//        if (dependentOrientationCalculator_ != nullptr) {
//            reference_frames::setAerodynamicDependentOrientationCalculatorClosure(
//                        dependentOrientationCalculator_, aerodynamicFlightConditions_->getAerodynamicAngleCalculator());
//        } else {
//            dependentOrientationCalculator_ = aerodynamicFlightConditions->getAerodynamicAngleCalculator();
//        }

        // Create closure between rotational ephemeris and aerodynamic angle calculator.
        if( rotationalEphemeris_ != nullptr && std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                    rotationalEphemeris_ ) == nullptr )
        {
            aerodynamicFlightConditions_->getAerodynamicAngleCalculator( )->setBodyFixedAngleInterface(
                        std::make_shared< reference_frames::FromGenericEphemerisAerodynamicAngleInterface >(
                            rotationalEphemeris_ ) );
        }
    }

    std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > > getBodyDeformationModels( )
    {
        return bodyDeformationModels_;
    }

    std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >& getBodyDeformationModelsReference( )
    {
        return bodyDeformationModels_;
    }

    void addBodyDeformationModel( const std::shared_ptr< basic_astrodynamics::BodyDeformationModel > deformationModel )
    {
        bodyDeformationModels_.push_back( deformationModel );
    }
    //! Function to set the radiation pressure interface of the body, for a single radiation source.
    /*!
     *  Function to set the radiation pressure interface of the body, for a single radiation source
     *  \param radiatingBody Name of body that is the source of the radiation.
     *  \param radiationPressureInterface Radiation pressure interface of the body.
     */
    void setRadiationPressureInterface(
            const std::string &radiatingBody,
            const std::shared_ptr<electromagnetism::RadiationPressureInterface>
            radiationPressureInterface) {
        // RP-OLD
        radiationPressureInterfaces_[radiatingBody] = radiationPressureInterface;
    }

    //! Function to set the radiation source model of the body.
    /*!
     *  Function to set the radiation source model of the body.
     *  \param radiationSourceModel Radiation source model of the body.
     */
    void setRadiationSourceModel(
            const std::shared_ptr<electromagnetism::RadiationSourceModel> radiationSourceModel)
    {
        radiationSourceModel_ = radiationSourceModel;
    }

    //! Function to set the radiation pressure target model of the body.
    /*!
     *  Function to set the radiation pressure target model of the body.
     *  \param radiationPressureTargetModel Radiation pressure target model of the body.
     */
    void setRadiationPressureTargetModel(
            const std::shared_ptr<electromagnetism::RadiationPressureTargetModel> radiationPressureTargetModel)
    {
        radiationPressureTargetModel_ = radiationPressureTargetModel;
    }

    //! Function to set object containing all variations in the gravity field of this body.
    /*!
     * Function to set object containing all variations in the gravity field of this body.
     * \param gravityFieldVariationSet Object containing all variations in the gravity field of this body.
     */
    void setGravityFieldVariationSet(
            const std::shared_ptr<gravitation::GravityFieldVariationsSet>
            gravityFieldVariationSet) {
        gravityFieldVariationSet_ = gravityFieldVariationSet;
    }

    //! Function to get the gravity field model of the body.
    /*!
     *  Function to get the gravity field model of the body.
     *  \return Gravity field model of the body.
     */
    std::shared_ptr<gravitation::GravityFieldModel> getGravityFieldModel() {
        return gravityFieldModel_;
    }

    double getGravitationalParameter( )
    {
        if( gravityFieldModel_ == nullptr )
        {
            throw std::runtime_error( "Error when retrieveing gravitational parameter from body " + bodyName_ +
                                      ", no gravity field model is defined" );
        }
        return gravityFieldModel_->getGravitationalParameter( );
    }

    //! Function to get the ephemeris of the body.
    /*!
     *  Function to get the ephemeris of the body.
     *  \return Ephemeris of the body.
     */
    std::shared_ptr<ephemerides::Ephemeris> getEphemeris() {
        return bodyEphemeris_;
    }

    //! Function to get the atmosphere model of the body.
    /*!
     *  Function to get the atmosphere model of the body.
     *  \return Atmosphere model of the body.
     */
    std::shared_ptr<aerodynamics::AtmosphereModel> getAtmosphereModel() {
        return atmosphereModel_;
    }

    //! Function to get the rotation model of the body.
    /*!
     *  Function to get the rotation model of the body.
     *  \return Rotation model of the body.
     */
    std::shared_ptr<ephemerides::RotationalEphemeris> getRotationalEphemeris() {
        return rotationalEphemeris_;
    }

//    //! Function to retrieve the model to compute the rotation of the body based on the current state of the environment.
//    /*!
//     * Function to retrieve the model to compute the rotation of the body based on the current state of the environment
//     * (model is only valid during propagation).
//     * \return Model to compute the rotation of the body based on the current state of the environment
//     */
//    std::shared_ptr<reference_frames::DependentOrientationCalculator> getDependentOrientationCalculator() {
//        return dependentOrientationCalculator_;
//    }

    //! Function to retrieve the shape model of body.
    /*!
     * Function to retrieve the shape model of body.
     * \return Shape model of body.
     */
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> getShapeModel() {
        return shapeModel_;
    }

    //! Function to retrieve the aerodynamic coefficient model of body.
    /*!
     * Function to retrieve the body aerodynamic coefficient model of body.
     * \return Aerodynamic coefficient model of body.
     */
    std::shared_ptr<aerodynamics::AerodynamicCoefficientInterface>
    getAerodynamicCoefficientInterface() {
        return aerodynamicCoefficientInterface_;
    }

    //! Function to retrieve the body flight conditions
    /*!
     * Function to retrieve the body flight conditions, which calculates current aerodynamic angles,
     * altitude, etc.
     * \return Body flight conditions
     */
    std::shared_ptr<aerodynamics::FlightConditions> getFlightConditions() {
        return aerodynamicFlightConditions_;
    }

    //! Function to retrieve the shape model of the body.
    /*!
     *  Function to retrieve the shape model of the body.
     *  \return Shape model of the body.
     */
    std::map<std::string, std::shared_ptr<electromagnetism::RadiationPressureInterface>>
    getRadiationPressureInterfaces() {
        // RP-OLD
        return radiationPressureInterfaces_;
    }

    //! Function to retrieve the radiation source model of the body.
    /*!
     *  Function to retrieve the radiation source model of the body.
     *  \return Radiation source model of the body.
     */
    const std::shared_ptr<electromagnetism::RadiationSourceModel> getRadiationSourceModel() const
    {
        return radiationSourceModel_;
    }

    //! Function to retrieve the radiation pressure target model of the body.
    /*!
     *  Function to retrieve the radiation pressure target model of the body.
     *  \return Radiation pressure target model of the body.
     */
    const std::shared_ptr<electromagnetism::RadiationPressureTargetModel> getRadiationPressureTargetModel() const
    {
        return radiationPressureTargetModel_;
    }

    //! Function to retrieve a single object describing variation in the gravity field of this body.
    /*!
     *  Function to retrieve a single object describing variation in the gravity field of this body.
     *  \param deformationType Type of gravity field variation.
     *  \param identifier Identifier of gravity field variation that is to be retrieved (empty by default; only required
     *  if multiple variations of same type are present)
     *  \return Object describing requested variation in the gravity field of this body.
     */
    std::pair<bool, std::shared_ptr<gravitation::GravityFieldVariations>>
    getGravityFieldVariation(
            const gravitation::BodyDeformationTypes &deformationType,
            const std::string identifier = "") {
        return gravityFieldVariationSet_->getGravityFieldVariation(deformationType, identifier);
    }

    //! Function to retrieve object containing all variations in the gravity field of this body.
    /*!
     * Function to retrieve object containing all variations in the gravity field of this body.
     * \return Object containing all variations in the gravity field of this body.
     */
    std::shared_ptr<gravitation::GravityFieldVariationsSet> getGravityFieldVariationSet() {
        return gravityFieldVariationSet_;
    }

    //! Function to retrieve container object with hardware systems present on/in body
    /*!
     * Function to retrieve container object with hardware systems present on/in body.
     * \return Container object with hardware systems present on/in body.
     */
    std::shared_ptr<system_models::VehicleSystems> getVehicleSystems()
    {
        if( vehicleSystems_ == nullptr )
        {
            vehicleSystems_ = std::make_shared< system_models::VehicleSystems >( );
        }
        return vehicleSystems_;
    }

    //! Function to set container object with hardware systems present on/in body
    /*!
     * Function to set container object with hardware systems present on/in body (typically only non-nullptr for a vehicle).
     * \param vehicleSystems Container object with hardware systems present on/in body.
     */
    void setVehicleSystems(const std::shared_ptr<system_models::VehicleSystems> vehicleSystems)
    {
        vehicleSystems_ = vehicleSystems;
    }

    std::shared_ptr< RigidBodyProperties > getMassProperties( )
    {
        return massProperties_;
    }

    void setMassProperties( const std::shared_ptr< RigidBodyProperties > massProperties )
    {
        if( gravityFieldModel_ != nullptr )
        {
            std::cerr<<"Warning, setting body mass distribution, but existing gravity field model and associated mass properties already found; overriding existing body mass properties"<<std::endl;
        }
        massProperties_ = massProperties;
    }

    //! Function to set the function returning body mass as a function of time
    /*!
     * Function to set the function returning body mass as a function of time
     * \param bodyMassFunction Function returning body mass as a function of time
     */
    void setBodyMassFunction(const std::function<double(const double)> bodyMassFunction)
    {
        if( massProperties_ == nullptr )
        {
            massProperties_ = std::make_shared< TimeDependentRigidBodyProperties >( bodyMassFunction );
        }
        else
        {
            if ( std::dynamic_pointer_cast<TimeDependentRigidBodyProperties>( massProperties_ ) != nullptr )
            {
                std::dynamic_pointer_cast<TimeDependentRigidBodyProperties>( massProperties_ )->setMassFunction(
                        bodyMassFunction );
            }
            else
            {
                throw std::runtime_error( "Error when resetting body mass function for " + bodyName_ +
                                          ", no compatible mass properties found" );
            }
        }
    }

    //! Function to set the body mass as being constant (i.e. time-independent)
    /*!
     * Function to set the body mass as being constant (i.e. time-independent)
     * \param bodyMass New constant body mass
     */
    void setConstantBodyMass(const double bodyMass)
    {
        if( massProperties_ == nullptr )
        {
            massProperties_ = std::make_shared< TimeDependentRigidBodyProperties >( bodyMass );
        }
        else
        {
            massProperties_->setCurrentMass( bodyMass );
        }
    }

    void setCurrentPropagatedBodyMass(const double bodyMass)
    {
        if( massProperties_ == nullptr )
        {
            massProperties_ = std::make_shared< TimeDependentRigidBodyProperties >( bodyMass );
        }
        else
        {
            massProperties_->setCurrentMass( bodyMass );
        }
    }

    //! Function to get the function returning body mass as a function of time
    /*!
     * Function to get the function returning body mass as a function of time
     * \return Function returning body mass as a function of time
     */
    std::function< double( const double ) > getBodyMassFunction( )
    {
        std::function<double(const double)> bodyMassFunction;
        if ( std::dynamic_pointer_cast<TimeDependentRigidBodyProperties>( massProperties_ ) != nullptr )
        {
            bodyMassFunction = std::dynamic_pointer_cast<TimeDependentRigidBodyProperties>( massProperties_ )->getMassFunction( );
        }
        else
        {
            throw std::runtime_error( "Error when getting body mass function for " + bodyName_ +
                                      ", no compatible mass properties found" );
        }
        return bodyMassFunction;
    }

    //! Function to update the body mass to the current time
    /*!
     * Function to update the body mass to the current time, using the bodyMassFunction_ function
     * \param time Current time
     */
    void updateMass(const double time)
    {
        if ( massProperties_ == nullptr)
        {
            throw std::runtime_error( "Error when updating body mass for " + bodyName_ +
                                      ", no mass properties found" );
        }
        else
        {
            massProperties_->updateMass( time );
        }
    }

    void updateMassDistribution(const double time)
    {
        if ( massProperties_ == nullptr)
        {
            throw std::runtime_error( "Error when updating body mass for " + bodyName_ +
                                     ", no mass properties found" );
        }
        else
        {
            massProperties_->updateMassDistribution( time );
        }
    }

    //! Function to retrieve the current body mass
    /*!
     * Function to retrieve the current body mass
     * \return Current body mass.
     */
    double getBodyMass( )
    {
        if( massProperties_ == nullptr )
        {
            throw std::runtime_error( "Error when retrieving mass of " + bodyName_ + ", no mass properties found" );
        }
        return massProperties_->getCurrentMass( );
    }

    Eigen::Vector3d getBodyFixedCenterOfMass( )
    {
        return massProperties_->getCurrentCenterOfMass( );
    }

    //! Function to retrieve the body moment-of-inertia tensor.
    /*!
     * Function to retrieve the body moment-of-inertia tensor.
     * \return  Body moment-of-inertia tensor.
     */
    Eigen::Matrix3d getBodyInertiaTensor( )
    {
        return massProperties_->getCurrentInertiaTensor( );
    }

    //! Function to (re)set the body moment-of-inertia tensor.
    /*!
     * Function to (re)set the body moment-of-inertia tensor.
     * \param bodyInertiaTensor Body moment-of-inertia tensor.
     */
    void setBodyInertiaTensor(const Eigen::Matrix3d &bodyInertiaTensor)
    {
        if( massProperties_ == nullptr )
        {
            throw std::runtime_error( "Error when setting inertia tensor; no body mass properties object found (add a body mass first)" );
        }
        else
        {
            if( std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ ) != nullptr )
            {
                std::dynamic_pointer_cast< TimeDependentRigidBodyProperties >( massProperties_ )->setInertiaTensor( bodyInertiaTensor );
            }
            else
            {
                throw std::runtime_error( "Error when trying to reset body inertia tensor, rigid body properties already exist, and are not of compatible type" );
            }
        }
    }

    //! Function to add a ground station to the body
    /*!
     * Function to add a ground station to the body
     * \param stationName Name of ground station
     * \param station Ground station object that is to be set
     */
    void addGroundStation(const std::string &stationName,
                          const std::shared_ptr<ground_stations::GroundStation> &station) {
        groundStationMap[stationName] = station;
    }

    //! Function to retrieve a ground station
    /*!
     * Function to retrieve a ground station
     * \param stationName Name of ground station
     * \return Ground station object that is retrieved
     */
    std::shared_ptr<ground_stations::GroundStation> getGroundStation(const std::string &stationName) const {
        if (groundStationMap.count(stationName) == 0)
        {
            std::vector< double > a;
            a.at( 0 );
            throw std::runtime_error("Error, station " + stationName + " does not exist");
        }

        return groundStationMap.at(stationName);
    }

    //! Function to retrieve full list of ground stations
    /*!
     * Function to retrieve full list of ground stations
     * \return Full list of ground stations
     */
    std::map<std::string, std::shared_ptr<ground_stations::GroundStation>> getGroundStationMap() const {
        return groundStationMap;
    }

    //! Function to recompute the internal variables of member variables that depend on the ephemerides bodies.
    /*!
     * Function to recompute the internal variables of member variables that depend on the ephemerides of this and other
     * bodies. This function is typically called after equations of motion have been computed and set in environment to
     * ensure full model consistency.
     */
    void updateConstantEphemerisDependentMemberQuantities() {
        if (std::dynamic_pointer_cast<gravitation::TimeDependentSphericalHarmonicsGravityField>(
                    gravityFieldModel_)
                != nullptr) {
            std::dynamic_pointer_cast<gravitation::TimeDependentSphericalHarmonicsGravityField>(
                        gravityFieldModel_)
                    ->updateCorrectionFunctions();
        }
    }

    //! Function to indicate that the state needs to be recomputed on next call to setStateFromEphemeris.
    /*!
     * Function to reset the time to which the state was last updated using setStateFromEphemeris function to nan, thereby
     * singalling that it needs to be recomputed upon next call.
     */
    void recomputeStateOnNextCall() {
        timeOfCurrentState_ = Time(TUDAT_NAN);
    }

    double getDoubleTimeOfCurrentState( )
    {
        return static_cast< double >( timeOfCurrentState_ );
    }


    //! Function to retrieve variable denoting whether this body is the global frame origin
    /*!
     * Function to retrieve variable denoting whether this body is the global frame origin
     * \return Variable denoting whether this body is the global frame origin
     */
    int getIsBodyGlobalFrameOrigin() {
        return bodyIsGlobalFrameOrigin_;
    }

    //! Function to set variable denoting whether this body is the global frame origin
    /*!
     * Function to set variable denoting whether this body is the global frame origin
     * \param bodyIsGlobalFrameOrigin Variable denoting whether this body is the global frame origin
     */
    void setIsBodyGlobalFrameOrigin(const int bodyIsGlobalFrameOrigin) {
        bodyIsGlobalFrameOrigin_ = bodyIsGlobalFrameOrigin;
    }

    //! Function to define whether the body is currently being propagated, or not
    /*!
     *  Function to define whether the body is currently being propagated, or not
     *  \param isBodyInPropagation Boolean defining whether the body is currently being propagated, or not
     */
    void setIsBodyInPropagation(const bool isBodyInPropagation);

//    void setSuppressDependentOrientationCalculatorWarning(const bool suppressDependentOrientationCalculatorWarning) {
//        suppressDependentOrientationCalculatorWarning_ = suppressDependentOrientationCalculatorWarning;
//    }

    std::string getBodyName( ){ return bodyName_; }

    void setBodyName( const std::string bodyName ){ bodyName_ = bodyName; }

protected:
private:
    //! Variable denoting whether this body is the global frame origin (1 if true, 0 if false, -1 if not yet set)
    int bodyIsGlobalFrameOrigin_;

    //! Current state.
    Eigen::Vector6d currentState_;

    //! Current state with long double precision.
    Eigen::Matrix<long double, 6, 1> currentLongState_;

    //! Current state.
    Eigen::Vector6d currentBarycentricState_;

    //! Current state with long double precision.
    Eigen::Matrix<long double, 6, 1> currentBarycentricLongState_;

    //! Time at which state was last set from ephemeris
    Time timeOfCurrentState_;

    //! Class returning the state of this body's ephemeris origin w.r.t. the global origin (as typically created by
    //! setGlobalFrameBodyEphemerides function).
    std::shared_ptr<BaseStateInterface> ephemerisFrameToBaseFrame_;

    Eigen::Quaterniond currentRotationToLocalFrame_;

    Eigen::Quaterniond currentRotationToGlobalFrame_;
    //! Current first derivative w.r.t. time of the rotation matrix from the global to the
    //! body-fixed frame.
    Eigen::Matrix3d currentRotationToLocalFrameDerivative_;

    Eigen::Vector3d currentAngularVelocityVectorInGlobalFrame_;

    //! Current angular velocity vector for body's rotation, expressed in the body-fixed frame.
    Eigen::Vector3d currentAngularVelocityVectorInLocalFrame_;

//    //! Mass of body (default set to zero, calculated from GravityFieldModel when it is set).
//    double currentMass_;

//    //! Function returning body mass as a function of time.
//    std::function<double(const double)> bodyMassFunction_;

//    //! Body moment-of-inertia tensor.
//    Eigen::Matrix3d bodyInertiaTensor_;

//    //! Body scaled mean moment of inertia
//    double scaledMeanMomentOfInertia_;

    //! Ephemeris of body.
    std::shared_ptr<ephemerides::Ephemeris> bodyEphemeris_;

    //! Gravity field model of body.
    std::shared_ptr<gravitation::GravityFieldModel> gravityFieldModel_;

    //! Object containing all variations in the gravity field of this body.
    std::shared_ptr<gravitation::GravityFieldVariationsSet> gravityFieldVariationSet_;

    //! Atmosphere model of body.
    std::shared_ptr<aerodynamics::AtmosphereModel> atmosphereModel_;

    //! Shape model of body.
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> shapeModel_;

    //! Aerodynamic coefficient model of body.
    std::shared_ptr<aerodynamics::AerodynamicCoefficientInterface> aerodynamicCoefficientInterface_;

    //! Object used for calculating current aerodynamic angles, altitude, etc.
    std::shared_ptr<aerodynamics::FlightConditions> aerodynamicFlightConditions_;

    //! Rotation model of body.
    std::shared_ptr<ephemerides::RotationalEphemeris> rotationalEphemeris_;

    std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > > bodyDeformationModels_;

    std::shared_ptr< RigidBodyProperties > massProperties_;

    //! List of radiation pressure models for the body, with the sources bodies as key
    // RP-OLD
    std::map<std::string, std::shared_ptr<electromagnetism::RadiationPressureInterface>>
            radiationPressureInterfaces_;

    //! Radiation source model of the body.
    std::shared_ptr<electromagnetism::RadiationSourceModel> radiationSourceModel_;

    //! Radiation pressure target model of the body.
    std::shared_ptr<electromagnetism::RadiationPressureTargetModel> radiationPressureTargetModel_;

    //! List of ground station objects on Body
    std::map<std::string, std::shared_ptr<ground_stations::GroundStation>> groundStationMap;

    //! Container object with hardware systems present on/in body (typically only non-nullptr for a vehicle).
    std::shared_ptr<system_models::VehicleSystems> vehicleSystems_;

    //!  Boolean defining whether the body is currently being propagated, or not
    bool isBodyInPropagation_ = false;

//    bool suppressDependentOrientationCalculatorWarning_ = false;

    std::string bodyName_;

    bool isStateSet_;

    bool isRotationSet_;
};


//! Typdef for a list of body objects (as unordered_map for efficiency reasons)
//typedef std::unordered_map< std::string, std::shared_ptr< Body > > SystemOfBodies;

std::shared_ptr< ephemerides::ReferenceFrameManager > createFrameManager(
        const std::unordered_map< std::string, std::shared_ptr< Body > > bodies );

//! Function to define the global origin and orientation of the reference frame
/*!
 * Function to define the global origin and orientation of the reference frame that is to be used in
 * the simulations.  This function checks the origin and orientation of the Ephemeris and
 * RotationalEphemeris, and checks whether their origin/orientation is the same as that
 * globalFrameOrigin and globalFrameOrientation provided as input.  In particular, this function
 * sets the ephemerisFrameToBaseFrameFunction_ anf ephemerisFrameToBaseFrameLongFunction_ variables
 * of the Body objects, which provide a time-dependent translation of the global origin to the
 * body's ephemeris origin. In case of an inconsistency in the current and requried frames, this
 * function throws an error.
 * \param bodies List of body objects that constitute the environment.
 * \param globalFrameOrigin Global reference frame origin.
 * \param globalFrameOrientation Global referencef frame orientation.
 */
template< typename StateScalarType = double, typename TimeType = double >
void setGlobalFrameBodyEphemerides( const std::unordered_map< std::string, std::shared_ptr< Body > > bodies,
                                    const std::string& globalFrameOrigin,
                                    const std::string& globalFrameOrientation )
{
    using namespace tudat::simulation_setup;
    std::string ephemerisFrameOrigin;
    std::string ephemerisFrameOrientation;
    std::string rotationModelFrame;

    std::vector< std::string > globalFrameOriginChain;

    // Get chain of ephemeris frame origins of global frame origin (if it is not SSB
    if( globalFrameOrigin != "SSB" )
    {
        if( bodies.count( globalFrameOrigin ) == 0 )
        {
            throw std::runtime_error(
                        "Error, body non-barycentric global frame origin selected, but this body " + globalFrameOrigin +
                        " is not found." );
        }
        else
        {
            std::string currentOrigin = globalFrameOrigin;
            while( currentOrigin != "SSB" )
            {
                std::shared_ptr< ephemerides::Ephemeris > currentEphemeris =
                        bodies.at( currentOrigin )->getEphemeris( );
                if( currentEphemeris == nullptr )
                {
                    throw std::runtime_error(
                                "Error, body non-barycentric global frame origin selected, but body " + currentOrigin +
                                " in chain has no ephemeris." );
                }
                else
                {
                    ephemerisFrameOrientation = currentEphemeris->getReferenceFrameOrientation( );
                    if( ephemerisFrameOrientation != globalFrameOrientation )
                    {
                        throw std::runtime_error(
                                    "Error, ephemeris orientation of body " + currentOrigin
                                    + " is not the same as global orientation " + ephemerisFrameOrientation
                                    + ", " + globalFrameOrientation );
                    }
                    currentOrigin = currentEphemeris->getReferenceFrameOrigin( );
                }

                if( std::find( globalFrameOriginChain.begin( ), globalFrameOriginChain.end( ), currentOrigin ) !=
                        globalFrameOriginChain.end( ) )
                {
                    throw std::runtime_error(
                                "Error, body non-barycentric global frame origin selected, but body " + currentOrigin +
                                " already found in origin chain." );
                }
                else
                {
                    globalFrameOriginChain.push_back( currentOrigin );
                }
            }
        }
    }

    // Iterate over all bodies
    for( auto bodyIterator : bodies )
    {
        // Check id body contains an ephemeris
        if( bodyIterator.second->getEphemeris( ) != nullptr )
        {
            // Retrieve ephemeris origin
            ephemerisFrameOrigin = bodyIterator.second->getEphemeris( )->getReferenceFrameOrigin( );

            // Check if ephemeris origin differs from global origin.
            if( ephemerisFrameOrigin != globalFrameOrigin )
            {
                // Make correction to SSB if it is global frame origin
                if( globalFrameOrigin == "SSB" )
                {
                    // Check if correction can be made
                    if( bodies.count( ephemerisFrameOrigin ) == 0 )
                    {
                        throw std::runtime_error(
                                    "Error, body " + bodyIterator.first + " has ephemeris in frame " +
                                    ephemerisFrameOrigin + ", but no conversion to frame " + globalFrameOrigin +
                                    " can be made" );
                    }
                    else
                    {
                        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                std::bind( &Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                                           bodies.at( ephemerisFrameOrigin ), std::placeholders::_1 );
                        std::shared_ptr< BaseStateInterface > baseStateInterface =
                                std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                    ephemerisFrameOrigin, stateFunction );
                        bodyIterator.second->setEphemerisFrameToBaseFrame( baseStateInterface );
                    }
                }
                // Make correction to global frame origin (if not SSB)
                else
                {
                    // Set barycentric state function of global frame origin
                    if( globalFrameOrigin == bodyIterator.first )
                    {
                        std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager =
                                createFrameManager( bodies );
                        frameManager->getEphemeris( globalFrameOrigin, "SSB" );

                        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                std::bind( &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< StateScalarType, TimeType >,
                                           frameManager->getEphemeris( globalFrameOrigin, "SSB" ), std::placeholders::_1 );

                        std::shared_ptr< BaseStateInterface > baseStateInterface =
                                std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                    globalFrameOrigin, stateFunction, true );
                        bodyIterator.second->setEphemerisFrameToBaseFrame( baseStateInterface );

                    }
                    // Set correction function if ephemeris origin is SSB
                    else if( ephemerisFrameOrigin == "SSB" )
                    {

                        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                std::bind( &Body::getGlobalFrameOriginBarycentricStateFromEphemeris< StateScalarType, TimeType >,
                                           bodies.at( globalFrameOrigin ), std::placeholders::_1 );
                        std::shared_ptr< BaseStateInterface > baseStateInterface =
                                std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                    globalFrameOrigin, stateFunction, true );
                        bodyIterator.second->setEphemerisFrameToBaseFrame( baseStateInterface );
                    }
                    else
                    {
                        // Check if correction can be made
                        if( bodies.count( ephemerisFrameOrigin ) == 0 )
                        {
                            throw std::runtime_error(
                                        "Error, body " + bodyIterator.first + " has ephemeris in frame " +
                                        ephemerisFrameOrigin + ", but no conversion to frame " + globalFrameOrigin +
                                        " can be made" );
                        }
                        else
                        {
                            // Set correction function from ephemeris origin to global frame origin
                            std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction =
                                    std::bind( &Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                                               bodies.at( ephemerisFrameOrigin ), std::placeholders::_1 );
                            std::shared_ptr< BaseStateInterface > baseStateInterface =
                                    std::make_shared< BaseStateInterfaceImplementation< TimeType, StateScalarType > >(
                                        ephemerisFrameOrigin, stateFunction, false );
                            bodyIterator.second->setEphemerisFrameToBaseFrame( baseStateInterface );
                        }
                    }
                }
            }

            // Retrieve ephemeris orientation
            ephemerisFrameOrientation = bodyIterator.second->getEphemeris( )->getReferenceFrameOrientation( );
            // If two are not equal, throw error.
            if( ephemerisFrameOrientation != globalFrameOrientation )
            {
                throw std::runtime_error(
                            "Error, ephemeris orientation of body " + bodyIterator.first
                            + " is not the same as global orientation " + ephemerisFrameOrientation
                            + ", " + globalFrameOrientation );
            }


        }

        // Set global frame origin identifiers
        if( globalFrameOrigin == bodyIterator.first )
        {
            bodyIterator.second->setIsBodyGlobalFrameOrigin( 1 );
        }
        else
        {
            bodyIterator.second->setIsBodyGlobalFrameOrigin( 0 );
        }

        // Check if body has rotational ephemeris.
        if( bodyIterator.second->getRotationalEphemeris( ) != nullptr )
        {
            // Check if rotational ephemeris base frame orienatation is equal to to global orientation.
            rotationModelFrame = bodyIterator.second->getRotationalEphemeris( )->getBaseFrameOrientation( );

            // Throw error if two frames are not equal.
            if( rotationModelFrame != globalFrameOrientation )
            {
                throw std::runtime_error(
                            "Error, rotation base orientation of body " + bodyIterator.first +
                            " is not the same as global orientation " + rotationModelFrame + ", " +
                            globalFrameOrientation );
            }
        }
    }

    // Set body state-dependent environment variables
    for( auto bodyIterator : bodies  )
    {
        bodyIterator.second->updateConstantEphemerisDependentMemberQuantities( );
    }

}


template< typename StateScalarType = double, typename TimeType = double >
void addEmptyEphemeris(
        const std::shared_ptr< Body > body,
        const std::string centralBody, const std::string frameOrientation,
        const bool ephemerisIsMultiArc = false, const bool overrideExisting = false )
{
    if( body->getEphemeris( ) != nullptr && !overrideExisting )
    {
        throw std::runtime_error( "Errror when adding empty default ephemeris; body already posseses ephemeris!" );
    }

    if( !ephemerisIsMultiArc )
    {
        body->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                                nullptr, centralBody, frameOrientation ) );
    }
    else
    {
        body->setEphemeris( std::make_shared< ephemerides::MultiArcEphemeris >(
                                nullptr, centralBody, frameOrientation ) );
    }
}


class SystemOfBodies
{
public:
    SystemOfBodies( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000",
                    const std::unordered_map< std::string, std::shared_ptr< Body > >& bodyMap =
            std::unordered_map< std::string, std::shared_ptr< Body > >( ) ):
        frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ), bodyMap_( bodyMap ){ }

    std::shared_ptr< Body > at( const std::string& bodyName ) const
    {
        if( bodyMap_.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving body " + bodyName + " from SystemOfBodies, no such body exists" );
        }
        return bodyMap_.at( bodyName );
    }

    std::shared_ptr< Body > getBody( const std::string& bodyName ) const
    {
        return at( bodyName );
    }

    int count( const std::string& bodyName ) const
    {
        return static_cast< int >( bodyMap_.count( bodyName ) );
    }

    bool doesBodyExist( const std::string& bodyName ) const
    {
        return ( this->count( bodyName ) != 0 );
    }

    std::vector< std::string > getListOfBodies( ) const
    {
        return utilities::createVectorFromUnorderedMapKeys( bodyMap_ );
    }

    int getNumberOfBodies( ) const
    {
        return static_cast< int >( bodyMap_.size( ) );
    }

    template< typename StateScalarType = double , typename TimeType = double >
    void createEmptyBody( const std::string bodyName, const bool processBody = true )
    {
        bodyMap_[ bodyName ] = std::make_shared< Body >( );
        bodyMap_[ bodyName ]->setBodyName( bodyName );
        if( processBody )
        {
            processBodyFrameDefinitions< StateScalarType, TimeType >( );
        }
    }

    template< typename StateScalarType = double , typename TimeType = double >
    void addBody( std::shared_ptr< Body > bodyToAdd, const std::string bodyName, const bool processBody = true )
    {
        bodyMap_[ bodyName ] = bodyToAdd;
        bodyMap_[ bodyName ]->setBodyName( bodyName );
        if( processBody )
        {
            processBodyFrameDefinitions< StateScalarType, TimeType >( );
        }
    }

    const std::unordered_map< std::string, std::shared_ptr< Body > >& getMap( ) const { return bodyMap_; }

    template< typename StateScalarType = double , typename TimeType = double >
    void processBodyFrameDefinitions( ) const
    {
        setGlobalFrameBodyEphemerides< StateScalarType, TimeType >( bodyMap_, frameOrigin_, frameOrientation_);

//        for( auto bodyIterator : bodyMap_ )
//        {
//            bodyIterator.second->setBaseFrameFunction(
//                        std::bind( &SystemOfBodies::processBodyFrameDefinitions, this ) );
//        }
    }


    std::string getFrameOrigin( ) const
    {
        return frameOrigin_;
    }

    std::string getFrameOrientation( ) const
    {
        return frameOrientation_;
    }

    std::unordered_map< std::string, std::shared_ptr< Body > > getMap( )
    {
        return bodyMap_;
    }

    void deleteBody( const std::string bodyName )
    {
        bodyMap_.at( bodyName ).reset( );
        bodyMap_.erase( bodyName );

    }
private:

    std::string frameOrigin_;

    std::string frameOrientation_;

    std::unordered_map< std::string, std::shared_ptr< Body > > bodyMap_;

};

double getBodyGravitationalParameter( const SystemOfBodies& bodies, const std::string bodyName );

//! Function ot retrieve the common global translational state origin of the environment
/*!
 * Function ot retrieve the common global translational state origin of the environment. This function throws an exception
 * if multiple bodies are found as the frame origin
 * \param bodies List of body objects.
 * \return Global translational state origin of the environment
 */
std::string getGlobalFrameOrigin(const SystemOfBodies &bodies);

//! Function to set whether the bodies are currently being propagated, or not
/*!
 * Function to set whether the bodies are currently being propagated, or not
 * \param bodies List of body objects.
 * \param areBodiesInPropagation Boolean defining whether the bodies are currently being propagated, or not
 */
void setAreBodiesInPropagation(const SystemOfBodies &bodies,
                               const bool areBodiesInPropagation);

bool isReferencePointGroundStation( const std::shared_ptr< Body > body,
                                    const std::string& referencePointName );

bool isReferencePointGroundStation( const SystemOfBodies &bodies,
                                    const std::string& bodyName,
                                    const std::string& referencePointName );

//! Function to compute the acceleration of a body, using its ephemeris and finite differences
/*!
 *  Function to compute the acceleration of a body, using its ephemeris and 8th order finite difference and 100 s time step
 *  \param bodyWithAcceleration Body for which acceleration is to be computed
 *  \param nominalEvalutationTime Time at which acceleration is to be evaluated.
 */
template<typename StateScalarType = double, typename TimeType = double>
Eigen::Matrix<StateScalarType, 3, 1> getBodyAccelerationInBaseFramefromNumericalDifferentiation(
        const std::shared_ptr<Body> bodyWithAcceleration,
        const TimeType nominalEvalutationTime) {
    std::function<Eigen::Matrix<StateScalarType, 6, 1>(const TimeType)> bodyStateFunction =
            std::bind(&Body::getStateInBaseFrameFromEphemeris<StateScalarType, TimeType>, bodyWithAcceleration, std::placeholders::_1);
    return numerical_derivatives::computeCentralDifferenceFromFunction(
                bodyStateFunction, nominalEvalutationTime, 100.0, numerical_derivatives::order8)
            .segment(3, 3);
}

}// namespace simulation_setup

}// namespace tudat

#endif// TUDAT_BODY_H
