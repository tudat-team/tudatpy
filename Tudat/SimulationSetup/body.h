/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h>
#include <Tudat/Astrodynamics/Ephemerides/ephemeris.h>
#include <Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h>
#include <Tudat/Astrodynamics/Gravitation/gravityFieldModel.h>
#include <Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h>
#include <Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h>
#include <Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h>

namespace tudat
{

namespace simulation_setup
{

//! Body class representing the properties of a celestial body (natural or artificial).
/*!
 *  Body class representing the properties of a celestial body (natural or artificial). By storing
 *  all properties of bodies (ephemeris, rotation, gravity, etc.) in a set of body objects,
 *  the simulation environment can be defined in a clear and modular way. To create body
 *  objects, the createBodies.h function provides a range of functionality. The
 *  createAccelerationModels.h file provides functions to use body objects to create acceleration
 *  objects.
 */
class Body
{
public:

    //! Constructor for a body
    /*!
     * Constructor for a body, sets current state (with zero default value).
     * \param state Current state of body at initialization (default = zeroes).
     */
    Body( const basic_mathematics::Vector6d& state =
            basic_mathematics::Vector6d::Zero( ) )
        : currentState_( state ),
          ephemerisFrameToBaseFrameFunction_( boost::lambda::constant( basic_mathematics::Vector6d::Zero( ) ) ),
          ephemerisFrameToBaseFrameLongFunction_( boost::lambda::constant( Eigen::Matrix< long double, 6, 1 >::Zero( ) ) ),
          currentRotationToLocalFrame_( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
          currentRotationToLocalFrameDerivative_( Eigen::Matrix3d::Zero( ) ),
          currentAngularVelocityVectorInGlobalFrame_( Eigen::Vector3d::Zero( ) ),
          bodyMassFunction_( NULL )
    {
        currentLongState_ = currentState_.cast< long double >( );
    }



    //! Set current state of body manually
    /*!
     * Set current state of body manually, which must be in the global frame. Note that this
     * function does not set the currentLongState_, use the setLongState when needing the use of the
     * long precision current state.
     * \param state Current state of the body that is set.
     */
    void setState( const basic_mathematics::Vector6d& state )
    {
        currentState_ = state;
    }

    //! Set current state of body manually in long double precision.
    /*!
     * Set current state of body manually in long double precision. State must be in the global
     * frame.  Note that this function sets both the currentState_ and currentLongState_ variables
     * (currentLongState_ directly and currentState_ by casting the input to double entries).
     * \param longState Current state of the body that is set, in long double precision.
     */
    void setLongState( const Eigen::Matrix< long double, 6, 1 >& longState )
    {
        currentLongState_ = longState;
        currentState_ = longState.cast< double >( );
    }

    //! Function to get the state of the current body from the ephemeris in the global frame
    /*!
     * Function to get the state of the current body from the ephemeris in the global frame. Calling
     * this function calls the bodyEphemeris_ object, not the currentState_ variable. In addition,
     * it adds the global frame origin w.r.t. the ephemeris origin using the
     * ephemerisFrameToBaseFrameFunction_ variable.
     * \param time Time at which the global state is to be retrieved.
     * \return Global state at current time, obtained from ephemeris and ephemerisFrameToBaseFrameFunction_
     */
    basic_mathematics::Vector6d getStateInBaseFrameFromEphemeris( const double& time )
    {
        return bodyEphemeris_->getCartesianStateFromEphemeris( time ) + ephemerisFrameToBaseFrameFunction_( time );
    }

    //! Function to get the long precisien state of the current body from ephemeris in the global frame
    /*!
     * Function to get the state of the current body, in long double precision, from the ephemeris
     * in the global frame.  Calling this function calls the bodyEphemeris_ object, not the
     * currentLongState_ variable.  In addition, it adds the global frame origin w.r.t. the
     * ephemeris origin using the ephemerisFrameToBaseFrameLongFunction_ variable.
     * \param time Time at which the global state is to be retrieved.
     * \return Global state in long double precision at current time, obtained from ephemeris and
     * ephemerisFrameToBaseFrameLongFunction_
     */
    Eigen::Matrix< long double, 6, 1 > getLongStateInBaseFrameFromEphemeris( const double& time )
    {
        return bodyEphemeris_->getCartesianLongStateFromEphemeris( time ) + ephemerisFrameToBaseFrameLongFunction_( time );
    }



    //! Function to set the state of the current body from the ephemeris in the global frame
    /*!
     * Function to set the state of the current body from the ephemeris in the global frame
     * (currentState_ variable).  Calling this function uses the bodyEphemeris_ object, and adds the
     * global frame origin w.r.t. the ephemeris origin using the ephemerisFrameToBaseFrameFunction_
     * variable.Note that this function does not set the currentLongState_, use the setLongState
     * when needing the use of the long precision current state.
     * \param time Time at which the global state is to be set.
     */
    void setStateFromEphemeris( const double& time )
    {
        currentState_ = bodyEphemeris_->getCartesianStateFromEphemeris( time ) + ephemerisFrameToBaseFrameFunction_( time );
    }

    //! Function to set long precision state of the current body from ephemeris in the global frame.
    /*!
     * Function to set the state of the current body from the ephemeris in the global frame in long
     * double precision (currentLongState_ variable). Calling this function uses the bodyEphemeris_
     * object, and adds the global frame origin w.r.t. the ephemeris origin using the
     * ephemerisFrameToBaseFrameLongFunction_ variable.  Note that this function sets both the
     * currentState_ and currentLongState_ variables (currentLongState_ directly and currentState_
     * by casting the input to double entries).
     * \param time Time at which the global state is to be set.
     */
    void setLongStateFromEphemeris( const double& time )
    {
        currentLongState_ = bodyEphemeris_->getCartesianLongStateFromEphemeris( time ) +
                ephemerisFrameToBaseFrameLongFunction_( time );
        currentState_ = currentLongState_.cast< double >( );
    }



    //! Templated function to set the state manually.
    /*!
     * Templated function to set the state manually, calls either setState or setLongState function.
     * \param state Current state of the body that is set, with StateScalarType precision.
     */
    template< typename StateScalarType >
    void setTemplatedState( const Eigen::Matrix< StateScalarType, 6, 1 >& state );

    //! Templated function to set the current state of the body from its ephemeris and
    //! global-to-ephemeris-frame function.
    /*!
     * Templated function to set the current state of the body from its ephemeris and
     * global-to-ephemeris-frame function. It calls either the setStateFromEphemeris or
     * setLongStateFromEphemeris function.
     * \param time Time at which the global state is to be set.
     */
    template< typename StateScalarType, typename TimeType >
    void setTemplatedStateFromEphemeris( const TimeType& time );

    //! Templated function to get the current state of the body from its ephemeris and
    //! global-to-ephemeris-frame function.
    /*!
     * Templated function to sgt the current state of the body from its ephemeris and
     * global-to-ephemeris-frame function.  It calls either the getStateInBaseFrameFromEphemeris or
     * getLongStateInBaseFrameFromEphemeris function.
     * \param time Time at which the global state is to be set.
     */
    template< typename StateScalarType, typename TimeType >
    Eigen::Matrix< StateScalarType, 6, 1 > getTemplatedStateInBaseFrameFromEphemeris( const TimeType& time );



    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    basic_mathematics::Vector6d getState( ) { return currentState_; }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getPosition( ) { return currentState_.segment( 0, 3 ); }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getVelocity( ) { return currentState_.segment( 3, 3 ); }


    //! Get current state, in long double precision
    /*!
     * Returns the internally stored current state vector, in long double precision
     * \return Current state, in long double precisio
     */
    Eigen::Matrix< long double, 6, 1 > getLongState( ) { return currentLongState_; }

    //! Get current position, in long double precision
    /*!
     * Returns the internally stored current position vector, in long double precision
     * \return Current position, in long double precision
     */
    Eigen::Matrix< long double, 3, 1 > getLongPosition( ) { return currentLongState_.segment( 0, 3 ); }

    //! Get current velocity, in long double precision.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity, in long double precision
     */
    Eigen::Matrix< long double, 3, 1 > getLongVelocity( ) { return currentLongState_.segment( 3, 3 ); }


    //! Function to set the rotation from global to body-fixed frame at given time
    /*!
     * Function to set the rotation from global to body-fixed frame at given time, using the
     * rotationalEphemeris_ member object
     * \param time Time at which the rotation is to be retrieved.
     */
    void setCurrentRotationToLocalFrameFromEphemeris( const double time )
    {
        if( rotationalEphemeris_!= NULL )
        {
            currentRotationToLocalFrame_ = rotationalEphemeris_->getRotationToTargetFrame( time );
        }
        else
        {
            throw std::runtime_error(
                        "Error, no rotationalEphemeris_ found in Body::setCurrentRotationToLocalFrameFromEphemeris" );
        }
    }

    //! Function to set the rotation matrix derivative from global to body-fixed frame at given time
    /*!
     * Function to set the rotation matrix derivative from global to body-fixed frame at given time,
     * using the rotationalEphemeris_ member object
     * \param time Time at which the rotation matrix derivative is to be retrieved.
     */
    void setCurrentRotationToLocalFrameDerivativeFromEphemeris( const double time )
    {
        if( rotationalEphemeris_!= NULL )
        {
            currentRotationToLocalFrameDerivative_
                    = rotationalEphemeris_->getDerivativeOfRotationToTargetFrame( time );
        }
        else
        {
            throw std::runtime_error(
                   "Error, no rotationalEphemeris_ found in Body::setCurrentRotationToLocalFrameDerivativeFromEphemeris" );
        }
    }

    //! Function to set the angular velocity vector in the global frame at given time
    /*!
     * Function to set the angular velocity vector in the global frame at given time, using the
     * rotationalEphemeris_ member object
     * \param time Time at which the angular velocity vector in the global frame is to be retrieved.
     */
    void setCurrentAngularVelocityVectorInGlobalFrame( const double time )
    {
        if( rotationalEphemeris_!= NULL )
        {
            currentAngularVelocityVectorInGlobalFrame_
                    = rotationalEphemeris_->getRotationalVelocityVectorInBaseFrame( time );
        }
        else
        {
            throw std::runtime_error(
                        "Error, no rotationalEphemeris_ found in Body::setCurrentAngularVelocityVectorInGlobalFrame" );
        }
    }

    //! Function to set the full rotational state at given time
    /*!
     * Function to set the full rotational state at (rotation from global to body-fixed frame
     * rotation matrix derivative from global to body-fixed frame and angular velocity vector in the
     * global frame) at given time, using the rotationalEphemeris_ member object.
     * \param time Time at which the angular velocity vector in the global frame is to be retrieved.
     */
    void setCurrentRotationalStateToLocalFrameFromEphemeris( const double time )
    {
        if( rotationalEphemeris_!= NULL )
        {
            rotationalEphemeris_->getFullRotationalQuantitiesToTargetFrame(
                        currentRotationToLocalFrame_, currentRotationToLocalFrameDerivative_,
                        currentAngularVelocityVectorInGlobalFrame_, time );
        }
        else
        {
            throw std::runtime_error(
                        "Error, no rotationalEphemeris_ found in Body::setCurrentRotationalStateToLocalFrameFromEphemeris" );
        }
    }


    //! Get current rotation from body-fixed to inertial frame.
    /*!
     *  Get current rotation from body-fixed to inertial frame, as set from the rotationalEphemeris_
     *  by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameFromEphemeris function.  If body has no rotational ephemeris,
     *  an identity quaternion (no rotation) is returned.
     *  \return Current rotation from body-fixed to inertial frame
     */
    Eigen::Quaterniond getCurrentRotationToGlobalFrame( )
    {
        return currentRotationToLocalFrame_.inverse( );
    }

    //! Get current rotation from inertial to body-fixed frame.
    /*!
     *  Get current rotation from inertial to body-fixed frame, as set from the rotationalEphemeris_
     *  by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameFromEphemeris function.  If body has no rotational ephemeris,
     *  an identity quaternion (no rotation) is returned.
     *  \return Current rotation from inertial to body-fixed frame
     */
    Eigen::Quaterniond getCurrentRotationToLocalFrame( )
    {
        return currentRotationToLocalFrame_;
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
        return currentRotationToLocalFrameDerivative_.transpose( );
    }

    //! Get current rotation matrix derivative from global to body-fixed frame.
    /*!
     *  Get current rotation matrix derivative from global to body-fixed frame, as set from the
     *  rotationalEphemeris_ by the setCurrentRotationalStateToLocalFrameFromEphemeris or
     *  setCurrentRotationToLocalFrameDerivativeFromEphemeris function. If body has no rotational
     *  ephemeris, an zero matrix (no rotation) is returned.
     *  \return Current otation matrix derivative from global to body-fixed frame.
     */
    Eigen::Matrix3d getCurrentRotationMatrixDerivativeToLocalFrame( )
    {
        return currentRotationToLocalFrameDerivative_;
    }

    //! Function to set the ephemeris of the body.
    /*!
     *  Function to set the ephemeris of the body, which is used to represent the (a priori)
     *  state history of the body.
     *  \param bodyEphemeris New ephemeris of the body.
     */
    void setEphemeris( const boost::shared_ptr< ephemerides::Ephemeris > bodyEphemeris )
    {
        bodyEphemeris_ = bodyEphemeris;
    }

    //! Function to set the function returning the state of this body's ephemeris origin
    /*!
     *  Function to set the function returning the state of this body's ephemeris origin w.r.t. the
     *  global origin (typically done by setGlobalFrameBodyEphemerides function).
     *  \param ephemerisFrameToBaseFrameFunction Function providing state of ephemeris origin w.r.t
     *  global origin.
     */
    void setBaseFrameFunction( const boost::function< basic_mathematics::Vector6d( const double& ) >
                               ephemerisFrameToBaseFrameFunction )
    {
        ephemerisFrameToBaseFrameFunction_ = ephemerisFrameToBaseFrameFunction;
    }

    //! Function to set the function returning the state (in long double) of body's ephemeris origin
    /*!
     *  Function to set the function returning the state of this body's ephemeris origin w.r.t. the
     *  global origin, with long double precision (typically done by setGlobalFrameBodyEphemerides
     *  function).
     *  \param ephemerisFrameToBaseFrameLongFunction Function providing state of ephemeris origin
     *  w.r.t global origin, with long double precision.
     */
    void setBaseFrameLongFunction( const boost::function< Eigen::Matrix< long double, 6, 1 >( const double& ) >
                                   ephemerisFrameToBaseFrameLongFunction )
    {
        ephemerisFrameToBaseFrameLongFunction_ = ephemerisFrameToBaseFrameLongFunction;
    }

    //! Function to set the gravity field of the body.
    /*!
     *  Function to set the gravity field of the body; input is also used to (re)set the mass of the
     *  body.
     *  \param gravityFieldModel New gravity field of the body.
     */
    void setGravityFieldModel(
            const boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel )
    {
        gravityFieldModel_ = gravityFieldModel;

        // Update current mass of body, provide warning
        if( bodyMassFunction_ != NULL )
        {
            std::cerr<<"Warning when settings gravity field model for body, mass function already found: resetting"<<std::endl;
        }
        currentMass_ = gravityFieldModel_->getGravitationalParameter( )
                       / physical_constants::GRAVITATIONAL_CONSTANT;
        bodyMassFunction_ = boost::lambda::constant( currentMass_ );
    }

    //! Function to set the atmosphere model of the body.
    /*!
     *  Function to set the atmosphere model of the body.
     *  \param atmosphereModel Atmosphere model of the body.
     */
    void setAtmosphereModel(
            const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel )
    {
        atmosphereModel_ = atmosphereModel;
    }

    //! Function to set the rotation model of the body.
    /*!
     *  Function to set the rotation model of the body.
     *  \param rotationalEphemeris Rotation model of the body.
     */
    void setRotationalEphemeris(
            const boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris )
    {
        rotationalEphemeris_ = rotationalEphemeris;
    }

    //! Function to set the shape model of the body.
    /*!
     *  Function to set the shape model of the body.
     *  \param shapeModel Shape model of the body.
     */
    void setShapeModel( const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel )
    {
        shapeModel_ = shapeModel;
    }

    //! Function to set the aerodynamic coefficient interface of the body.
    /*!
     *  Function to set the aerodynamic coefficient interface of the body.
     *  \param aerodynamicCoefficientInterface Aerodynamic coefficient interface of the body.
     */
    void setAerodynamicCoefficientInterface(
            const boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
            aerodynamicCoefficientInterface)
    {
        aerodynamicCoefficientInterface_ = aerodynamicCoefficientInterface;
    }

    //! Function to set the body flight conditions
    /*!
     * Function to set the body flight conditions, which calculates current aerodynamic angles,
     * altitude, etc.
     * \param aerodynamicFlightConditions Body flight conditions
     */
    void setFlightConditions(
            const boost::shared_ptr< aerodynamics::FlightConditions > aerodynamicFlightConditions )
    {
        aerodynamicFlightConditions_ = aerodynamicFlightConditions;
    }

    //! Function to set the radiation pressure interface of the body, for a single radiation source.
    /*!
     *  Function to set the radiation pressure interface of the body, for a single radiation source
     *  \param radiatingBody Name of body that is the source of the radiation.
     *  \param radiationPressureInterface Radiation pressure interface of the body.
     */
    void setRadiationPressureInterface(
            const std::string& radiatingBody,
            const boost::shared_ptr< electro_magnetism::RadiationPressureInterface >
                radiationPressureInterface )
    {
        radiationPressureInterfaces_[ radiatingBody ] = radiationPressureInterface;
    }

    void setGravityFieldVariationSet(
            const boost::shared_ptr< gravitation::GravityFieldVariationsSet >
                gravityFieldVariationSet )
    {
        gravityFieldVariationSet_ = gravityFieldVariationSet;
    }

    //! Function to get the gravity field model of the body.
    /*!
     *  Function to get the gravity field model of the body.
     *  \return Gravity field model of the body.
     */
    boost::shared_ptr< gravitation::GravityFieldModel > getGravityFieldModel( )
    {
        return gravityFieldModel_;
    }

    //! Function to get the ephemeris of the body.
    /*!
     *  Function to get the ephemeris of the body.
     *  \return Ephemeris of the body.
     */
    boost::shared_ptr< ephemerides::Ephemeris > getEphemeris( )
    {
        return bodyEphemeris_;
    }

    //! Function to get the atmosphere model of the body.
    /*!
     *  Function to get the atmosphere model of the body.
     *  \return Atmosphere model of the body.
     */
    boost::shared_ptr< aerodynamics::AtmosphereModel > getAtmosphereModel( )
    {
        return atmosphereModel_;
    }

    //! Function to get the rotation model of the body.
    /*!
     *  Function to get the rotation model of the body.
     *  \return Rotation model of the body.
     */
    boost::shared_ptr< ephemerides::RotationalEphemeris > getRotationalEphemeris( )
    {
        return rotationalEphemeris_;
    }

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > getShapeModel( )
    {
        return shapeModel_;
    }

    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
    getAerodynamicCoefficientInterface( )
    {
        return aerodynamicCoefficientInterface_;
    }

    //! Function to retrieve the body flight conditions
    /*!
     * Function to retrieve the body flight conditions, which calculates current aerodynamic angles,
     * altitude, etc.
     * \return Body flight conditions
     */
    boost::shared_ptr< aerodynamics::FlightConditions > getFlightConditions( )
    {
        return aerodynamicFlightConditions_;
    }

    //! Function to retrieve the shape model of the body.
    /*!
     *  Function to retrieve the shape model of the body.
     *  \return Shape model of the body.
     */
    std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > >
    getRadiationPressureInterfaces( )
    {
        return radiationPressureInterfaces_;
    }

    std::pair< bool, boost::shared_ptr< gravitation::GravityFieldVariations > >
            getGravityFieldVariation(
                const gravitation::BodyDeformationTypes& deformationType,
                const std::string identifier = "" )
    {
        return gravityFieldVariationSet_->getGravityFieldVariation( deformationType, identifier );
    }

    boost::shared_ptr< gravitation::GravityFieldVariationsSet > getGravityFieldVariationSet( )
    {
        return gravityFieldVariationSet_;
    }



    //! Function to set the function returning body mass as a function of time
    /*!
     * Function to set the function returning body mass as a function of time
     * \param bodyMassFunction Function returning body mass as a function of time
     */
    void setBodyMassFunction( const boost::function< double( const double ) > bodyMassFunction )
    {
        bodyMassFunction_ = bodyMassFunction;
    }

    //! Function to set the body mass as being constant (i.e. time-independent)
    /*!
     * Function to set the body mass as being constant (i.e. time-independent)
     * \param bodyMass New constant body mass
     */
    void setConstantBodyMass( const double bodyMass )
    {
        bodyMassFunction_ = boost::lambda::constant( bodyMass );
        currentMass_ = bodyMass;
    }

    //! Function to get the function returning body mass as a function of time
    /*!
     * Function to get the function returning body mass as a function of time
     * \return Function returning body mass as a function of time
     */
    boost::function< double( const double ) > getBodyMassFunction( )
    {
        return bodyMassFunction_;
    }

    //! Function to update the body mass to the current time
    /*!
     * Function to update the body mass to the current time, using the bodyMassFunction_ function
     * \param time Current time
     */
    void updateMass( const double time )
    {
        if( bodyMassFunction_ != NULL )
        {
            currentMass_ = bodyMassFunction_( time );
        }
        else
        {
            throw std::runtime_error( "Error when updating body mass, no mass function is set" );
        }
    }

    //! Function to retrieve the current body mass
    /*!
     * Function to retrieve the current body mass
     * \return Current body mass.
     */
    double getBodyMass( )
    {
        return currentMass_;
    }


    void updateConstantEphemerisDependentMemberQuantities( )
    {
        if( boost::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    gravityFieldModel_ ) != NULL )
        {
            //std::cerr<<"Error, time. dep. update disabled due to circular dependency"<<std::endl;
            boost::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        gravityFieldModel_ )->updateCorrectionFunctions( );
        }
    }

protected:

private:


    //! Current state.
    basic_mathematics::Vector6d currentState_;

    //! Current state with long double precision.
    Eigen::Matrix< long double, 6, 1 > currentLongState_;




    //! Function returning the state of this body's ephemeris origin w.r.t. the global origin
    //! (as set by setGlobalFrameBodyEphemerides function).
    boost::function< basic_mathematics::Vector6d( const double& ) >
            ephemerisFrameToBaseFrameFunction_;

    //! Function returning the state of this body's ephemeris origin w.r.t. the global origin
    //! (as set by setGlobalFrameBodyEphemerides function), with long double precision
    boost::function< Eigen::Matrix< long double, 6, 1 >( const double& ) >
            ephemerisFrameToBaseFrameLongFunction_;




    //! Current rotation from the global to the body-fixed frame.
    Eigen::Quaterniond currentRotationToLocalFrame_;

    //! Current first derivative w.r.t. time of the rotation matrix from the global to the
    //! body-fixed frame.
    Eigen::Matrix3d currentRotationToLocalFrameDerivative_;

    //! Current angular velocity vector for body's rotation, expressed in the global frame.
    Eigen::Vector3d currentAngularVelocityVectorInGlobalFrame_;


    //! Mass of body (default set to zero, calculated from GravityFieldModel when it is set).
    double currentMass_;

    //! Function returning body mass as a function of time.
    boost::function< double( const double ) > bodyMassFunction_;


    //! Ephemeris of body.
    boost::shared_ptr< ephemerides::Ephemeris > bodyEphemeris_;

    //! Gravity field model of body.
    boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

    boost::shared_ptr< gravitation::GravityFieldVariationsSet > gravityFieldVariationSet_;

    //! Atmosphere model of body.
    boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel_;

    //! Shape model of body.
    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel_;

    //! Aerodynamic coefficient model of body.
    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
            aerodynamicCoefficientInterface_;

    //! Object used for calculating current aerodynamic angles, altitude, etc.
    boost::shared_ptr< aerodynamics::FlightConditions > aerodynamicFlightConditions_;

    //! Rotation model of body.
    boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris_;

    //! List of radiation pressure models for the body, with the sources bodies as key
    std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > >
            radiationPressureInterfaces_;

    //! Predefined iterator for efficiency purposes.
    std::map< std::string,
              boost::shared_ptr< electro_magnetism::RadiationPressureInterface > >::iterator
    radiationPressureIterator_;
};

typedef std::unordered_map< std::string, boost::shared_ptr< Body > > NamedBodyMap;

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_BODY_H
