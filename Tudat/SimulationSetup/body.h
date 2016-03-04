/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      121030    K. Kumar          File created.
 *      130225    K. Kumar          Updated include-guard and namespace names; updated Vector6d
 *                                  references to use Tudat definition.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT__BODY_H
#define TUDAT__BODY_H

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
     * Constructor for a body, sets current time, state, rotation and mass values
     * (all with default parameters). The input state is used internally to
     * set the current position (taken as a segment of the input state given by the indices
     * (0, 3)) and the current velocity (taken as a segment of the input state given by the indices
     * (3, 3).
     * \param state Current state of body at initialization (default = zeroes).
     * \param time Current time of body at initialization (default = zeroes).
     * \param bodyMass Current mass of body at initialization (default = zeroes).
     * \param currentRotationToGlobalFrame Current rotation of from body-fixed to inertial frames
     *  at initialization (default = identity)
     */
    Body( const basic_mathematics::Vector6d& state =
            basic_mathematics::Vector6d::Zero( ),
          const double time = 0.0, const double bodyMass = 0.0,
          const Eigen::Quaterniond currentRotationToLocalFrame =
            Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) )
        : currentState( state ),
          currentTime( time ),
          currentRotationToLocalFrame_( currentRotationToLocalFrame ),
          bodyMass_( bodyMass ),
          ephemerisFrameToBaseFrameFunction_( boost::lambda::constant( basic_mathematics::Vector6d::Zero( ) ) ),
          ephemerisFrameToBaseFrameLongFunction_( boost::lambda::constant( Eigen::Matrix< long double, 6, 1 >::Zero( ) ) )
    { }

      void setState( const basic_mathematics::Vector6d& state )
    {
        currentState = state; // Must be in global frame.
        //currentLongState = state.cast< long double >( );
    }

    void setLongState( const Eigen::Matrix< long double, 6, 1 >& longState )
    {
        currentLongState = longState; // Must be in global frame.
        currentState = longState.cast< double >( );
    }

    template< typename ScalarStateType >
    void setTemplatedState( const Eigen::Matrix< ScalarStateType, 6, 1 >& state );

    basic_mathematics::Vector6d getStateInBaseFrameFromEphemeris( const double& time )
    {
        return bodyEphemeris_->getCartesianStateFromEphemeris( time ) + ephemerisFrameToBaseFrameFunction_( time );
    }

    Eigen::Matrix< long double, 6, 1 > getLongStateInBaseFrameFromEphemeris( const double& time )
    {
        return bodyEphemeris_->getCartesianLongStateFromEphemeris( time ) + ephemerisFrameToBaseFrameLongFunction_( time );
    }

    void setStateFromEphemeris( const double& time )
    {
        currentState = bodyEphemeris_->getCartesianStateFromEphemeris( time ) + ephemerisFrameToBaseFrameFunction_( time );
        //currentLongState = state.cast< long double >( );

    }

    template< typename StateScalarType, typename TimeType >
    void setTemplatedStateFromEphemeris( const TimeType& time );

    template< typename ScalarStateType, typename TimeType >
    Eigen::Matrix< ScalarStateType, 6, 1 > getTemplatedStateInBaseFrameFromEphemeris( const TimeType& time );

    void setCurrentRotationToLocalFrameFromEphemeris( const double time )
    {
        currentRotationToLocalFrame_ = rotationalEphemeris_->getRotationToTargetFrame( time );
    }

    void setCurrentRotationToLocalFrameDerivativeFromEphemeris( const double time )
    {
        currentRotationToLocalFrameDerivative_ = rotationalEphemeris_->getDerivativeOfRotationToTargetFrame( time );
    }

    void setCurrentAngularVelocityVectorInGlobalFrame( const double time )
    {
        currentAngularVelocityVectorInGlobalFrame_ = rotationalEphemeris_->getRotationalVelocityVectorInBaseFrame( time );
    }

    void setCurrentRotationalStateToLocalFrameFromEphemeris( const double time )
    {
        rotationalEphemeris_->getFullRotationalQuantitiesToTargetFrame(
                    currentRotationToLocalFrame_, currentRotationToLocalFrameDerivative_, currentAngularVelocityVectorInGlobalFrame_, time );
        //setCurrentRotationToLocalFrameFromEphemeris( time );
        //setCurrentRotationToLocalFrameDerivativeFromEphemeris( time );
        //setCurrentAngularVelocityVectorInGlobalFrame( time );
    }

    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    basic_mathematics::Vector6d getState( ) { return currentState; }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getPosition( ) { return currentState.segment( 0, 3 ); }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getVelocity( ) { return currentState.segment( 3, 3 ); }

    //! Get current time.
    /*!
     * Returns the internally stored current time.
     * \return Current time.
     */
    double getCurrentTime( ) { return currentTime; }

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

    void setBaseFrameFunction( boost::function< basic_mathematics::Vector6d( const double& ) > ephemerisFrameToBaseFrameFunction )
    {
        ephemerisFrameToBaseFrameFunction_ = ephemerisFrameToBaseFrameFunction;
    }

    void setBaseFrameLongFunction( boost::function< Eigen::Matrix< long double, 6, 1 >( const double& ) > ephemerisFrameToBaseFrameLongFunction )
    {
        ephemerisFrameToBaseFrameLongFunction_ = ephemerisFrameToBaseFrameLongFunction;
    }

    //! Function to set the gravity field of the body.
    /*!
     *  Function to set the gravity field of the body; input is also used to (re)set the mass
     *  of the body.
     *  \param gravityFieldModel New gravity field of the body.
     */
    void setGravityFieldModel(
            const boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel )
    {
        gravityFieldModel_ = gravityFieldModel;
        bodyMass_ = gravityFieldModel_->getGravitationalParameter( );
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

    void setShapeModel( const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel )
    {
        shapeModel_ = shapeModel;
    }

    void setAerodynamicCoefficientInterface(
            const boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
            aerodynamicCoefficientInterface)
    {
        aerodynamicCoefficientInterface_ = aerodynamicCoefficientInterface;
    }

    void setFlightConditions(
            const boost::shared_ptr< aerodynamics::FlightConditions > aerodynamicFlightConditions )
    {
       aerodynamicFlightConditions_ = aerodynamicFlightConditions;
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

    boost::shared_ptr< aerodynamics::FlightConditions > getFlightConditions( )
    {
       return aerodynamicFlightConditions_;
    }


    //! Get current rotation from body-fixed to inertial frame.
    /*!
     *  Get current rotation from body-fixed to inertial frame, as set from the rotationalEphemeris_
     *  by the setCurrentTimeAndState function. If body has no rotational ephemeris, an identity
     *  quaternion (no rotation) is returned.
     *  \return Current rotation from body-fixed to inertial frame
     */
    Eigen::Quaterniond getCurrentRotationToGlobalFrame( )
    {
        return currentRotationToLocalFrame_.inverse( );
    }

    //! Get current rotation from inertial to body-fixed frame.
    /*!
     *  Get current rotation from inertial to body-fixed frame, as set from the rotationalEphemeris_
     *  by the setCurrentTimeAndState function. If body has no rotational ephemeris, an identity
     *  quaternion (no rotation) is returned.
     *  \return Current rotation from inertial to body-fixed frame
     */
    Eigen::Quaterniond getCurrentRotationToLocalFrame( )
    {
        return currentRotationToLocalFrame_;
    }

    Eigen::Matrix3d getCurrentRotationMatrixDerivativeToGlobalFrame( )
    {
        return currentRotationToLocalFrameDerivative_.inverse( );
    }


    Eigen::Matrix3d getCurrentRotationMatrixDerivativeToLocalFrame( )
    {
        return currentRotationToLocalFrameDerivative_;
    }

    double getBodyMass( )
    {
        return bodyMass_;
    }

    void updateMass( const double bodyMass )
    {
        bodyMass_ = bodyMass;
    }

    std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > >
    getRadiationPressureInterfaces( )
    {
        return radiationPressureInterfaces_;
    }

    void setRadiationPressureInterface(
            const std::string radiatingBody,
            const boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface )
    {
        radiationPressureInterfaces_[ radiatingBody ] = radiationPressureInterface;
    }



protected:

private:

    //! Current state.
    basic_mathematics::Vector6d currentState;

    Eigen::Matrix< long double, 6, 1 > currentLongState;

    //! Current time.
    double currentTime;

    Eigen::Quaterniond currentRotationToLocalFrame_;

    Eigen::Matrix3d currentRotationToLocalFrameDerivative_;

    Eigen::Vector3d currentAngularVelocityVectorInGlobalFrame_;

    //! Mass of body (default set to zero, calculated from GravityFieldModel when it is set).
    double bodyMass_;

    boost::function< basic_mathematics::Vector6d( const double& ) > ephemerisFrameToBaseFrameFunction_;

    boost::function< Eigen::Matrix< long double, 6, 1 >( const double& ) > ephemerisFrameToBaseFrameLongFunction_;

    //! Ephemeris of body.
    boost::shared_ptr< ephemerides::Ephemeris > bodyEphemeris_;

    //! Gravity field model of body.
    boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

    //! Atmosphere model of body.
    boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel_;

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel_;

    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficientInterface_;

    //! Body's aerodynamic coefficient interface object.
    boost::shared_ptr< aerodynamics::FlightConditions > aerodynamicFlightConditions_;

    //! Rotation model of body.
    boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris_;

    std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > > radiationPressureInterfaces_;

    std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > >::iterator radiationPressureIterator_;
};

typedef std::map< std::string, boost::shared_ptr< Body > > NamedBodyMap;


}

}

#endif // SATELLITE_PROPAGATOR_EXAMPLES_BODY_H
