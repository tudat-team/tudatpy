#ifndef PARTIALS_H
#define PARTIALS_H

#include <string>
#include <map>
#include <Eigen/Core>

#include <boost/bind.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


//! Base class for objects calculating partial derivatives of accelerations, as used in orbit determination.
/*!
 *  Base class for objects calculating partial derivatives of accelerations, as used in orbit determination.
 *  Derived classes implement derivative-calculating models for specific acceleration models, so that the calculation
 *  of all partials of a single type acceleration model is encompassed in a single derived class.
 */
class AccelerationPartial: public StateDerivativePartial
{

public:
    //! Base class constructor.
    /*!
     *  Constructor of base class, sets the base class member variables identifying the body undergoing and exerting the
     *  acceleration.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     *  \param accelerationType Type of acceleration w.r.t. which partial is taken.
     */
    AccelerationPartial( const std::string& acceleratedBody, const std::string& acceleratingBody,
                         const basic_astrodynamics::AvailableAcceleration accelerationType ):
        StateDerivativePartial( propagators::transational_state, std::make_pair( acceleratedBody, "" ) ),
        acceleratedBody_( acceleratedBody ), acceleratingBody_( acceleratingBody ),accelerationType_( accelerationType ) { }

    //! Base class destructor.
    /*!
     *  Base class destructor.
     */
    virtual ~AccelerationPartial( ) { }

    std::pair< boost::function< Eigen::MatrixXd( ) >, int >  getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunction = std::make_pair( boost::function< Eigen::MatrixXd( ) >( ), 0 );

        if( integratedStateType == propagators::transational_state )
        {
            if( stateReferencePoint.second != "" )
            {
                std::cerr<<"Error when getting state derivative partial acceleration model, cannot have reference point on body for dynamics"<<std::endl;
            }
            else if( stateReferencePoint.first == acceleratedBody_ )
            {
                partialFunction = std::make_pair( boost::bind( &AccelerationPartial::wrtStateOfAcceleratedBody, this ), 3 );
            }
            else if( stateReferencePoint.first == acceleratingBody_ )
            {
                partialFunction = std::make_pair( boost::bind( &AccelerationPartial::wrtStateOfAcceleratingBody, this ), 3 );
            }
            else if( isAccelerationPartialWrtAdditionalBodyNonNull( stateReferencePoint.first ) )
            {
                partialFunction = std::make_pair( boost::bind( &AccelerationPartial::wrtStateOfAdditionalBody, this, stateReferencePoint.first ), 3 );
            }
            else
            {
                partialFunction = std::make_pair( boost::function< Eigen::MatrixXd( ) >( ), 0 );
            }
        }

        return partialFunction;
    }

    bool isStateDerivativeDependentOnIntegratedState(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        bool isDependent = 0;
        if( integratedStateType == propagators::transational_state )
        {
            if( stateReferencePoint.second != "" )
            {
                std::cerr<<"Error when checking state derivative partial dependency of acceleration model, cannot have reference point on body for dynamics"<<std::endl;
            }
            else if( stateReferencePoint.first == acceleratedBody_ || stateReferencePoint.first == acceleratingBody_ ||
                     isAccelerationPartialWrtAdditionalBodyNonNull( stateReferencePoint.first ) )
            {
                isDependent = 1;
            }
        }
        return isDependent;
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body.
     *  Variables required for its calculation should be set (as shared pointers) in derived class.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    virtual Eigen::Matrix3d wrtPositionOfAcceleratedBody( ) = 0;

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the accelerated body.
     *  Variables required for its calculation should be set (as shared pointers) in derived class.
     *  \return Partial derivative of acceleration w.r.t. velocity of body undergoing acceleration.
     */
    virtual Eigen::Matrix3d wrtVelocityOfAcceleratedBody( ) = 0;

    Eigen::Matrix< double, 3, 6 > wrtStateOfAcceleratedBody( )
    {
        return ( Eigen::Matrix< double, 3, 6 >( )<<wrtPositionOfAcceleratedBody( ), wrtVelocityOfAcceleratedBody( ) ).finished( );
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the body exerting acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the body exerting acceleration.
     *  Variables required for its calculation should be set (as shared pointers) in derived class.
     *  \return Partial derivative of acceleration w.r.t. position of body exerting acceleration.
     */
    virtual Eigen::Matrix3d wrtPositionOfAcceleratingBody( ) = 0;

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the body exerting acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the body exerting acceleration.
     *  Variables required for its calculation should be set (as shared pointers) in derived class.
     *  \return Partial derivative of acceleration w.r.t. velocity of body exerting acceleration.
     */
    virtual Eigen::Matrix3d wrtVelocityOfAcceleratingBody( ) = 0;

    Eigen::Matrix< double, 3, 6 > wrtStateOfAcceleratingBody( )
    {
        return ( Eigen::Matrix< double, 3, 6 >( )<<wrtPositionOfAcceleratingBody( ), wrtVelocityOfAcceleratingBody( ) ).finished( );
    }

    virtual Eigen::Matrix3d wrtPositionOfAdditionalBody( const std::string& bodyName )
    {
        return Eigen::Matrix3d::Zero( );
    }

    virtual Eigen::Matrix3d wrtVelocityOfAdditionalBody( const std::string& bodyName )
    {
        return Eigen::Matrix3d::Zero( );
    }

    Eigen::Matrix< double, 3, 6 > wrtStateOfAdditionalBody( const std::string& bodyName )
    {
        return ( Eigen::Matrix< double, 3, 6 >( )<<wrtPositionOfAdditionalBody( bodyName ), wrtVelocityOfAdditionalBody( bodyName ) ).finished( );
    }

    virtual bool isAccelerationPartialWrtAdditionalBodyNonNull( const std::string& bodyName )
    {
        return 0;
    }

    //! Function to retrieve the name of the body undergoing acceleration
    /*!
     *  Function to retrieve the name of the body undergoing acceleration
     */
    std::string getAcceleratedBody( ) { return acceleratedBody_; }

    //! Function to retrieve the name of the body exerting acceleration
    /*!
     *  Function to retrieve the name of the body exerting acceleration
     */
    std::string getAcceleratingBody( ) { return acceleratingBody_; }

    basic_astrodynamics::AvailableAcceleration getAccelerationType( )
    {
        return accelerationType_;
    }

protected:    
    //! Name of the body undergoing acceleration
    /*!
     *  Name of the body undergoing acceleration
     */
    std::string acceleratedBody_;

    //! Name of the body exerting acceleration
    /*!
     *  Name of the body exerting acceleration
     */
    std::string acceleratingBody_;

    //! Type of acceleration w.r.t. which partial is taken.
    /*!
     *  Type of acceleration w.r.t. which partial is taken.
     */
    basic_astrodynamics::AvailableAcceleration accelerationType_;
};


}

}

}
#endif // PARTIALS_H
