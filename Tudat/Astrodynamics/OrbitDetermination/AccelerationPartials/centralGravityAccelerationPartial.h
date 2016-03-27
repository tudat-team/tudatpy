#ifndef CENTRALGRAVITYACCELERATIONPARTIALS_H
#define CENTRALGRAVITYACCELERATIONPARTIALS_H

#include <iostream>

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

using namespace gravitation;

//! Calculates the partial derivative of a point mass gravitational acceleration wrt the position of the body being accelerated.
/*!
 *  Calculates the partial derivative of a point mass gravitational acceleration wrt the position
 *  of the body being accelerated.
 *  \param stateOfGravitatedBody Cartesian state of body being accelerated.
 *  \param stateOfGravitatingBody Cartesian state of body exerting acceleration.
 *  \param gravitationalParameter Gravitational parameter of gravitating body.
 *  \return Matrix with the Jacobian of the acceleration vector w.r.t. the position vector.
 */
Eigen::Matrix3d calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
        const Eigen::Vector3d stateOfGravitatedBody,
        const Eigen::Vector3d stateOfGravitatingBody,
        double gravitationalParameter );

//! Calculates the partial derivative of a point mass gravitational acceleration wrt the gravitational parameter of the central body
/*!
 *  Calculates the partial derivative of a point mass gravitational acceleration wrt the gravitational parameter of the central body
 *  \param stateOfGravitatedBody Cartesian state of body being accelerated.
 *  \param stateOfGravitatingBody Cartesian state of body exerting acceleration.
 *  \param gravitationalParameter Gravitational parameter of gravitating body.
 *  \return Vector with the partial of the acceleration vector w.r.t. ational parameter of the central body.
 */
Eigen::Vector3d computePartialOfCentralGravityWrtGravitationalParameter( const Eigen::Vector3d acceleratedBodyPosition,
                                                                         const Eigen::Vector3d acceleratingBodyPosition);


//! Class to calculate the partials of teh central gravitational acceleration w.r.t. parameters and states.
class CentralGravitationPartial: public AccelerationPartial
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param gravitationalAcceleration Central gravitational acceleration w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    CentralGravitationPartial(
            const boost::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > gravitationalAcceleration,
            const std::string acceleratedBody,
            const std::string acceleratingBody );

    CentralGravitationPartial(
            const boost::shared_ptr< CentralGravitationPartial > originalAccelerationPartial );


    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtPositionOfAcceleratedBody( )
    {
        return currentPositionPartial;
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtVelocityOfAcceleratedBody( )
    {
        return Eigen::Matrix3d::Zero( );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtPositionOfAcceleratingBody( )
    {
        return -currentPositionPartial;
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtVelocityOfAcceleratingBody( )
    {
        return Eigen::Matrix3d::Zero( );
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for central gravity acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in return function (0 for no dependency, 1 otherwise).
     */
    virtual std::pair< boost::function< Eigen::MatrixXd( ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    virtual std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< Eigen::MatrixXd( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for updating partial to current state.
    /*!
     *  Function for updating partial to current state (computes position partial)
     */
    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentPositionPartial = calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                        acceleratedBodyState_( ),
                        centralBodyState_( ),
                        gravitationalParameterFunction_( ) );
            currentTime_ = currentTime;
        }
    }

    boost::function< double( ) > getGravitationalParameterFunction( )
    {
        return gravitationalParameterFunction_;
    }

    bool getAccelerationUsesMutualAttraction( )
    {
        return accelerationUsesMutualAttraction_;
    }

    boost::function< Eigen::Vector3d( ) > getPositionFunctionOfBodyExertingAcceleration( )
    {
        return centralBodyState_;
    }

    boost::function< Eigen::Vector3d( ) > getPositionFunctionOfBodyUndergoingAcceleration( )
    {
        return acceleratedBodyState_;
    }

    Eigen::Matrix3d getCurrentPositionPartial( )
    {
        return currentPositionPartial;
    }


protected:

    std::pair< boost::function< Eigen::MatrixXd( ) >, int > getGravitationalParameterPartialFunction(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterId );

    //! Function to calculate central gravity partial w.r.t. central body gravitational parameter
    /*!
     *  Function to calculate central gravity partial w.r.t. central body gravitational parameter
     */
    Eigen::Vector3d wrtGravitationalParameterOfCentralBody( );

    //! Function to retrieve current gravitational parameter of central body.
    /*!
     *  Function to retrieve current gravitational parameter of central body.
     */
    boost::function< double( ) > gravitationalParameterFunction_;

    //! Function to retrieve current state of body exerting acceleration.
    /*!
     *  Function to retrieve current state of body exerting acceleration.
     */
    boost::function< Eigen::Vector3d( ) > centralBodyState_;

    //! Function to retrieve current state of body undergoing acceleration.
    /*!
     *  Function to retrieve current state of body undergoing acceleration.
     */
    boost::function< Eigen::Vector3d( ) > acceleratedBodyState_;

    bool accelerationUsesMutualAttraction_;




    //! Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
    /*!
     *  Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. position of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentPositionPartial;

};

}

}

}
#endif // CENTRALGRAVITYACCELERATIONPARTIALS_H
