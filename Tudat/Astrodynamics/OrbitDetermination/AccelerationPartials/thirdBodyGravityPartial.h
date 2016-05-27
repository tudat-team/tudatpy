/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THIRDBODYGRAVITYPARTIAL_H
#define TUDAT_THIRDBODYGRAVITYPARTIAL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Function to get the third body acceleration type from the direct acceleration partial object.
/*!
 *  Function to get the third body acceleration type from the direct acceleration partial object.
 *  \param directGravityPartial Partial derivative of direct acceleration.
 *  \return Type of acceleration for third body acceleration for which direct acceleration partial is given by
 *  directGravityPartial (i.e. third_body_central_gravity if input is of type CentralGravitationPartial).
 */
template< typename DirectGravityPartial >
basic_astrodynamics::AvailableAcceleration getAccelerationTypeOfThirdBodyGravity(
        const boost::shared_ptr< DirectGravityPartial > directGravityPartial )
{
    using namespace basic_astrodynamics;
    AvailableAcceleration accelerationType;

    // Check type of direct partial derivative.
    if( boost::dynamic_pointer_cast< CentralGravitationPartial >( directGravityPartial ) != NULL )
    {
        accelerationType = third_body_central_gravity;
    }
    else
    {
        throw std::runtime_error( "Error when getting third body partial type, type not identified" );
    }
    return accelerationType;
}

//! Class to calculate the partials of a third-body gravitational acceleration w.r.t. parameters and states.
/*!
 *  Class to calculate the partials of a third-body gravitational acceleration w.r.t. parameters and states. This class may
 *  be used for any direct gravitational acceleration (central, spherical harmonic, mutual spherical harmonic, etc.),
 *  providin a generic third-body partial interface. The template parameter is the derived class of AccelerationPartial
 *  for the associated direct acceleration partial.
 */
template< typename DirectGravityPartial >
class ThirdBodyGravityPartial: public AccelerationPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param partialOfDirectGravityOnBodyUndergoingAcceleration Partial derivative of direct acceleration from
     * acceleratingBody on acceleratedBody.
     * \param partialOfDirectGravityOnCentralBody Partial derivative of direct acceleration from
     * centralBodyName on acceleratedBody.
     * \param acceleratedBody Name of body undergoing acceleration
     * \param acceleratingBody Name of body exerting acceleration (third-body)
     * \param centralBodyName Name of central body w.r.t. which the acceleration is computed.
     */
    ThirdBodyGravityPartial(
            const boost::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnBodyUndergoingAcceleration,
            const boost::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnCentralBody,
            const std::string& acceleratedBody, const std::string& acceleratingBody,
            const std::string& centralBodyName ):
        AccelerationPartial( acceleratedBody, acceleratingBody, getAccelerationTypeOfThirdBodyGravity(
                                 partialOfDirectGravityOnBodyUndergoingAcceleration ) ),
        partialOfDirectGravityOnBodyUndergoingAcceleration_( partialOfDirectGravityOnBodyUndergoingAcceleration ),
        partialOfDirectGravityOnCentralBody_( partialOfDirectGravityOnCentralBody ), centralBodyName_( centralBodyName ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAcceleratedBody(
                    partialMatrix, addContribution, startRow, startColumn );

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( acceleratedBody_ ) == 1 )
        {
            partialOfDirectGravityOnCentralBody_->wrtPositionOfAdditionalBody(
                        acceleratedBody_, partialMatrix, addContribution, startRow, startColumn   );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. velocity of body undergoing acceleration.
     */
    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAcceleratedBody(
                    partialMatrix, addContribution, startRow, startColumn  );

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( acceleratedBody_ ) == 1 )
        {
            partialOfDirectGravityOnCentralBody_->wrtVelocityOfAdditionalBody(
                        acceleratedBody_, partialMatrix, addContribution, startRow, startColumn  );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body exerting acceleration.
     */
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )

    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAcceleratingBody(
                    partialMatrix, addContribution, startRow, startColumn  );
        partialOfDirectGravityOnCentralBody_->wrtPositionOfAcceleratingBody(
                    partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn  );

    }


    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. velocity of body exerting acceleration.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAcceleratingBody(
                    partialMatrix, addContribution, startRow, startColumn );
        partialOfDirectGravityOnCentralBody_->wrtVelocityOfAcceleratingBody(
                    partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn );
    }

    void wrtPositionOfAdditionalBody( const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                      const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( bodyName == centralBodyName_ )
        {
            partialOfDirectGravityOnCentralBody_->wrtPositionOfAcceleratedBody(
                        partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn  );
        }

        if( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAdditionalBody(
                        bodyName, partialMatrix, addContribution, startRow, startColumn  );
        }

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            partialOfDirectGravityOnCentralBody_->wrtPositionOfAdditionalBody(
                        bodyName, partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn );
        }
    }

    void wrtVelocityOfAdditionalBody( const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                      const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        if( bodyName == centralBodyName_ )
        {
            partialOfDirectGravityOnCentralBody_->wrtVelocityOfAcceleratedBody(
                        partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn  );
        }

        if( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAdditionalBody(
                        bodyName, partialMatrix, addContribution, startRow, startColumn );
        }

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            partialOfDirectGravityOnCentralBody_->wrtVelocityOfAdditionalBody(
                        bodyName, partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn );
        }
    }

    bool isAccelerationPartialWrtAdditionalBodyNonNull( const std::string& bodyName )
    {
        bool isAccelerationDependentOnBody = 0;
        if( bodyName == centralBodyName_ )
        {
            isAccelerationDependentOnBody = 1;
        }

        if( ( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) == 1 ) ||
                ( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) == 1 ) )
        {
            isAccelerationDependentOnBody = 1;
        }

        return isAccelerationDependentOnBody;
    }

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );

        return createMergedParameterPartialFunction( partialFunctionFromDirectGravity,
                                                     partialFunctionFromCentralGravity );
    }

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );

        return createMergedParameterPartialFunction( partialFunctionFromDirectGravity,
                                                     partialFunctionFromCentralGravity );
    }

    int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromDirectGravity.second > 0 )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );
        if( partialFunctionFromCentralGravity.second > 0 )
        {
            partialOfDirectGravityOnCentralBody_->setParameterPartialUpdateFunction( parameter );
        }

        if( partialFunctionFromCentralGravity.second > 0 || partialFunctionFromDirectGravity.second > 0 )
        {
            parameterDoublePartialFunctions_[ parameter ] =
                    getCombinedCurrentDoubleParameterFunction(
                        partialOfDirectGravityOnBodyUndergoingAcceleration_,
                        partialOfDirectGravityOnCentralBody_,
                        parameter, partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second, 1 );
            isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
            currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
        }
        return std::max( partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second );
    }


    int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromDirectGravity.second > 0 )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );
        if( partialFunctionFromCentralGravity.second > 0 )
        {
            partialOfDirectGravityOnCentralBody_->setParameterPartialUpdateFunction( parameter );
        }

        if( partialFunctionFromCentralGravity.second > 0 || partialFunctionFromDirectGravity.second > 0 )
        {
            parameterVectorPartialFunctions_[ parameter ] =
                    getCombinedCurrentVectorParameterFunction(
                        partialOfDirectGravityOnBodyUndergoingAcceleration_,
                        partialOfDirectGravityOnCentralBody_,
                        parameter, partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second, 1 );
            isCurrentVectorParameterPartialSet_[ parameter ] = 0;
            currentVectorParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, parameter->getParameterSize( ) );

        }
        return std::max( partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second );
    }


    void update( const double currentTime )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->update( currentTime );
        partialOfDirectGravityOnCentralBody_->update( currentTime );

        currentTime_ = currentTime;
    }

    boost::shared_ptr< DirectGravityPartial > getPartialOfDirectGravityOnBodyUndergoingAcceleration( )
    {
        return partialOfDirectGravityOnBodyUndergoingAcceleration_;
    }

    boost::shared_ptr< DirectGravityPartial > getPartialOfDirectGravityOnCentralBody( )
    {
        return partialOfDirectGravityOnCentralBody_;
    }

    std::string getCentralBodyName( )
    {
        return centralBodyName_;
    }


protected:


    void resetTimeOfMemberObjects( )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->resetTime( currentTime_ );
        partialOfDirectGravityOnCentralBody_->resetTime( currentTime_ );
    }

    void updateParameterPartialsOfMemberObjects( )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->updateParameterPartials( );
        partialOfDirectGravityOnCentralBody_->updateParameterPartials( );
    }

private:
    boost::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnBodyUndergoingAcceleration_;

    boost::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnCentralBody_;

    std::string centralBodyName_;

};

}

}

}

#endif // TUDAT_THIRDBODYGRAVITYPARTIAL_H
