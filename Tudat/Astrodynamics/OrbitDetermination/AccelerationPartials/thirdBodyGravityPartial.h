#ifndef THIRDBODYGRAVITYPARTIAL_H
#define THIRDBODYGRAVITYPARTIAL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

template< typename DirectGravityPartial >
basic_astrodynamics::AvailableAcceleration getAccelerationTypeOfThirdBodyGravity(
        const boost::shared_ptr< DirectGravityPartial > directGravityPartial )
{
    using namespace basic_astrodynamics;
    AvailableAcceleration accelerationType;

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

template< typename DirectGravityPartial >
class ThirdBodyGravityPartial: public AccelerationPartial
{
public:
    ThirdBodyGravityPartial( const boost::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnBodyUndergoingAcceleration,
                             const boost::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnCentralBody,
                             const std::string& acceleratedBody, const std::string& acceleratingBody,
                             const std::string& centralBodyName ):
        AccelerationPartial( acceleratedBody, acceleratingBody, getAccelerationTypeOfThirdBodyGravity(
                                 partialOfDirectGravityOnBodyUndergoingAcceleration ) ),
        partialOfDirectGravityOnBodyUndergoingAcceleration_( partialOfDirectGravityOnBodyUndergoingAcceleration ),
        partialOfDirectGravityOnCentralBody_( partialOfDirectGravityOnCentralBody ), centralBodyName_( centralBodyName ){ }

    Eigen::Matrix3d wrtPositionOfAcceleratedBody( )
    {
        Eigen::Matrix3d partial = partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAcceleratedBody( );

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( acceleratedBody_ ) == 1 )
        {
            partial += partialOfDirectGravityOnCentralBody_->wrtPositionOfAdditionalBody( acceleratedBody_ );
        }

        return partial;
    }

    Eigen::Matrix3d wrtVelocityOfAcceleratedBody( )
    {
        Eigen::Matrix3d partial = partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAcceleratedBody( );

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( acceleratedBody_ ) == 1 )
        {
            partial += partialOfDirectGravityOnCentralBody_->wrtVelocityOfAdditionalBody( acceleratedBody_ );
        }

        return partial;
    }

    Eigen::Matrix3d wrtPositionOfAcceleratingBody( )
    {
        Eigen::Matrix3d partial = partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAcceleratingBody( ) -
                partialOfDirectGravityOnCentralBody_->wrtPositionOfAcceleratingBody( );
        return partial;

    }

    Eigen::Matrix3d wrtVelocityOfAcceleratingBody( )
    {
        return partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAcceleratingBody( ) -
                partialOfDirectGravityOnCentralBody_->wrtVelocityOfAcceleratingBody( );
    }

    Eigen::Matrix3d wrtPositionOfAdditionalBody( const std::string& bodyName )
    {
        Eigen::Matrix3d positionPartial = Eigen::Matrix3d::Zero( );

        if( bodyName == centralBodyName_ )
        {
            positionPartial -= partialOfDirectGravityOnCentralBody_->wrtPositionOfAcceleratedBody( );
        }

        if( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            positionPartial += partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAdditionalBody( bodyName );
        }

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            positionPartial -= partialOfDirectGravityOnCentralBody_->wrtPositionOfAdditionalBody( bodyName );
        }

        return positionPartial;
    }

    Eigen::Matrix3d wrtVelocityOfAdditionalBody( const std::string& bodyName )
    {
        Eigen::Matrix3d velocityPartial = Eigen::Matrix3d::Zero( );

        if( bodyName == centralBodyName_ )
        {
            velocityPartial -= partialOfDirectGravityOnCentralBody_->wrtVelocityOfAcceleratedBody( );
        }

        if( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            velocityPartial += partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAdditionalBody( bodyName );
        }

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonNull( bodyName ) )
        {
            velocityPartial -= partialOfDirectGravityOnCentralBody_->wrtVelocityOfAdditionalBody( bodyName );
        }

        return velocityPartial;
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

    std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );

        return createMergedParameterPartialFunction( partialFunctionFromDirectGravity,
                                                     partialFunctionFromCentralGravity );
    }

    std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );

        return createMergedParameterPartialFunction( partialFunctionFromDirectGravity,
                                                     partialFunctionFromCentralGravity );
    }

    int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromDirectGravity.second > 0 )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromCentralGravity =
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
            doesCurrentDoubleParameterPartialExist_[ parameter ] = 0;
        }
        return std::max( partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second );
    }


    int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromDirectGravity.second > 0 )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionFromCentralGravity =
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
            doesCurrentVectorParameterPartialExist_[ parameter ] = 0;

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
        //std::cout<<"Updating parameter partials (third body)"<<std::endl;

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

#endif // THIRDBODYGRAVITYPARTIAL_H
