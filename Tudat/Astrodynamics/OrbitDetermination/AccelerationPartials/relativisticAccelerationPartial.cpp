
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/relativisticAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

Eigen::Matrix3d computePartialOfSchwardschildAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d relativeState, Eigen::Vector3d currentAcceleration,
        const double gravitationalParameter, const double ppnParameterGamma, const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    Eigen::Matrix3d partial = 2.0 * ( ppnParameterGamma + ppnParameterBeta ) * gravitationalParameter / distance * (
                Eigen::Matrix3d::Identity( ) - position * position.transpose( ) / ( distance * distance ) );
    partial -= ppnParameterGamma * velocity.dot( velocity ) * Eigen::Matrix3d::Identity( );
    partial += 2.0 * ( 1.0 + ppnParameterGamma ) * velocity * velocity.transpose( );
    partial *= gravitationalParameter *  physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance );
    partial -= 3.0 * currentAcceleration * position.transpose( ) / ( distance * distance );
    return  partial;
}

Eigen::Matrix3d computePartialOfSchwardschildAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter, const double ppnParameterGamma )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    Eigen::Matrix3d partial = gravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( - 2.0 * ppnParameterGamma * position * velocity.transpose( ) +
              2.0 * ( 1.0 + ppnParameterGamma ) * (
                  position.dot( velocity ) * Eigen::Matrix3d::Identity( ) + velocity * position.transpose( ) ) );

    return partial;
}

void computePartialOfSchwardschildAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partial, const double ppnParameterGamma, const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partial = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( -ppnParameterGamma * ( velocity.dot( velocity ) ) * position +
              2.0 * ( 1.0 + ppnParameterGamma ) *
              ( position.dot( velocity ) ) * velocity
              + 4.0 * gravitationalParameter * ( ppnParameterGamma + ppnParameterBeta ) *
              position / distance );
}

Eigen::Vector3d computePartialOfSchwardschildAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    return physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) * (
                ( 2.0 * gravitationalParameter / distance - velocity.dot( velocity ) ) * position + 2.0 * position.dot( velocity ) * velocity );
}

Eigen::Vector3d computePartialOfSchwardschildAccelerationCorrectionWrtPpnParameterBeta(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    double distance = position.norm( );

    return physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) * (
                ( 2.0 * gravitationalParameter / distance ) * position );
}

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
RelativisticAccelerationPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;
    if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::gravitational_parameter:

            partialFunction = boost::bind( &RelativisticAccelerationPartial::wrtGravitationalParameterOfCentralBody,
                                           this, _1 );
            numberOfRows = 1;

            break;
        default:
            break;
        }
    }
//    else if( parameter->getParameterName( ).second.first == "global_metric"  )
//    {
//        switch( parameter->getParameterName( ).first )
//        {
//        case estimatable_parameters::ppn_parameter_gamma:

//            std::cout<<"PPN gamma "<<std::endl;
//            partialFunction = boost::bind( &RelativisticAccelerationPartial::wrtPpnParameterGamma,
//                                           this );
//            numberOfRows = 1;

//            break;
//        case estimatable_parameters::ppn_parameter_beta:

//            std::cout<<"PPN beta "<<std::endl;
//            partialFunction = boost::bind( &RelativisticAccelerationPartial::wrtPpnParameterBeta,
//                                           this );
//            numberOfRows = 1;

//            break;
//        default:
//            break;
//        }
//    }
    return std::make_pair( partialFunction, numberOfRows );
}

void RelativisticAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRelativeState_ = ( acceleratedBodyState_( ) - centralBodyState_( ) );
        currentAcceleration_ = currentAccelerationFunction_( );
        currentPartialWrtPosition_ = computePartialOfSchwardschildAccelerationCorrectionWrtPosition(
                    currentRelativeState_, currentAcceleration_, centralBodyGravitationalParameterFunction_( ),
                    ppnGammaParameterFunction_( ), ppnBetaParameterFunction_( ) );

        currentPartialWrtVelocity_ = computePartialOfSchwardschildAccelerationCorrectionWrtVelocity(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ),
                    ppnGammaParameterFunction_( ) );
        currentTime_ = currentTime;
    }
}

}

}


