#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/empiricalAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/empiricalAccelerationCoefficients.h"

namespace tudat
{

namespace acceleration_partials
{

using namespace gravitation;

//! Function determine the numerical partial derivative of the true anomaly wrt the elements of the Cartesian state
Eigen::Matrix< double, 1, 6 > calculateNumericalPartialOfTrueAnomalyWrtState(
        const Eigen::Vector6d& cartesianElements, const double gravitationalParameter,
        const Eigen::Vector6d& perturbations )
{
    using namespace tudat::orbital_element_conversions;

    // Initialize partial to zero
    Eigen::Matrix< double, 1, 6 > partial = Eigen::Matrix< double, 1, 6 >::Zero( );

    // Iterate over all six elements and calculate partials.
    Eigen::Vector6d perturbedCartesianElements;
    double upPerturbedTrueAnomaly, downPerturbedTrueAnomaly;
    for( int i = 0; i < 6; i++ )
    {
        // Calculate true anomaly at up-perturbed cartesian element i
        perturbedCartesianElements = cartesianElements;
        perturbedCartesianElements( i ) += perturbations( i );
        upPerturbedTrueAnomaly = convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter )( trueAnomalyIndex );

        // Check validity of result
        if( upPerturbedTrueAnomaly != upPerturbedTrueAnomaly )
        {
            std::cout<<"Error 1 in partial of true anomaly wrt cartesian state, element"<<i<<
                       ", keplerian state: "<<std::endl<<
                       convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter ).transpose( )<<std::endl<<
                       perturbedCartesianElements.transpose( )<<std::endl;
        }

        // Calculate true anomaly at down-perturbed cartesian element i
        perturbedCartesianElements = cartesianElements;
        perturbedCartesianElements( i ) -= perturbations( i );
        downPerturbedTrueAnomaly = convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter )( trueAnomalyIndex );

        // Check validity of result
        if( downPerturbedTrueAnomaly != downPerturbedTrueAnomaly )
        {
            std::cout<<"Error 2 in partial of true anomaly wrt cartesian state, element"<<i<<
                       ", keplerian state: "<<std::endl<<
                       convertCartesianToKeplerianElements( perturbedCartesianElements, gravitationalParameter ).transpose( )<<std::endl<<
                       perturbedCartesianElements.transpose( )<<std::endl;
        }

        // Calculate central difference of term i
        partial( 0, i ) = ( upPerturbedTrueAnomaly - downPerturbedTrueAnomaly ) / ( 2.0 * perturbations( i ) );

        // Check validity of result
        if( partial( 0, i ) != partial( 0, i ) )
        {
            std::cout<<"Error 2 in partial of true anomaly wrt cartesian state, element"<<i<<std::endl;
        }
    }

    return partial;
}

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > EmpiricalAccelerationPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case empirical_acceleration_coefficients:
        {
            partialFunction = boost::bind(
                        &EmpiricalAccelerationPartial::wrtEmpiricalAccelerationCoefficientFromIndices, this, parameter->getParameterSize( ),
                        boost::dynamic_pointer_cast< EmpiricalAccelerationCoefficientsParameter >( parameter )->getIndices( ), _1 );
            numberOfRows = parameter->getParameterSize( );
            break;
        }
        case arc_wise_empirical_acceleration_coefficients:
        {
            partialFunction = boost::bind(
                        &EmpiricalAccelerationPartial::wrtArcWiseEmpiricalAccelerationCoefficient, this,
                        boost::dynamic_pointer_cast< ArcWiseEmpiricalAccelerationCoefficientsParameter >( parameter ), _1 );
            numberOfRows = parameter->getParameterSize( );
            break;
        }
        default:
            break;
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}

void EmpiricalAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {

        using namespace tudat::basic_mathematics;
        using namespace tudat::linear_algebra;

        Eigen::Vector6d currentState = empiricalAcceleration_->getCurrentState( );
        Eigen::Vector3d angularMomentumVector = Eigen::Vector3d( currentState.segment( 0, 3 ) ).cross(
                    Eigen::Vector3d( currentState.segment( 3, 3 ) ) );
        Eigen::Vector3d crossVector = angularMomentumVector.cross( Eigen::Vector3d( currentState.segment( 0, 3 ) ) );

        Eigen::Matrix3d normPositionWrtPosition = calculatePartialOfNormalizedVector( Eigen::Matrix3d::Identity( ), currentState.segment( 0, 3 ) );
        Eigen::Matrix3d angularMomentumWrtPosition = -getCrossProductMatrix( currentState.segment( 3, 3 ) );
        Eigen::Matrix3d angularMomentumWrtVelocity = getCrossProductMatrix( currentState.segment( 0, 3 ) );
        Eigen::Matrix3d normAngularMomentumWrtPosition = calculatePartialOfNormalizedVector( angularMomentumWrtPosition, angularMomentumVector );
        Eigen::Matrix3d normAngularMomentumWrtVelocity = calculatePartialOfNormalizedVector( angularMomentumWrtVelocity, angularMomentumVector );
        Eigen::Matrix3d crossVectorWrtPosition = getCrossProductMatrix( angularMomentumVector ) -
                angularMomentumWrtVelocity * angularMomentumWrtPosition;
        Eigen::Matrix3d crossVectorWrtVelocity = -getCrossProductMatrix( currentState.segment( 0, 3 ) ) *
                getCrossProductMatrix( currentState.segment( 0, 3 ) );
        Eigen::Matrix3d normCrossVectorWrtPosition = calculatePartialOfNormalizedVector( crossVectorWrtPosition, crossVector );
        Eigen::Matrix3d normCrossVectorWrtVelocity = calculatePartialOfNormalizedVector( crossVectorWrtVelocity, crossVector );

        Eigen::Vector3d localAcceleration = empiricalAcceleration_->getCurrentLocalAcceleration( );

        currentPositionPartial_ = localAcceleration.x( ) * normPositionWrtPosition + localAcceleration.y( ) * normCrossVectorWrtPosition +
                localAcceleration.z( ) * normAngularMomentumWrtPosition;
        currentVelocityPartial_ = localAcceleration.y( ) * normCrossVectorWrtVelocity + localAcceleration.z( ) * normAngularMomentumWrtVelocity;

        Eigen::Matrix< double, 1, 6 > localTrueAnomalyPartial = calculateNumericalPartialOfTrueAnomalyWrtState(
                    empiricalAcceleration_->getCurrentState( ),
                    empiricalAcceleration_->getCurrentGravitationalParameter( ), 0.1 * perturbations );

        Eigen::Matrix< double, 1, 6 > trueAnomalyPartial = localTrueAnomalyPartial;

        currentPositionPartial_ += empiricalAcceleration_->getCurrentToInertialFrame( ) * (
                    empiricalAcceleration_->getAccelerationComponent( basic_astrodynamics::sine_empirical ) * std::cos( empiricalAcceleration_->getCurrentTrueAnomaly( ) ) -
                    empiricalAcceleration_->getAccelerationComponent( basic_astrodynamics::cosine_empirical ) * std::sin( empiricalAcceleration_->getCurrentTrueAnomaly( ) )
                    ) * trueAnomalyPartial.block( 0, 0, 1, 3 );

        currentVelocityPartial_ += ( empiricalAcceleration_->getCurrentToInertialFrame( ) * (
                                         empiricalAcceleration_->getAccelerationComponent( basic_astrodynamics::sine_empirical ) * std::cos( empiricalAcceleration_->getCurrentTrueAnomaly( ) ) -
                                         empiricalAcceleration_->getAccelerationComponent( basic_astrodynamics::cosine_empirical ) * std::sin( empiricalAcceleration_->getCurrentTrueAnomaly( ) )
                                         ) ) * trueAnomalyPartial.block( 0, 3, 1, 3 );
        currentTime_ = currentTime;

        if( currentPositionPartial_ != currentPositionPartial_ )
        {
            std::cerr<<"Error 2: "<<std::endl;
            std::cout<<trueAnomalyPartial<<std::endl;
            std::cout<<localTrueAnomalyPartial<<std::endl;
        }

        if( currentVelocityPartial_ != currentVelocityPartial_ )
        {
            std::cerr<<"Error 3: "<<std::endl;
            std::cout<<trueAnomalyPartial<<std::endl;
            std::cout<<localTrueAnomalyPartial<<std::endl;

        }
    }
}

}

}
