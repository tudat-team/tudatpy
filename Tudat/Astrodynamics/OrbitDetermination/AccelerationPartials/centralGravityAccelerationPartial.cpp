#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"


namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Calculates the partial derivative of a point mass gravitational acceleration wrt the position of the body being accelerated.
Eigen::Matrix3d calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
        const Eigen::Vector3d stateOfGravitatedBody,
        const Eigen::Vector3d stateOfGravitatingBody,
        double gravitationalParameter )
{
    // Calculate relative position
    Eigen::Vector3d relativePosition = stateOfGravitatedBody - stateOfGravitatingBody;

    // Calculate partial (Montenbruck & Gill, Eq. 7.56)
    double relativePositionNorm = relativePosition.norm( );
    double invSquareOfPositionNorm = 1.0 / ( relativePositionNorm * relativePositionNorm );
    double invCubeOfPositionNorm = invSquareOfPositionNorm / relativePositionNorm;
    Eigen::Matrix3d partialMatrix = -1.0 * gravitationalParameter *
            ( Eigen::Matrix3d::Identity( ) * invCubeOfPositionNorm -
              ( 3.0 * invSquareOfPositionNorm * invCubeOfPositionNorm ) * relativePosition * relativePosition.transpose( ) );

    return partialMatrix;
}

//! Calculates the partial derivative of a point mass gravitational acceleration wrt the gravitational parameter of the central body
Eigen::Vector3d computePartialOfCentralGravityWrtGravitationalParameter( const Eigen::Vector3d acceleratedBodyPosition,
                                                                         const Eigen::Vector3d acceleratingBodyPosition)
{
    // Calculate relative position
    Eigen::Vector3d relativePosition = acceleratingBodyPosition - acceleratedBodyPosition;

    // Calculate partial (Montenbruck & Gill, Eq. 7.76)
    double positionNorm = relativePosition.norm( );
    Eigen::Vector3d partialMatrix = relativePosition / ( positionNorm * positionNorm * positionNorm );
    return partialMatrix;
}

//! Constructor
CentralGravitationPartial::CentralGravitationPartial(
        const boost::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > gravitationalAcceleration,
        const std::string acceleratedBody,
        const std::string acceleratingBody ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::central_gravity )
{
    gravitationalParameterFunction_ = gravitationalAcceleration->getGravitationalParameterFunction( );
    centralBodyState_ = gravitationalAcceleration->getStateFunctionOfBodyExertingAcceleration( );
    acceleratedBodyState_ = gravitationalAcceleration->getStateFunctionOfBodyUndergoingAcceleration( );
    accelerationUsesMutualAttraction_ = gravitationalAcceleration->getIsMutualAttractionUsed( );
}

CentralGravitationPartial::CentralGravitationPartial(
        const boost::shared_ptr< CentralGravitationPartial > originalAccelerationPartial ):
    AccelerationPartial( originalAccelerationPartial->getAcceleratedBody( ),
                         originalAccelerationPartial->getAcceleratingBody( ), basic_astrodynamics::central_gravity )
{
    gravitationalParameterFunction_ = originalAccelerationPartial->getGravitationalParameterFunction( );
    centralBodyState_ = originalAccelerationPartial->getPositionFunctionOfBodyExertingAcceleration( );
    acceleratedBodyState_ = originalAccelerationPartial->getPositionFunctionOfBodyUndergoingAcceleration( );
    accelerationUsesMutualAttraction_ = originalAccelerationPartial->getAccelerationUsesMutualAttraction( );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< boost::function< Eigen::MatrixXd( ) >, int >
CentralGravitationPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

{
    std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionPair;

    if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter )
    {
        partialFunctionPair = getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
    }
    else
    {
        partialFunctionPair = std::make_pair( boost::function< Eigen::MatrixXd( ) >( ), 0 );
    }

    return partialFunctionPair;
}

std::pair< boost::function< Eigen::MatrixXd( ) >, int >
CentralGravitationPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    boost::function< Eigen::MatrixXd( ) > partialFunction;
    int numberOfColumns = 0;

    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check for dependency
        if( parameterId.second.first == acceleratingBody_ )
        {
            partialFunction = boost::bind( &CentralGravitationPartial::wrtGravitationalParameterOfCentralBody,
                                           this );
            numberOfColumns = 1;

        }

        if( parameterId.second.first == acceleratedBody_ )
        {

            if( accelerationUsesMutualAttraction_ )
            {
                partialFunction = boost::bind( &CentralGravitationPartial::wrtGravitationalParameterOfCentralBody,
                                               this );
                numberOfColumns = 1;
            }
        }
    }
    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
Eigen::Vector3d CentralGravitationPartial::wrtGravitationalParameterOfCentralBody( )
{
    return computePartialOfCentralGravityWrtGravitationalParameter(
                acceleratedBodyState_( ).segment( 0, 3 ), centralBodyState_( ).segment( 0, 3 ) );
}



}

}

}
