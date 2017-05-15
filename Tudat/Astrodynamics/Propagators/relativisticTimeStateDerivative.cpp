#include "Tudat/Astrodynamics/Propagators/relativisticTimeStateDerivative.h"

namespace tudat
{

namespace propagators
{

PostNewtonianRelativisticTimeStateDerivative::PostNewtonianRelativisticTimeStateDerivative(
        const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
        const std::pair< std::string, std::string >& referencePoint,
        const std::vector< std::string >& externalBodies,
        const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType,
        const boost::function< double( const double ) > timeVariableConversionFunction,
        const double distanceScalingFactor ):RelativisticTimeStateDerivative(
        referencePoint, relativisticStateDerivativeType ),
    timeVariableConversionFunction_( timeVariableConversionFunction ),distanceScalingFactor_( distanceScalingFactor )
{
    // Check if central body exists.
    if( referencePoint.first == "SSB" )
    {
        centralBodyStateFunction_ = boost::lambda::constant( basic_mathematics::Vector6d::Zero( ) );
    }
    else if( bodyMap.count( referencePoint.first ) == 0 )
    {
        std::cerr<<"Error, could not find centeal body in body map when making 1st order relativistic time converter"<<std::endl;
    }
    else
    {
        // Set central body state function.
        centralBodyStateFunction_ = boost::bind( &bodies::Body::getState, bodyMap.at( referencePoint.first ) );
    }

    // Iterate over all external bodies and set their state and gravitational parameter functions.
    for( unsigned int i = 0; i < externalBodies.size( ); i++ )
    {
        // Check if current body exists.
        if( bodyMap.count( externalBodies.at( i ) ) == 0 )
        {
            std::cerr<<"Error, could not find body "<<externalBodies.at( i )<<" in body map when making 1st order relativistic time converter"<<std::endl;
        }
        // Check if current body is a CelestialBody
        else if( boost::dynamic_pointer_cast< bodies::CelestialBody >( bodyMap.at( externalBodies.at( i ) ) ) == NULL )
        {
            std::cerr<<"Error, could not find celestial body "<<externalBodies.at( i )<<" in body map when making 1st order relativistic time converter"<<std::endl;
        }
        // Check if current body has a gravity field.
        else if( boost::dynamic_pointer_cast< bodies::CelestialBody >( bodyMap.at( externalBodies.at( i ) ) )->getGravityFieldModel( ) == NULL )
        {
            std::cerr<<"Error, celestial body "<<externalBodies.at( i )<<" has no gravity field model when making 1st order relativistic time converter"<<std::endl;
        }
        else
        {
            // Set state and gravitational parameter functions of current body.
            boost::shared_ptr< bodies::CelestialBody > celestialBody =
                    boost::dynamic_pointer_cast< bodies::CelestialBody >( bodyMap.at( externalBodies.at( i ) ) );
            externalBodyStateFunctions_.push_back( boost::bind( &bodies::CelestialBody::getTemplatedState< double >, celestialBody ) );
            externalBodyGravitationalParameterFunctions_.push_back( boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                                                 celestialBody->getGravityFieldModel( ) ) );
        }
    }
}


//! Class contructor, sets required properties of bodies.
FirstOrderBarycentricToBodyCentricTimeStateDerivative::FirstOrderBarycentricToBodyCentricTimeStateDerivative(
        const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap, const std::string& centralBody,
        const std::vector< std::string >& externalBodies,
        const std::map< std::string, std::pair< int, int > > sphericalHarmonicGravityExpansions,
        const boost::function< double( const double ) > timeVariableConversionFunction,
        const double distanceScalingFactor,
        const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType ):
    PostNewtonianRelativisticTimeStateDerivative( bodyMap, std::make_pair( centralBody, "" ), externalBodies, relativisticStateDerivativeType,
                                                  timeVariableConversionFunction, distanceScalingFactor ),
    externalBodies_( externalBodies )
{
    // Set sizes of vectors which are to hold the current states of the system.
    currentExternalBodyDistances_.resize( externalBodyStateFunctions_.size( ) );
    currentExternalBodyGravitationalParameters_.resize( externalBodyGravitationalParameterFunctions_.size( ) );
    currentExternalBodyStates_.resize( externalBodyStateFunctions_.size( ) );

    for( unsigned int i = 0; i < externalBodies.size( ); i++ )
    {
        if( sphericalHarmonicGravityExpansions.count( externalBodies.at( i ) ) > 0 )
        {
            boost::shared_ptr< bodies::CelestialBody > currentCelestialBody = boost::dynamic_pointer_cast< bodies::CelestialBody >(
                        bodyMap.at( externalBodies.at( i ) ) );
            if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( currentCelestialBody->getGravityFieldModel( ) ) == NULL )
            {
                std::cerr<<"Error when making sh correction for relativistic time conversion, body "<<externalBodies.at( i )<<" has no"<<
                           " sh gravity field."<<std::endl;
            }
            else
            {
                basic_mathematics::LegendreCache* legendreCache = new basic_mathematics::LegendreCache(
                            sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).first,
                            sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).second );
                boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > shGravityField =
                        boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( currentCelestialBody->getGravityFieldModel( ) );
                higherOrderGravityFieldPotentialFunctions_[ i ] = boost::bind(
                            &gravitation::SphericalHarmonicsGravityField::getGravitationalPotential, shGravityField, _1,
                            sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).first,
                            sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).second, legendreCache, 2, 0 );
            }
        }
    }
}

//! Base class function for evaluating the integrand which is required for the conversion.
Eigen::Matrix< double, Eigen::Dynamic, 1 > FirstOrderBarycentricToBodyCentricTimeStateDerivative::calculateSystemStateDerivative(
        const double currentGlobalTime, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& integratedValue )
{
    // Update bodies to current time.
    //updateStateDerivativeModel( currentGlobalTime );

    // Calculate and return value of integrand
    return ( Eigen::Matrix< double, 1, 1 >( ) << calculateFirstOrderTcbToTcgIntegrand(
                 currentVelocity_, currentExternalPotential_ ) ).finished( );
}

//! Function to update all environment variables to current time.
void FirstOrderBarycentricToBodyCentricTimeStateDerivative::updateStateDerivativeModel( double baseFrameTime )
{
    // Update variables relevant for both this class and its derived class (2nd order conversion)
    updateBaseVariables( baseFrameTime );

    // Calculate state, and distance to central body, of each external body.
    for( unsigned int i = 0; i < externalBodyStateFunctions_.size( ); i++ )
    {
        currentExternalBodyStates_[ i ] = externalBodyStateFunctions_[ i ]( );
        currentExternalBodyStates_[ i ].segment( 0, 3 ) = currentExternalBodyStates_[ i ].segment( 0, 3 ) * distanceScalingFactor_;

        currentExternalBodyDistances_[ i ] = ( currentCentralBodyState_.segment( 0, 3 ) - currentExternalBodyStates_[ i ].segment( 0, 3 ) ).norm( );
    }



    // Calculate 1st order external potential.
    if( higherOrderGravityFieldPotentialFunctions_.size( ) == 0 )
    {
        currentExternalPotential_ = calculateFirstOrderExternalScalarPotential(
                    currentExternalBodyGravitationalParameters_, currentExternalBodyDistances_ );
    }
    else
    {
        for( shPotentialIterator_ = higherOrderGravityFieldPotentialFunctions_.begin( ); shPotentialIterator_ !=
             higherOrderGravityFieldPotentialFunctions_.end( ); shPotentialIterator_++ )
        {
            std::cerr<<"Implement this potential stuff"<<std::endl;
        }
    }


}

//! Function to update all environment variables common to this class, and its (2nd order) derived class to current time.
void FirstOrderBarycentricToBodyCentricTimeStateDerivative::updateBaseVariables( double baseFrameTime )
{
    double convertedTime = timeVariableConversionFunction_( baseFrameTime );

    std::map< propagators::IntegratedStateType, Eigen::VectorXd > currentState;
    currentState[ propagators::proper_time ] = ( Eigen::VectorXd( 1 ) << 0.0 ).finished( );
    //environmentUpdater_->updateEnvironment( convertedTime, currentState );

    // Calculate state of central body.
    currentCentralBodyState_ = centralBodyStateFunction_( );
    currentCentralBodyState_.segment( 0, 3 ) = currentCentralBodyState_.segment( 0, 3 ) * distanceScalingFactor_;

    currentVelocity_ = currentCentralBodyState_.segment( 3, 3 ).norm( );
    // Retrieve gravitational parameters of all bodies.
    for( unsigned int i = 0; i < externalBodyStateFunctions_.size( ); i++ )
    {
        currentExternalBodyGravitationalParameters_[ i ] = externalBodyGravitationalParameterFunctions_[ i ]( ) * distanceScalingFactor_;
    }
}

//! Constructor of bodycentric<->topocentric time calculator for a single ground station on a celestial body.
FirstOrderBodyCentricToTopoCentricTimeCalculator::FirstOrderBodyCentricToTopoCentricTimeCalculator(
        const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
        const std::string& centralBody,
        const std::vector< std::string >& externalBodies,
        const std::string& groundStation,
        const int maximumSphericalHarmonicDegree,
        const bool useAccelerationTerm,
        const bool useTimeDependentBodyFixedPosition ): PostNewtonianRelativisticTimeStateDerivative(
        bodyMap, std::make_pair( centralBody, groundStation ), externalBodies, propagators::first_order_bodycentric_to_topocentric )
{
    // Check if central body has a rotation to body-fixed frame.
    if( bodyMap.at( centralBody )->getRotationalEphemeris( ) == NULL )
    {
        std::cerr<<"Error, could not find rotation model of central body when making PCRS to TPRS time transformation."<<std::endl;
    }
    else
    {
        // Get rotation vector and rotation matrix functions.
        centralBodyRotationVector_ = boost::bind( &ephemerides::RotationalEphemeris::getRotationalVelocityVectorInBaseFrame,
                                                  bodyMap.at( centralBody )->getRotationalEphemeris( ), _1 );
        toInertialFrameTransformation_ = boost::bind( &ephemerides::RotationalEphemeris::getRotationToBaseFrame,
                                                      bodyMap.at( centralBody )->getRotationalEphemeris( ), _1 );
    }

    // Check if central body is a celestial body
    if( boost::dynamic_pointer_cast< bodies::CelestialBody >( bodyMap.at( centralBody ) ) == NULL )
    {
        std::cerr<<"Error, central body "<<centralBody<<" is not a celestial body when making PCRS to TPRS time transformation."<<std::endl;
    }
    else
    {
        boost::shared_ptr< bodies::CelestialBody > centralCelestialBody =
                boost::dynamic_pointer_cast< bodies::CelestialBody >( bodyMap.at( centralBody ) );

        if( centralCelestialBody->getGroundStationMap( ).count( groundStation ) == 0 )
        {
            std::cerr<<"Error, station "<<groundStation<<" not found on "<<centralBody<<" when making PCRS to TPRS time transformation."<<std::endl;
        }
        else
        {
            // Get ground station position function.
            boost::shared_ptr< NominalGroundStationState > groundStationState = centralCelestialBody->getGroundStation(
                        groundStation )->getNominalStationState( );
            if( !useTimeDependentBodyFixedPosition )
            {
                pointPositionFunctionInPcrs_ = boost::bind( &NominalGroundStationState::getNominalCartesianPosition, groundStationState );
            }
            else
            {
                pointPositionFunctionInPcrs_ = boost::bind( &NominalGroundStationState::getCartesianPositionInTime, groundStationState, _1,
                                                            basic_astrodynamics::JULIAN_DAY_ON_J2000 );
            }

            // Get central body potential function.
            if( centralCelestialBody->getGravityFieldModel( ) == NULL )
            {
                std::cerr<<"Error, could not find gravity field model of central body when making PCRS to TPRS time transformation."<<std::endl;
            }
            else
            {
                // If spherical harmonic expansion is requested, check if required field exists.
                if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( centralCelestialBody->getGravityFieldModel( ) ) == NULL
                        && maximumSphericalHarmonicDegree > 0 )
                {
                    std::cerr<<"Error, requested spherical harmonic terms when making PCRS to TPRS time transformation, but no such model found."<<std::endl;
                }
                // If no spherical harmonic expansion is request, check if field is purely central, and set function.
                else if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( centralCelestialBody->getGravityFieldModel( ) ) ==
                         NULL && maximumSphericalHarmonicDegree == 0 )
                {
                    localCentralBodyPotentialFunction_ = boost::bind( &gravitation::GravityFieldModel::getGravitationalPotential,
                                                                      centralCelestialBody->getGravityFieldModel( ), _1 );
                }
                else
                {
                    // Retrieve spherical harmonic potential function from gravtiy field object.
                    boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicField =
                            boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( centralCelestialBody->getGravityFieldModel( ) );
                    localCentralBodyPotentialFunction_ = boost::bind(
                                &gravitation::SphericalHarmonicsGravityField::getGravitationalPotential,
                                sphericalHarmonicField, _1, maximumSphericalHarmonicDegree, maximumSphericalHarmonicDegree,
                                new basic_mathematics::LegendreCache( maximumSphericalHarmonicDegree, maximumSphericalHarmonicDegree ), 0, 0 );
                }

            }
        }
    }

    // Add acceleration function, if requested.
    if( useAccelerationTerm )
    {
        std::cerr<<"Error, acceleration term not yet implemented for PCRS to TPRS time transformation."<<std::endl;
    }
    else
    {
        centralBodyAccelerationFunction_ = boost::lambda::constant( Eigen::Vector3d::Zero( ) );
    }

    currentExternalBodyGravitationalParameters_.resize( externalBodyStateFunctions_.size( ) );
    currentExternalBodyRelativePositions_.resize( externalBodyStateFunctions_.size( ) );
    currentExternalBodyDistances_.resize( externalBodyStateFunctions_.size( ) );
}

//! Function for evaluating the integrand which is required for the conversion.
Eigen::Matrix< double, Eigen::Dynamic, 1 > FirstOrderBodyCentricToTopoCentricTimeCalculator::calculateSystemStateDerivative(
        const double currentGlobalTime, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& integratedValue )
{
    //updateStateDerivativeModel( currentGlobalTime );
    return ( Eigen::Matrix< double, 1, 1 >( ) << calculateFirstOrderPlanetocentricToTopocentricConversion(
                 currentPointPositionInPcrs_, currentPointVelocityInPcrs_, currentLocalPotential_,
                 currentBarycentricAccelerationOfCentralBody_, currentExternalBodyRelativePositions_,
                 currentExternalBodyGravitationalParameters_ ) ).finished( );
}

//! Function to update all environment variables to current time.
void FirstOrderBodyCentricToTopoCentricTimeCalculator::updateStateDerivativeModel( const double baseFrameTime )
{
    std::map< propagators::IntegratedStateType, Eigen::VectorXd > currentState;
    currentState[ propagators::proper_time ] = ( Eigen::VectorXd( 1 ) << 0.0 ).finished( );
    //environmentUpdater_->updateEnvironment( baseFrameTime, currentState );
    // Get central body state
    Eigen::Matrix< double, 6, 1 > currentCentralBodyState = centralBodyStateFunction_( );
    Eigen::Vector3d currentCentralBodyPosition = currentCentralBodyState.segment( 0, 3 );
    currentCentralBodyBarycentricVelocity_ = currentCentralBodyState.segment( 3, 3 );

    // Get gravitational parameter and positions of all external bodies.
    for( unsigned int i = 0; i < externalBodyStateFunctions_.size( ); i++ )
    {
        currentExternalBodyGravitationalParameters_[ i ] = externalBodyGravitationalParameterFunctions_[ i ]( );
        currentExternalBodyRelativePositions_[ i ] = externalBodyStateFunctions_[ i ]( ).segment( 0, 3 ) - currentCentralBodyPosition;
        currentExternalBodyDistances_[ i ] = currentExternalBodyRelativePositions_[ i ].norm( );
    }

    // Get current position of ground station in body-fixed frame.
    currentPointPositionInPcrs_ = pointPositionFunctionInPcrs_( baseFrameTime );

    // Evaluate gravitational potential of central body at current postion of ground station.
    currentLocalPotential_ = localCentralBodyPotentialFunction_( currentPointPositionInPcrs_ );

    // Convert ground station position to non-co-rotating frame.
    currentPointPositionInPcrs_ = toInertialFrameTransformation_( baseFrameTime ) * currentPointPositionInPcrs_;

    // Get current velocity of ground station in non-co-rotating frame.
    currentPointVelocityInPcrs_ = centralBodyRotationVector_( baseFrameTime ).cross( currentPointPositionInPcrs_ );

    // Calculate barycentric acceleration of central body.
    currentBarycentricAccelerationOfCentralBody_ = centralBodyAccelerationFunction_( baseFrameTime );

}

}

}

