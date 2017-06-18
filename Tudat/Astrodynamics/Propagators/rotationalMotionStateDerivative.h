#ifndef ROTATIONALMOTIONSTATEDERIVATIVE_H
#define ROTATIONALMOTIONSTATEDERIVATIVE_H


#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace propagators
{

Eigen::Vector3d evaluateRotationalEquationsOfMotion(
        const Eigen::Matrix3d& inertiaTensor, const Eigen::Vector3d& totalTorque,
        const Eigen::Vector3d& rotationVector, const Eigen::Matrix3d& inertiaTimeDerivative = Eigen::Matrix3d::Zero( ) );

Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& currentBodyFixedRotationRate );

Eigen::Vector4d calculateQuaternionDerivative(
        const Eigen::Vector4d& currentQuaternion, const Eigen::Vector3d& currentBodyFixedRotationRate );


template< typename StateScalarType = double, typename TimeType = double >
class RotationalMotionStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;

    RotationalMotionStateDerivative(
            const basic_astrodynamics::TorqueModelMap& torqueModelsPerBody,
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::vector< std::string >& bodiesToPropagate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::rotational_state ),
        torqueModelsPerBody_( torqueModelsPerBody ),
        bodyMap_( bodyMap ),
        bodiesToPropagate_( bodiesToPropagate )
    {
        // Check whether the bodies that are to e integrated exist in bodyMap
        for( unsigned int i = 0; i < bodiesToPropagate_.size( ); i++ )
        {
            if( bodyMap.count( bodiesToPropagate_[ i ] ) == 0 )
            {
                std::cerr<<"Warning when creating RotationalMotionStateDerivative, body "<<bodiesToPropagate_.at( i )<<std::endl;
                std::cerr<<" not present in provided ody map, cannot integrate body!"<<std::endl;
            }


            bodyInertiaTensorFunctions_.push_back( boost::bind( &simulation_setup::Body::getBodyInertiaTensor,
                                                                bodyMap.at( bodiesToPropagate_.at( i ) ) ) );
            bodyInertiaTensorTimeDerivativeFunctions_.push_back( boost::lambda::constant( Eigen::Matrix3d::Zero( ) ) );
        }
    }

    ~RotationalMotionStateDerivative( ){ }

    void updateStateDerivativeModel( const TimeType currentTime )
    {
        for( torqueModelMapIterator = torqueModelsPerBody_.begin( );
             torqueModelMapIterator != torqueModelsPerBody_.end( ); torqueModelMapIterator++ )
        {
            for( innerTorqueIterator  = torqueModelMapIterator->second.begin( ); innerTorqueIterator !=
                 torqueModelMapIterator->second.end( ); innerTorqueIterator++ )
            {
                for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                {
                    innerTorqueIterator->second[ j ]->updateMembers( currentTime );
                }
            }
        }
    }


    virtual void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > inertialToBodyFixedRotationalState )
    {
        inertialToBodyFixedRotationalState = internalSolution;
        for( unsigned int i = 0; i < bodiesToPropagate_.size( ); i++ )
        {
            inertialToBodyFixedRotationalState.block( 7 * i, 0, 4, 1 ) =
                    internalSolution.block( 7 * i, 0, 4, 1 ).normalized( );
        }
    }

    void clearStateDerivativeModel( )
    {
        for( torqueModelMapIterator = torqueModelsPerBody_.begin( );
             torqueModelMapIterator != torqueModelsPerBody_.end( ); torqueModelMapIterator++ )
        {
            for( innerTorqueIterator  = torqueModelMapIterator->second.begin( ); innerTorqueIterator !=
                 torqueModelMapIterator->second.end( ); innerTorqueIterator++ )
            {
                for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                {
                    innerTorqueIterator->second[ j ]->resetTime( TUDAT_NAN );
                }
            }
        }
    }

    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( stateOfSystemToBeIntegrated.rows( ), 1 );
        std::vector< Eigen::Vector3d > torquesActingOnBodies = sumTorquesPerBody( );

        for( unsigned int i = 0; i < torquesActingOnBodies.size( ); i++ )
        {
            Eigen::Matrix< StateScalarType, 4, 1 > currentQuaternion = ( stateOfSystemToBeIntegrated.block( 7 * i, 0, 4, 1 ) ).normalized( );
            Eigen::Matrix< StateScalarType, 3, 1 > currentBodyFixedRotationRate = stateOfSystemToBeIntegrated.block( 7 * i + 4, 0, 3, 1 );

            stateDerivative.block( 7 * i, 0, 4, 1 ) = calculateQuaternionDerivative(
                        currentQuaternion.template cast< double >( ), currentBodyFixedRotationRate.template cast< double >( ) ).
                    template cast< StateScalarType >( );
            stateDerivative.block( 7 * i + 4, 0, 3, 1 ) = evaluateRotationalEquationsOfMotion(
                        bodyInertiaTensorFunctions_.at( i )( ), torquesActingOnBodies.at( i ),
                        currentBodyFixedRotationRate.template cast< double >( ),
                        bodyInertiaTensorTimeDerivativeFunctions_.at( i )( ) ).template cast< StateScalarType >( );

        }
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }

    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSoluton )
    {
        currentLocalSoluton = internalSolution;
    }


    std::string getNameOfPropagatedBody( int bodyIndex )
    {
        return bodiesToPropagate_[ bodyIndex ];
    }

    std::vector< std::string > getBodiesToBeIntegratedNumerically( )
    {
        return bodiesToPropagate_;
    }

    basic_astrodynamics::TorqueModelMap getAccelerationsMap( )
    {
        return torqueModelsPerBody_;
    }

    int getStateSize( )
    {
        return 7 * bodiesToPropagate_.size( );
    }

protected:

    std::vector< Eigen::Vector3d > sumTorquesPerBody( )
    {
        using namespace basic_astrodynamics;

        std::vector< Eigen::Vector3d > torques;
        torques.resize( bodiesToPropagate_.size( ) );

        for( unsigned int i = 0; i < bodiesToPropagate_.size( ); i++ )
        {
            torques[ i ].setZero( );

            if( torqueModelsPerBody_.count( bodiesToPropagate_[ i ] ) != 0 )
            {
                for( innerTorqueIterator  = torqueModelsPerBody_[ bodiesToPropagate_[ i ] ].begin( );
                     innerTorqueIterator != torqueModelsPerBody_[ bodiesToPropagate_[ i ] ].end( );
                     innerTorqueIterator++ )
                {
                    for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                    {
                        torques[ i ] += ( innerTorqueIterator->second[ j ]->getTorque( ) );
                    }
                }
            }
        }

        return torques;
    }

    std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions_;

    std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions_;

    basic_astrodynamics::TorqueModelMap torqueModelsPerBody_;

    simulation_setup::NamedBodyMap bodyMap_;

    std::vector< std::string > bodiesToPropagate_;

    basic_astrodynamics::TorqueModelMap::iterator torqueModelMapIterator;

    basic_astrodynamics::SingleBodyTorqueModelMap::iterator innerTorqueIterator;

};

}

}

#endif // ROTATIONALMOTIONSTATEDERIVATIVE_H
