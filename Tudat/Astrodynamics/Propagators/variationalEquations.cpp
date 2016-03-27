#include <map>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Propagators/variationalEquations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"


namespace tudat
{

namespace propagators
{

using orbit_determination::partial_derivatives::StateDerivativePartialsMap;
using namespace tudat::estimatable_parameters;

//std::map< IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap > makeStateDerivativePartialList(
//        const orbit_determination::partial_derivatives::AccelerationPartialsMap& accelerationPartials )
//{
//    //std::vector< std::map< std::string, std::vector< boost::shared_ptr< partial_derivatives::AccelerationPartial > > > >
//    orbit_determination::partial_derivatives::StateDerivativePartialsMap stateDerivativePartials;
//    for( unsigned int i = 0; i < accelerationPartials.size( ); i++ )
//    {
//        std::map< std::string, std::vector< boost::shared_ptr< orbit_determination::partial_derivatives::StateDerivativePartial > > >
//                currentStateDerivativePartialList;
//        for( std::map< std::string, std::vector< boost::shared_ptr< orbit_determination::partial_derivatives::AccelerationPartial > > >::const_iterator
//             partialIterator = accelerationPartials.at( i ).begin( ); partialIterator != accelerationPartials.at( i ).end( );
//             partialIterator++ )
//        {
//            std::vector< boost::shared_ptr< orbit_determination::partial_derivatives::StateDerivativePartial > > currentStateDerivativePartials;
//            for( unsigned int j = 0; j < partialIterator->second.size( ); j++ )
//            {
//                currentStateDerivativePartials.push_back( partialIterator->second.at( j ) );
//            }
//            currentStateDerivativePartialList[ partialIterator->first ] = currentStateDerivativePartials;
//        }
//        stateDerivativePartials.push_back( currentStateDerivativePartialList );
//    }
//    std::map< IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap > stateDerivativePartialList;
//    stateDerivativePartialList[ transational_state ] = stateDerivativePartials;
//    return stateDerivativePartialList;

//}


//! Calculates matrix containing partial derivatives of accelerarion w.r.t. body state.
Eigen::MatrixXd& VariationalEquations::getBodyInitialStatePartialMatrix( )
{
    using namespace orbit_determination::partial_derivatives;

    // Initialize partial matrix
    variationalMatrix_.setZero( );

    if( dynamicalStatesToEstimate_.count( propagators::transational_state ) > 0 )
    {
        int startIndex = stateTypeStartIndices_.at( propagators::transational_state );
        for( unsigned int i = 0; i < dynamicalStatesToEstimate_.at( propagators::transational_state ).size( ); i++ )
        {
            variationalMatrix_.block( startIndex + i * 6, startIndex + i * 6 + 3, 3, 3 ).setIdentity( );
        }
    }

    // Iterate over all bodies undergoing accelerations for which initial condition is to be estimated.
    for( std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >, boost::function< Eigen::MatrixXd( ) > > > >::iterator
         typeIterator = statePartialList_.begin( ); typeIterator != statePartialList_.end( ); typeIterator++ )
    {
        int startIndex = stateTypeStartIndices_.at( typeIterator->first );
        int currentStateSize = getSingleIntegrationSize( typeIterator->first );
        int entriesToSkipPerEntry = currentStateSize - currentStateSize / getSingleIntegrationDifferentialEquationOrder( typeIterator->first );
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            // Iterate over all bodies exerting an acceleration on this body.
            for( statePartialIterator_ = typeIterator->second.at( i ).begin( ); statePartialIterator_ != typeIterator->second.at( i ).end( );
                 statePartialIterator_++ )
            {
                variationalMatrix_.block(
                            startIndex + entriesToSkipPerEntry + i * currentStateSize, statePartialIterator_->first.first,
                            currentStateSize - entriesToSkipPerEntry, statePartialIterator_->first.second ) += statePartialIterator_->second( );

            }
        }
    }

   for( unsigned int i = 0; i < statePartialAdditionIndices_.size( ); i++ )
   {
       //std::cout<<statePartialAdditionIndices_.at( i ).first<<" "<<statePartialAdditionIndices_.at( i ).second<<std::endl;
       variationalMatrix_.block( 0, statePartialAdditionIndices_.at( i ).second, totalDynamicalStateSize_, 3 ) +=
               variationalMatrix_.block( 0, statePartialAdditionIndices_.at( i ).first, totalDynamicalStateSize_, 3 );
   }

    //variationalMatrix_.block( 0, 6, 12, 3 ) += variationalMatrix_.block( 0, 0, 12, 3 );

    return variationalMatrix_;
}

//! Calculates matrix containing partial derivatives of accelerarion w.r.t. parameters.
Eigen::MatrixXd& VariationalEquations::getParameterPartialMatrix( const double ephemerisTime )
{
    // Initialize matrix to zeros
    variationalParameterMatrix_.setZero( );

    if( estimatedUnintegratedBodiesVectorSize_ > 0 )
    {
        std::cerr<<"Error, unintegrated body partials disabled 2"<<std::endl;
        //        for( unsigned int i = 0; i < unintegratedBodyPartialList_.size( ); i++ )
        //        {
        //            // Iterate over all parameter partial functions determined by setParameterPartialFunctionList( )
        //            for( stateFunctionIterator = unintegratedBodyPartialList_[ i ].begin( );
        //                 stateFunctionIterator != unintegratedBodyPartialList_[ i ].end( );
        //                 stateFunctionIterator++ )
        //            {
        //                // Add parameter partial to matrix.
        //                variationalMatrix.block( 3 + 6 * i, stateFunctionIterator->first, 3, 3 ) = stateFunctionIterator->second( );
        //            }
        //        }
        //        variationalMatrix.block( 0, totalDynamicalStateSize_, totalDynamicalStateSize_, estimatedUnintegratedBodiesVectorSize_ ) *=
        //                estimatedUnintegratedBodiesStateTransitionMatrixFunction_( ephemerisTime );
    }

    // Iterate over all bodies undergoing accelerations for which initial condition is to be estimated.
    for( std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >, boost::function< Eigen::MatrixXd& ( ) > > > >::iterator typeIterator =
         parameterPartialList_.begin( ); typeIterator != parameterPartialList_.end( ); typeIterator++ )
    {
        int startIndex = stateTypeStartIndices_.at( typeIterator->first );
        int currentStateSize = getSingleIntegrationSize( typeIterator->first );
        int entriesToSkipPerEntry = currentStateSize - currentStateSize / getSingleIntegrationDifferentialEquationOrder( typeIterator->first );

        // Iterate over all bodies being estimated.
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            // Iterate over all parameter partial functions determined by setParameterPartialFunctionList( )
            for( functionIterator = typeIterator->second[ i ].begin( ); functionIterator != typeIterator->second[ i ].end( );
                 functionIterator++ )
            {
                // Add parameter partial to matrix.
                variationalParameterMatrix_.block(
                            startIndex + entriesToSkipPerEntry + currentStateSize * i,
                            functionIterator->first.first,
                            currentStateSize - entriesToSkipPerEntry, functionIterator->first.second ) += functionIterator->second( );
            }
        }

    }

    return variationalParameterMatrix_;
}


//! Evaluates the variational equations.
template< >
Eigen::MatrixXd& VariationalEquations::evaluateVariationalEquations< double >(
        const double ephemerisTime,
        const Eigen::MatrixXd& stateTransitionAndSensitivityMatrices )
{
    // Initialize return matrix to zero.
    currentMatrixDerivative_.setZero( );

    // Add partials of body positions and velocities.
    currentMatrixDerivative_.block( 0, 0, totalDynamicalStateSize_, numberOfParameterValues_ ) +=
            getBodyInitialStatePartialMatrix( ) * stateTransitionAndSensitivityMatrices;

    // Add partials of parameters.
    currentMatrixDerivative_.block( 0, 0, totalDynamicalStateSize_, numberOfParameterValues_ ) +=
            getParameterPartialMatrix( ephemerisTime );

    return currentMatrixDerivative_;
}

template< >
Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& VariationalEquations::evaluateVariationalEquations< long double >(
        const double ephemerisTime,
        const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& stateTransitionAndSensitivityMatrices )
{
    currentLongMatrixDerivative_  =
            evaluateVariationalEquations< double >( ephemerisTime, stateTransitionAndSensitivityMatrices.cast< double >( ) ).
                cast< long double >( );

   //std::cout<<"state trans: "<<std::endl<<stateTransitionAndSensitivityMatrices<<std::endl<<std::endl;
    return currentLongMatrixDerivative_;
}

//! This function updates the total state of each body, acceleration and acceleration partial in the simulation at the given time and state
//! of bodies that are integrated numerically.
void VariationalEquations::updatePartials( const double currentTime )
{
    for( std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator
         stateDerivativeTypeIterator = stateDerivativePartialList_.begin( ); stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator++ )
    {
        for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
        {
            for( unsigned int j = 0; j < stateDerivativeTypeIterator->second.at( i ).size( ); j++ )
            {
                stateDerivativeTypeIterator->second.at( i ).at( j )->resetTime( TUDAT_NAN );
            }

        }
    }

    // Update all acceleration partials to current state and time. Information is passed indirectly from here, through
    // (function) pointers set in acceleration partial classes
    for( std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator
         stateDerivativeTypeIterator = stateDerivativePartialList_.begin( ); stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator++ )
    {
        for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
        {
            for( unsigned int j = 0; j < stateDerivativeTypeIterator->second.at( i ).size( ); j++ )
            {
                stateDerivativeTypeIterator->second.at( i ).at( j )->update( currentTime );
            }

        }
    }

//    std::cout<<"Updating partials"<<std::endl;
    for( std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator
         stateDerivativeTypeIterator = stateDerivativePartialList_.begin( ); stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator++ )
    {
        for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
        {
            for( unsigned int j = 0; j < stateDerivativeTypeIterator->second.at( i ).size( ); j++ )
            {
                stateDerivativeTypeIterator->second.at( i ).at( j )->updateParameterPartials( );
            }

        }
    }
}

void VariationalEquations::setStatePartialFunctionList( )
{
    using namespace orbit_determination::partial_derivatives;

    std::pair< boost::function< Eigen::MatrixXd( ) >, int > currentDerivativeFunction;

    for( std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator
         stateDerivativeTypeIterator = stateDerivativePartialList_.begin( ); stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator++ )
    {
        // Iterate over all bodies undergoing accelerations for which initial condition is to be estimated.
        for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
        {
            std::multimap< std::pair< int, int >, boost::function< Eigen::MatrixXd( ) > > currentBodyPartialList;
            // Iterate over all accelerations from single body on other single body
            for( unsigned int j = 0; j < stateDerivativeTypeIterator->second.at( i ).size( ); j++ )
            {
                for( std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >::iterator
                     estimatedStateIterator = dynamicalStatesToEstimate_.begin( ); estimatedStateIterator != dynamicalStatesToEstimate_.end( );
                     estimatedStateIterator++ )
                {
                    // Iterate over all bodies to see if body exerting acceleration is also to be estimated (cross-terms)
                    for( unsigned int k = 0; k < estimatedStateIterator->second.size( ); k++ )
                    {
                        currentDerivativeFunction = stateDerivativeTypeIterator->second.at( i ).at( j )->
                                getDerivativeFunctionWrtStateOfIntegratedBody(
                                    estimatedStateIterator->second.at( k ), estimatedStateIterator->first );

                        if( currentDerivativeFunction.second != 0 )
                        {
                            currentBodyPartialList.insert(
                                        std::make_pair(
                                            std::make_pair( k * getSingleIntegrationSize( estimatedStateIterator->first ) +
                                                            stateTypeStartIndices_.at( estimatedStateIterator->first ),
                                                            getSingleIntegrationSize( estimatedStateIterator->first ) ),
                                            currentDerivativeFunction.first ) );
                        }
                    }
                }
            }
            statePartialList_[ stateDerivativeTypeIterator->first ].push_back( currentBodyPartialList );
        }
    }
}


void VariationalEquations::setUnintegratedBodyPartialList( )
{
    using namespace tudat::orbit_determination::partial_derivatives;

    for( std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator
         stateDerivativeTypeIterator = stateDerivativePartialList_.begin( ); stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
         stateDerivativeTypeIterator++ )
    {
        // Iterate over all bodies undergoing accelerations for which initial condition is to be estimated.
        for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
        {
            std::multimap< int, boost::function< Eigen::Matrix3d( ) > > functionListOfBody;

            // Iterate over all bodies to see if body exerting acceleration is also to be estimated (cross-terms)
            for( unsigned int k = 0; k < estimatedUnintegratedBodies_.size( ); k++ )
            {
                // Iterate over all bodies exerting an acceleration on this body.
                for( unsigned int j = 0; j < stateDerivativeTypeIterator->second.at( i ).size( ); j++ )
                {
                    if( stateDerivativeTypeIterator->second.at( i ).at( j )->getIntegrationReferencePoint( ).first ==
                            estimatedUnintegratedBodies_[ k ] )
                    {

                        std::cerr<<"Error, unintegrated body partials disabled"<<std::endl;
                        //                        functionListOfBody.insert(
                        //                                    std::pair< int, boost::function< Eigen::Matrix3d( ) > >
                        //                                    ( totalDynamicalStateSize_ + k * 6,
                        //                                      boost::bind( &AccelerationPartial::wrtPositionOfAcceleratingBody,
                        //                                                   accelerationIterator->second[ j ] ) ) );
                        //                        functionListOfBody.insert(
                        //                                    std::pair< int, boost::function< Eigen::Matrix3d( ) > >
                        //                                    ( totalDynamicalStateSize_ + k * 6 + 3,
                        //                                      boost::bind( &AccelerationPartial::wrtVelocityOfAcceleratingBody,
                        //                                                   accelerationIterator->second[ j ] ) ) );
                    }
                }
            }
        }
        //unintegratedBodyPartialList_.push_back( functionListOfBody );
    }
}


}

}
