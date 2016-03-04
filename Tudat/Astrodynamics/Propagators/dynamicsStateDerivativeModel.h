#ifndef HYBRIDSTATEDERIVATIVEMODEL_H
#define HYBRIDSTATEDERIVATIVEMODEL_H


#include <map>
#include <utility>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdater.h"

namespace tudat
{
namespace propagators
{



//! Class for the numerical integrations of general dynamical equations.
/*!
 *  Class for the numerical integrations of general dynamical equations.
 */
template< typename TimeType = double, typename StateScalarType = double >
class DynamicsStateDerivativeModel
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > StateType;
    typedef std::map< IntegratedStateType, std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >
    StateDerivativeCalculatorList;
    //! Derivative model constructor.
    /*!
     *  Derivative model constructor. Takes state derivative model and environment updater. Constructor checks whether all models use the same environment updater.
     *  \param stateDerivativeModels Vector of state derivative models, with one entry for each type of dynamical equation.
     *  \param environmentUpdater Object which is used to update time-dependent environment models to current time and state,
     *  must be consistent with member environment updaters of stateDerivativeModels entries.
     */
    DynamicsStateDerivativeModel(
            const std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels,
            const boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater ):
        environmentUpdater_( environmentUpdater )
    {
        std::vector< IntegratedStateType > stateTypeList;
        totalStateSize_ = 0;

        // Iterate over vector of state derivative models, check validity, set member variable map stateDerivativeModels_ and size indices.
        for( unsigned int i = 0; i < stateDerivativeModels.size( ); i++ )
        {
            if( i > 0 )
            {
                if( ( std::find( stateTypeList.begin( ), stateTypeList.end( ), stateDerivativeModels.at( i )->getIntegratedStateType( ) ) !=
                        stateTypeList.end( ) ) &&
                        ( stateDerivativeModels.at( i )->getIntegratedStateType( ) !=
                          stateDerivativeModels.at( i - 1 )->getIntegratedStateType( ) ) )
                {
                    std::cerr<<"Warning when making hybrid state derivative models, state type "<<
                               stateDerivativeModels.at( i )->getIntegratedStateType( )<<
                               " entries are non-contiguous"<<std::endl;
                }
            }

            // Check uniqueness of state derivative type calculator in list
            if( std::find( stateTypeList.begin( ), stateTypeList.end( ), stateDerivativeModels.at( i )->getIntegratedStateType( ) ) ==
                    stateTypeList.end( ) )
            {
                stateTypeSize_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] = 0;
                stateTypeStartIndex_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] = totalStateSize_;
            }
            stateTypeList.push_back( stateDerivativeModels.at( i )->getIntegratedStateType( ) );

            // Set state part sizes
            stateIndices_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ].push_back(
                        std::make_pair( totalStateSize_, stateDerivativeModels.at( i )->getStateSize( ) ) );
            totalStateSize_ += stateDerivativeModels.at( i )->getStateSize( );

            stateTypeSize_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ] += stateDerivativeModels.at( i )->getStateSize( );

            // Set current model in member map.
            stateDerivativeModels_[ stateDerivativeModels.at( i )->getIntegratedStateType( ) ].push_back( stateDerivativeModels.at( i ) );
        }
    }


    //! Function to calculate the system state derivative
    /*!
     *  Function to calculate the system state derivative, with settings as by last call to setPropagationSettings function. Dimensions of
     *  state must be consistent with these settings. Depending on the settings, this function may calculate the dynamical equations and/or variational
     *  equations for a subset of the dynamical equation types that are set in the stateDerivativeModels_ map.
     *  \param time Current time.
     *  \param state Current complete state.
     *  \return Calculated state derivative.
     */
    StateType computeStateDerivative( const TimeType time, const StateType& state )
    {
        // Initialize state derivative
        StateType stateDerivative = StateType::Zero( state.rows( ), state.cols( ) );

        // If dynamical equations are integrated, update the environment with the current state.
        if( evaluateDynamicsEquations_ )
        {
            std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > statesPerTypeInConventionalRepresentation =
                    getStatesPerTypeInConventionalRepresentation( state, time );
            environmentUpdater_->updateEnvironment( time, statesPerTypeInConventionalRepresentation, integratedStatesFromEnvironment_ );
        }
        else
        {
            environmentUpdater_->updateEnvironment( time, std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( ),
                                                    integratedStatesFromEnvironment_ );
            //std::cerr<<"Error, no evaluation of dynamical equations not yet implemented"<<std::endl;
        }



        std::pair< int, int > currentIndices;


        // If dynamical equations are integrated, evaluate dynamics state derivatives.
        if( evaluateDynamicsEquations_ )
        {
            // Iterate over all types of equations.
            for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( ); stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
                 stateDerivativeModelsIterator_++ )

            {
                for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
                {
                    // Update state derivative models
                    stateDerivativeModelsIterator_->second.at( i )->updateStateDerivativeModel( time );

                    // Evaluate and set current dynamical state derivative
                    currentIndices = stateIndices_.at( stateDerivativeModelsIterator_->first ).at( i );

                    stateDerivativeModelsIterator_->second.at( i )->calculateSystemStateDerivative(
                        time, state.block( currentIndices.first, startColumn_, currentIndices.second, 1 ) );

                    stateDerivative.block( currentIndices.first, startColumn_, currentIndices.second, 1 ) =
                            stateDerivativeModelsIterator_->second.at( i )->calculateSystemStateDerivative(
                                time, state.block( currentIndices.first, startColumn_, currentIndices.second, 1 ) );
                }
            }
        }

        return stateDerivative;
    }

    Eigen::MatrixXd computeStateDoubleDerivative(
            const double time, const Eigen::MatrixXd& state )
    {
        return computeStateDerivative( static_cast< TimeType >( time ), state.template cast< StateScalarType >( ) ).template cast< double >( );
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& outputState, const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > internalState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( outputState.rows( ), 1 );
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( ); stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            std::vector< std::pair< int, int > > currentStateIndices = stateIndices_.at( stateDerivativeModelsIterator_->first );
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                internalState.segment( currentStateIndices.at( i ).first, currentStateIndices.at( i ).second ) =
                        stateDerivativeModelsIterator_->second.at( i )->convertFromOutputSolution(
                            outputState.segment( currentStateIndices.at( i ).first, currentStateIndices.at( i ).second ), time );
            }
        }
        return internalState;
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > outputState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( internalSolution.rows( ), 1 );
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( ); stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            std::vector< std::pair< int, int > > currentStateIndices = stateIndices_.at( stateDerivativeModelsIterator_->first );
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                outputState.segment( currentStateIndices.at( i ).first, currentStateIndices.at( i ).second ) =
                        stateDerivativeModelsIterator_->second.at( i )->convertToOutputSolution(
                            internalSolution.segment( currentStateIndices.at( i ).first, currentStateIndices.at( i ).second ), time );
            }
        }
        return outputState;
    }


    void setPropagationSettings(
            const std::vector< IntegratedStateType >& stateTypesToNotIntegrate,
            const bool evaluateDynamicsEquations )
    {
        integratedStatesFromEnvironment_ = stateTypesToNotIntegrate;
        evaluateDynamicsEquations_ = evaluateDynamicsEquations;
        startColumn_ = 0;
    }

    std::map< IntegratedStateType, std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >
    getStateDerivativeModels( )
    {
        return stateDerivativeModels_;
    }

    std::map< IntegratedStateType, int > getStateTypeStartIndices( )
    {
        return stateTypeStartIndex_;
    }

    void updateStateDerivativeModelSettings(
            const simulation_setup::NamedBodyMap& bodyMap,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates,
            const int currentStateArcIndex )
    {
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( ); stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            switch( stateDerivativeModelsIterator_->first )
            {
            case transational_state:
                break;
            default:
                break;
            }
        }
    }


private:

    std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getStatesPerTypeInConventionalRepresentation(
            const StateType& state, const TimeType& time )
    {
        int startColumn = 0;

        std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitConventionalStates;

        std::pair< int, int > currentIndices;
        for( stateDerivativeModelsIterator_ = stateDerivativeModels_.begin( ); stateDerivativeModelsIterator_ != stateDerivativeModels_.end( );
             stateDerivativeModelsIterator_++ )
        {
            splitConventionalStates[ stateDerivativeModelsIterator_->first ] = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                        stateTypeSize_.at( stateDerivativeModelsIterator_->first ), 1 );

            int currentStateTypeSize = 0;
            for( unsigned int i = 0; i < stateDerivativeModelsIterator_->second.size( ); i++ )
            {
                currentIndices = stateIndices_.at( stateDerivativeModelsIterator_->first ).at( i );
                splitConventionalStates[ stateDerivativeModelsIterator_->first ].block( currentStateTypeSize, 0, currentIndices.second, 1 ) =
                        stateDerivativeModelsIterator_->second.at( i )->convertCurrentStateToGlobalRepresentation(
                            state.block( currentIndices.first, startColumn, currentIndices.second, 1 ), time );
                currentStateTypeSize += currentIndices.second;
            }
        }
        return splitConventionalStates;
    }

    boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater_;

    std::map< IntegratedStateType, std::vector< std::pair< int, int > > > stateIndices_;

    std::map< IntegratedStateType, int > stateTypeSize_;

    std::map< IntegratedStateType, int > stateTypeStartIndex_;

    std::map< IntegratedStateType, std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >
    stateDerivativeModels_;

    typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >::iterator
    stateDerivativeModelsIterator_;

    int totalStateSize_;

    std::vector< IntegratedStateType > integratedStatesFromEnvironment_;

    bool evaluateDynamicsEquations_;

    int startColumn_;


};

template< typename TimeType = double, typename StateScalarType = double,
          typename ConversionClassType = DynamicsStateDerivativeModel< TimeType, StateScalarType > >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > convertNumericalStateSolutionsToOutputSolutions(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& rawSolution,
        boost::shared_ptr< ConversionClassType > converterClass )
{
    // Initialize converted solution.
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > convertedSolution;

    // Iterate over all times.
    for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator stateIterator =
         rawSolution.begin( ); stateIterator != rawSolution.end( ); stateIterator++ )
    {
        // Convert solution at this time to output (typically ephemeris frame of given body) solution
        convertedSolution[ stateIterator->first ] = converterClass->convertToOutputSolution( stateIterator->second, stateIterator->first );
    }
    return convertedSolution;
}

} // namespace state_derivative_models
} // namespace tudat


#endif // HYBRIDSTATEDERIVATIVEMODEL_H
