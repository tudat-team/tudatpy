//
// Created by dominic on 22-12-22.
//

#ifndef TUDAT_PROPAGATIONRESULTS_H
#define TUDAT_PROPAGATIONRESULTS_H

#include <map>
#include <string>

#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

namespace tudat
{

    template<typename A, typename TimeType = double>
    std::map <TimeType, TimeType > getTimeStepHistory( const std::map< TimeType, A > stateHistory )
    {
        std::map <TimeType, TimeType > timeStepHistory;
        auto it_pre = stateHistory.begin( );
        auto it_post = stateHistory.begin( );
        it_post++;

        while( it_post != stateHistory.end() )
        {
            timeStepHistory[ it_pre->first ] = it_post->first - it_pre->first;
            it_pre++;
            it_post++;
        }
        return timeStepHistory;
    }

    namespace propagators
    {
        template<typename StateScalarType = double, typename TimeType = double>
        class SimulationResults
        {
        public:
            SimulationResults( ) 
            {}

            virtual ~SimulationResults( ) 
            {}

            virtual std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( ) = 0;

        };


        template<typename StateScalarType, typename TimeType>
        class SingleArcDynamicsSimulator;

        template <template<class, class> class SingleArcResults, class StateScalarType, class TimeType>
        class MultiArcSimulationResults;

        //! Object that holds the numerical results of the propagation of a single-arc propagation of dynamics
        /*
         *  Object that holds the numerical results of the propagation of a single-arc propagation of dynamics,
         *  this object contains the raw (as propagated) results, as well as the processed results (e.g. Cartesian
         *  elements for translational state), in addition to number of function evaluations, cumulative propagation
         *  time etc., if requested by user
         */
        template<typename StateScalarType = double, typename TimeType = double >
        class SingleArcSimulationResults : public SimulationResults<StateScalarType, TimeType>
        {
        public:


            static const bool is_variational = false;
            static const int number_of_columns = 1;

            SingleArcSimulationResults(const std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > > integratedStateAndBodyList,
                                       const std::shared_ptr <SingleArcPropagatorProcessingSettings> &outputSettings,
                                       const std::function< void ( std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&,
                                                                   const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& ) > rawSolutionConversionFunction,
                                       const std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > dependentVariableInterface,
                                       const bool sequentialPropagation = true ) :
                    SimulationResults<StateScalarType, TimeType>(),
                    processedStateIds_(getProcessedStateStrings( integratedStateAndBodyList ) ),
                    propagatedStateIds_( getPropagatedStateStrings( integratedStateAndBodyList ) ),
                    integratedStateAndBodyList_( integratedStateAndBodyList ),
                    outputSettings_(outputSettings),
                    dependentVariableInterface_( dependentVariableInterface ),
                    sequentialPropagation_( sequentialPropagation ),
                    rawSolutionConversionFunction_( rawSolutionConversionFunction ),
                    propagationIsPerformed_(false),
                    solutionIsCleared_( false ),
                    onlyProcessedSolutionSet_( false ),
                    propagationTerminationReason_(
                            std::make_shared<PropagationTerminationDetails>(propagation_never_run))
            {
            }

            //! Function that resets the state of this object, typically to signal that a new propagation is to be performed.
            void reset( ) 
            {
                clearSolutionMaps( );
                propagationIsPerformed_ = false;
                solutionIsCleared_ = false;
                onlyProcessedSolutionSet_ = false;
                propagationTerminationReason_ = std::make_shared<PropagationTerminationDetails>(propagation_never_run);
            }
            
            void manuallySetSecondaryData( const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > resultsToCopy )
            {
                dependentVariableHistory_ = resultsToCopy->getDependentVariableHistory( );
                cumulativeComputationTimeHistory_ =  resultsToCopy->getCumulativeComputationTimeHistory( );
                cumulativeNumberOfFunctionEvaluations_ =  resultsToCopy->getCumulativeNumberOfFunctionEvaluations( );
                propagationTerminationReason_ = resultsToCopy->getPropagationTerminationReason( );
                propagationIsPerformed_ = true;
            }

            //! Function that sets new numerical results of a propagation, after the propagation of the dynamics
            void reset(
                    const std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >>& equationsOfMotionNumericalSolutionRaw,
                    const std::map <TimeType, Eigen::VectorXd>& dependentVariableHistory,
                    const std::map<TimeType, double>& cumulativeComputationTimeHistory,
                    const std::map<TimeType, unsigned int>& cumulativeNumberOfFunctionEvaluations,
                    std::shared_ptr <PropagationTerminationDetails> propagationTerminationReason )
            {
                if ( sequentialPropagation_ || !isPropagationOngoing_ )
                {
                    reset( );
                    equationsOfMotionNumericalSolutionRaw_ = equationsOfMotionNumericalSolutionRaw;
                    dependentVariableHistory_ = dependentVariableHistory;
                    cumulativeComputationTimeHistory_ = cumulativeComputationTimeHistory;
                    cumulativeNumberOfFunctionEvaluations_ = cumulativeNumberOfFunctionEvaluations;

                    if ( !sequentialPropagation_ )
                    {
                        isPropagationOngoing_ = true;
                    }
                }
                else if ( !sequentialPropagation_ && isPropagationOngoing_ )
                {
                    equationsOfMotionNumericalSolutionRaw_.insert( equationsOfMotionNumericalSolutionRaw.begin( ), equationsOfMotionNumericalSolutionRaw.end( ) );
                    dependentVariableHistory_.insert( dependentVariableHistory.begin( ), dependentVariableHistory.end( ) );
                    cumulativeComputationTimeHistory_.insert ( cumulativeComputationTimeHistory.begin( ), cumulativeComputationTimeHistory.end( ) );
                    cumulativeNumberOfFunctionEvaluations_.insert( cumulativeNumberOfFunctionEvaluations.begin( ), cumulativeNumberOfFunctionEvaluations.end( ) );
                    isPropagationOngoing_ = false;
                }
                rawSolutionConversionFunction_( equationsOfMotionNumericalSolution_, equationsOfMotionNumericalSolutionRaw_ );
                propagationTerminationReason_ = propagationTerminationReason;

            }

            //! Function to clear all maps with numerical results, but *not* signal that a new propagation will start,
            //! this is typically done to save memory usage (and is called using the clearNumericalSolution setting
            //! of the PropagatorProcessingSettings
            void clearSolutionMaps( )
            {
                equationsOfMotionNumericalSolution_.clear();
                equationsOfMotionNumericalSolutionRaw_.clear();
                dependentVariableHistory_.clear();
                cumulativeComputationTimeHistory_.clear();
                cumulativeNumberOfFunctionEvaluations_.clear();
                solutionIsCleared_ = true;
            }

            //! Get initial and final propagation time from raw results
            std::pair< TimeType, TimeType > getArcInitialAndFinalTime( )
            {
                if( equationsOfMotionNumericalSolutionRaw_.size( ) == 0 )
                {
                    throw std::runtime_error( "Error when getting single-arc dynamics initial and final times; no results set" );
                }
                return std::make_pair( equationsOfMotionNumericalSolutionRaw_.begin( )->first, equationsOfMotionNumericalSolutionRaw_.rbegin( )->first );
            }

            //! Function to signal that propagation is finished, and add number of function evaluations
            void finalizePropagation( const std::map<TimeType, unsigned int> cumulativeNumberOfFunctionEvaluations )
            {
                cumulativeNumberOfFunctionEvaluations_ = cumulativeNumberOfFunctionEvaluations;
                propagationIsPerformed_ = true;
            }

            //! Manually set processed numerical solution, to be used when this object *is not* used in the
            //! propagation loop, but *is* used to store the numerical results
            void setEquationsOfMotionNumericalSolution(
                    const std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>> & equationsOfMotionNumericalSolution )
            {
                onlyProcessedSolutionSet_ = true;
                equationsOfMotionNumericalSolution_ = equationsOfMotionNumericalSolution;
            }

            //! Function to check if output map that is requested is available
            void checkAvailabilityOfSolution( const std::string& dataToRetrieve, const bool checkEomOnly = true )
            {
                if( !propagationIsPerformed_ )
                {
                    throw std::runtime_error( "Error when retrieving " + dataToRetrieve + ", propagation is not yet performed." );
                }
                else if( solutionIsCleared_ )
                {
                    throw std::runtime_error( "Error when retrieving " + dataToRetrieve + ", propagation has been performed, but results have been cleared." );
                }
                else if( checkEomOnly && onlyProcessedSolutionSet_ )
                {
                    throw std::runtime_error( "Error when retrieving " + dataToRetrieve + ", propagation has been performed using other object; current object only holds equations of motion solution." );
                }
            }

            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>> &
            getEquationsOfMotionNumericalSolution( )
            {
                if( !onlyProcessedSolutionSet_ )
                {
                    checkAvailabilityOfSolution( "equations of motion numerical solution", false );
                }
                return equationsOfMotionNumericalSolution_;
            }

            std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolutionDouble( )
            {
                return utilities::staticCastMapKeys< double, TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( equationsOfMotionNumericalSolution_ );
            }

            std::pair< std::vector< double >, std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionDoubleSplit( )
            {
                std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > stateMap = getEquationsOfMotionNumericalSolutionDouble( );
                return std::make_pair( utilities::createVectorFromMapKeys( stateMap ), utilities::createVectorFromMapValues( stateMap ) );
            }

            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>> &
            getEquationsOfMotionNumericalSolutionRaw( )
            {
                checkAvailabilityOfSolution( "equations of motion unprocessed numerical solution" );
                return equationsOfMotionNumericalSolutionRaw_;
            }

            std::map <TimeType, Eigen::VectorXd> &getDependentVariableHistory( )
            {
                checkAvailabilityOfSolution( "dependent variable history", false );
                return dependentVariableHistory_;
            }

            std::map<TimeType, double> &getCumulativeComputationTimeHistory( )
            {
                checkAvailabilityOfSolution( "cumulative computation time history", false );
                return cumulativeComputationTimeHistory_;
            }

            double getTotalComputationRuntime( )
            {
                checkAvailabilityOfSolution( "cumulative computation time history", false );
                return std::max( cumulativeComputationTimeHistory_.begin( )->second,
                                 cumulativeComputationTimeHistory_.rbegin( )->second );
            }

            std::map<TimeType, unsigned int> &getCumulativeNumberOfFunctionEvaluations( )
            {
                checkAvailabilityOfSolution( "cumulative number of function evaluations", false );
                return cumulativeNumberOfFunctionEvaluations_;
            }

            double getTotalNumberOfFunctionEvaluations( )
            {
                checkAvailabilityOfSolution( "cumulative number of function evaluations", false );
                return std::max( cumulativeNumberOfFunctionEvaluations_.begin( )->second,
                                 cumulativeNumberOfFunctionEvaluations_.rbegin( )->second );
            }

            std::shared_ptr <PropagationTerminationDetails> getPropagationTerminationReason( )
            {
                return propagationTerminationReason_;
            }

            bool integrationCompletedSuccessfully( ) const
            {
                return (propagationTerminationReason_->getPropagationTerminationReason() ==
                        termination_condition_reached);
            }


            std::map <std::pair<int, int>, std::string> getDependentVariableId( )  const
            {
                if( dependentVariableInterface_ != nullptr )
                {
                    return dependentVariableInterface_->getDependentVariableIds( );
                }
                else
                {
                    return std::map <std::pair<int, int>, std::string>( );
                }
            }

            std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > getOrderedDependentVariableSettings( ) const
            {
                if( dependentVariableInterface_ != nullptr )
                {
                    return dependentVariableInterface_->getOrderedDependentVariableSettings( );
                }
                else
                {
                    return std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > >( );
                }
            }

            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > getOriginalDependentVariableSettings( ) const
            {
                if( dependentVariableInterface_ != nullptr )
                {
                    return dependentVariableInterface_->getDependentVariablesSettings( );
                }
                else
                {
                    return std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( );
                }
            }

            std::map <std::pair<int, int>, std::string> getProcessedStateIds( )
            {
                return processedStateIds_;
            }

            std::map <std::pair<int, int>, std::string> getPropagatedStateIds( )
            {
                return propagatedStateIds_;
            }

            int getPropagatedStateSize( )
            {
                return propagatedStateIds_.rbegin( )->first.first + propagatedStateIds_.rbegin( )->first.second;
            }

            std::shared_ptr <SingleArcPropagatorProcessingSettings> getOutputSettings( )
            {
                return outputSettings_;
            }

            bool getPropagationIsPerformed( )
            {
                return propagationIsPerformed_;
            }

            bool getSolutionIsCleared( )
            {
                return solutionIsCleared_;
            }

            bool isPropagatedAndProcessedStateEqual( )
            {
                bool areEqual = true;
                for( auto it : integratedStateAndBodyList_ )
                {
                    for( unsigned int i = 0; i < it.second.size( ); i++ )
                    {
                        if( it.first == translational_state && std::get< 2 >( it.second.at( i ) ).translationalPropagatorType_ != cowell )
                        {
                            areEqual = false;
                        }
                        else if( it.first == rotational_state && std::get< 2 >( it.second.at( i ) ).rotationalPropagatorType_ != quaternions )
                        {
                            areEqual = false;
                        }
                    }
                }
                return areEqual;
            }

            void updateDependentVariableInterface( )
            {
                if( dependentVariableHistory_.size( ) > 0 && dependentVariableInterface_ != nullptr )
                {
                    std::shared_ptr< interpolators::LagrangeInterpolator< TimeType, Eigen::VectorXd > > dependentVariablesInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< TimeType, Eigen::VectorXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::VectorXd, TimeType >( dependentVariableHistory_ ),
                                    utilities::createVectorFromMapValues< Eigen::VectorXd, TimeType >( dependentVariableHistory_ ), 8 );
                    dependentVariableInterface_->updateDependentVariablesInterpolator( dependentVariablesInterpolator );
                }
            }

            std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > getSingleArcDependentVariablesInterface( )
            {
                return dependentVariableInterface_;
            }

            std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( )
            {
                return getSingleArcDependentVariablesInterface( );
            }


        private:

            //! Map of state history of numerically integrated bodies.
            /*!
             *  Map of state history of numerically integrated bodies, i.e. the result of the numerical integration, transformed
             *  into the 'conventional form' (\sa SingleStateTypeDerivative::convertToOutputSolution). Key of map denotes time,
             *  values are concatenated vectors of integrated body states (order defined by propagatorSettings_).
             *  NOTE: this map is empty if clearNumericalSolutions_ is set to true.
             */
            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution_;

            //! Map of state history of numerically integrated bodies.
            /*!
            *  Map of state history of numerically integrated bodies, i.e. the result of the numerical integration, in the
            *  original propagation coordinates. Key of map denotes time, values are concatenated vectors of integrated body
            * states (order defined by propagatorSettings_).
            *  NOTE: this map is empty if clearNumericalSolutions_ is set to true.
            */
            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>> equationsOfMotionNumericalSolutionRaw_;

            //! Map of dependent variable history that was saved during numerical propagation.
            std::map <TimeType, Eigen::VectorXd> dependentVariableHistory_;

            //! Map of cumulative computation time history that was saved during numerical propagation.
            std::map<TimeType, double> cumulativeComputationTimeHistory_;

            //! Map of cumulative number of function evaluations that was saved during numerical propagation.
            std::map<TimeType, unsigned int> cumulativeNumberOfFunctionEvaluations_;

            std::map <std::pair<int, int>, std::string> processedStateIds_;

            std::map <std::pair<int, int>, std::string> propagatedStateIds_;

            std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > > integratedStateAndBodyList_;

            std::shared_ptr <SingleArcPropagatorProcessingSettings> outputSettings_;

            std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > dependentVariableInterface_;



            //! Bool denoting whether the propagation is sequential or bidirectional (default is true).
            bool sequentialPropagation_;

            //! Function to convert the propagated solution to conventional solution (see DynamicsStateDerivativeModel::convertToOutputSolution)
            const std::function< void ( std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&,
                                        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& ) > rawSolutionConversionFunction_;

            bool propagationIsPerformed_;

            bool solutionIsCleared_;

            bool onlyProcessedSolutionSet_;

            //! Event that triggered the termination of the propagation
            std::shared_ptr <PropagationTerminationDetails> propagationTerminationReason_;

            friend class SingleArcDynamicsSimulator<StateScalarType, TimeType>;

//            friend class MultiArcSimulationResults<StateScalarType, TimeType, NumberOfStateColumns >;

            //! Boolean denoting whether the full propagation has been fully completed or is ongoing (for non sequential propagations only)
            bool isPropagationOngoing_ = false;

        };


        template<typename StateScalarType = double, typename TimeType = double >
        class SingleArcVariationalSimulationResults: public SimulationResults< StateScalarType, TimeType >
        {
        public:

            static const bool is_variational = true;
            static const int number_of_columns = Eigen::Dynamic;

            SingleArcVariationalSimulationResults( const std::shared_ptr <SingleArcSimulationResults< StateScalarType, TimeType > >  singleArcDynamicsResults,
                                                   const int stateTransitionMatrixSize,
                                                   const int sensitivityMatrixSize ):
            singleArcDynamicsResults_( singleArcDynamicsResults ),
            stateTransitionMatrixSize_( stateTransitionMatrixSize ),
            sensitivityMatrixSize_( sensitivityMatrixSize )
            {

            }

            void reset( )
            {
                clearSolutionMaps( );
                singleArcDynamicsResults_->reset( );
            }

            void reset(
                    std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >>& fullSolution,
                    const std::map <TimeType, Eigen::VectorXd>& dependentVariableHistory,
                    const std::map<TimeType, double>& cumulativeComputationTimeHistory,
                    const std::map<TimeType, unsigned int>& cumulativeNumberOfFunctionEvaluations,
                    std::shared_ptr <PropagationTerminationDetails> propagationTerminationReason )
            {
                std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >> equationsOfMotionNumericalSolutionRaw;
                splitSolution( fullSolution, equationsOfMotionNumericalSolutionRaw );
                singleArcDynamicsResults_->reset(
                        equationsOfMotionNumericalSolutionRaw,
                        dependentVariableHistory,
                        cumulativeComputationTimeHistory,
                        cumulativeNumberOfFunctionEvaluations,
                        propagationTerminationReason );
            }

            void manuallySetSecondaryData( const std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > resultsToCopy )
            {
                singleArcDynamicsResults_->manuallySetSecondaryData( resultsToCopy->getDynamicsResults( ) );
            }

            //! Function to split the full numerical solution into the solution for state transition matrix, sensitivity matrix, and unprocessed dynamics solution
            void splitSolution(
                    const std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >>& fullSolution,
                    std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >>& equationsOfMotionNumericalSolutionRaw )
            {
                for( auto it : fullSolution )
                {
                    stateTransitionSolution_[ static_cast< double >( it.first ) ] = it.second.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).template cast< double >( );
                    sensitivitySolution_[ static_cast< double >( it.first ) ] = it.second.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ).template cast< double >( );
                    equationsOfMotionNumericalSolutionRaw[ static_cast< double >( it.first ) ] = it.second.block( 0, stateTransitionMatrixSize_ + sensitivityMatrixSize_, stateTransitionMatrixSize_, 1 );
                }
            }

            void clearSolutionMaps( )
            {
                singleArcDynamicsResults_->clearSolutionMaps( );
                stateTransitionSolution_.clear( );
                sensitivitySolution_.clear( );
            }

            void finalizePropagation( const std::map<TimeType, unsigned int> cumulativeNumberOfFunctionEvaluations )
            {
                singleArcDynamicsResults_->finalizePropagation( cumulativeNumberOfFunctionEvaluations );
            }

            std::map < double, Eigen::MatrixXd >& getStateTransitionSolution( )
            {
                return stateTransitionSolution_;
            }

            std::map < double, Eigen::MatrixXd >& getSensitivitySolution( )
            {
                return sensitivitySolution_;
            }

            const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > getDynamicsResults( )
            {
                return singleArcDynamicsResults_;
            }

            std::pair< TimeType, TimeType > getArcInitialAndFinalTime( )
            {
                if( stateTransitionSolution_.size( ) == 0 )
                {
                    throw std::runtime_error( "Error when getting single-arc variational initial and final times; no results set" );
                }
                return std::make_pair( stateTransitionSolution_.begin( )->first, stateTransitionSolution_.rbegin( )->first );
            }

            int getStateTransitionMatrixSize( )
            {
                return stateTransitionMatrixSize_;
            }

            int getSensitivityMatrixSize( )
            {
                return sensitivityMatrixSize_;
            }

            std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > getSingleArcDependentVariablesInterface( )
            {
                return singleArcDynamicsResults_->getSingleArcDependentVariablesInterface( );
            }

            std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( )
            {
                return getSingleArcDependentVariablesInterface( );
            }


        protected:
            const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > singleArcDynamicsResults_;

            const int stateTransitionMatrixSize_;

            const int sensitivityMatrixSize_;

            std::map < double, Eigen::MatrixXd > stateTransitionSolution_;

            std::map < double, Eigen::MatrixXd > sensitivitySolution_;
        };

        template< typename SimulationResults, typename StateScalarType = double, typename TimeType = double >
        class SingleArcResultsRetriever
        {
        public:

            static std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > getSingleArcSimulationResults(
                    const std::shared_ptr< SimulationResults > simulationResults );

        };

        template< typename StateScalarType, typename TimeType >
        class SingleArcResultsRetriever< SingleArcSimulationResults< StateScalarType, TimeType >, StateScalarType, TimeType >
        {
        public:

            static std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > getSingleArcSimulationResults(
                    const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > simulationResults )
            {
                return simulationResults;
            }

        };

        template< typename StateScalarType, typename TimeType >
        class SingleArcResultsRetriever< SingleArcVariationalSimulationResults< StateScalarType, TimeType >, StateScalarType, TimeType >
        {
        public:

            static std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > getSingleArcSimulationResults(
                    const std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > simulationResults )
            {
                return simulationResults->getDynamicsResults( );
            }

        };
        //! Class that holds numerical results for multi-arc simulations. This class may
        //! hold results for dynamics-only or variational+dynamics results. For the former,
        //! the SingleArcResults template argument is SingleArcSimulationResults, for the latter it is
        //! SingleArcVariationalSimulationResults
        template <template<class, class> class SingleArcResults, class StateScalarType, class TimeType>
        class MultiArcSimulationResults : public SimulationResults<StateScalarType, TimeType> {
        public:

            using single_arc_type = SingleArcResults< StateScalarType, TimeType >;

            MultiArcSimulationResults(
                    const std::vector< std::shared_ptr< SingleArcResults< StateScalarType, TimeType > > > singleArcResults,
                    const std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > dependentVariableInterface ):
                    singleArcResults_( singleArcResults ), propagationIsPerformed_( false ), solutionIsCleared_( false ),
                    dependentVariableInterface_( dependentVariableInterface )
                    {
                        if( dependentVariableInterface_ == nullptr )
                        {
                            throw std::runtime_error( "Error when creating MultiArcSimulationResults, dependentVariableInterface_ is NULL " );
                        }
                    }

            ~MultiArcSimulationResults( )
            {}

            bool getPropagationIsPerformed( )
            {
                return propagationIsPerformed_;
            }

            void restartPropagation( )
            {
                propagationIsPerformed_ = false;
                solutionIsCleared_ = false;
                arcStartTimes_.clear( );
                arcEndTimes_.clear( );

                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    singleArcResults_.at( i )->reset( );
                }
            }

            void setPropagationIsPerformed( )
            {
                propagationIsPerformed_ = true;
                arcStartTimes_.clear( );
                arcEndTimes_.clear( );
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    auto initialAndFinalTime = singleArcResults_.at( i )->getArcInitialAndFinalTime( );
                    arcStartTimes_.push_back( initialAndFinalTime.first );
                    arcEndTimes_.push_back( initialAndFinalTime.second );
                }
            }

//            void setPropagationIsPerformed( const std::vector< double >& arcStartTimes )
//            {
//                propagationIsPerformed_ = true;
//                arcStartTimes_ = arcStartTimes;
//            }

            void manuallySetPropagationResults(
                    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > numericalMultiArcSolution )
            {
                for( unsigned int i = 0; i < numericalMultiArcSolution.size( ); i++ )
                {
                    singleArcResults_.at( i )->setEquationsOfMotionNumericalSolution( numericalMultiArcSolution.at( i ) );
                }
                propagationIsPerformed_ = true;

                arcStartTimes_.clear( );
                arcEndTimes_.clear( );
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    arcStartTimes_.push_back( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolution( ).begin( )->first );
                    arcEndTimes_.push_back( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolution( ).rbegin( )->first );
                }
            }

            void manuallySetSecondaryData( const std::shared_ptr< MultiArcSimulationResults< SingleArcResults, StateScalarType, TimeType > > resultsToCopy )
            {
                if( resultsToCopy->getSingleArcResults( ).size( ) != singleArcResults_.size( ) )
                {
                    throw std::runtime_error( "Error when manually resetting multi-arc secondary data; arc sizes are incompatible" );
                }

                for( unsigned int i = 0; i < resultsToCopy->getSingleArcResults( ).size( ); i++ )
                {
                    singleArcResults_.at( i )->manuallySetSecondaryData( resultsToCopy->getSingleArcResults( ).at( i ) );
                }
            }

            std::vector< std::shared_ptr< SingleArcResults< StateScalarType, TimeType > > > getSingleArcResults( )
            {
                return singleArcResults_;
            }

            std::vector< double > getArcStartTimes( )
            {
                if( !propagationIsPerformed_ )
                {
                    throw std::runtime_error( "Error when getting multi-arc initial times; propagation not yet performed" );
                }
                return arcStartTimes_;
            }

            std::vector< double > getArcEndTimes( )
            {
                if( !propagationIsPerformed_ )
                {
                    throw std::runtime_error( "Error when getting multi-arc final times; propagation not yet performed" );
                }
                return arcEndTimes_;
            }

            bool getSolutionIsCleared( )
            {
                return solutionIsCleared_;
            }

            bool integrationCompletedSuccessfully() const
            {
                bool isSuccesful = true;
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    if( !singleArcResults_.at( i )->integrationCompletedSuccessfully( ) )
                    {
                        isSuccesful = false;
                        break;
                    }
                }
                return isSuccesful;
            }

            void clearSolutionMaps( )
            {
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    singleArcResults_.at( i )->clearSolutionMaps( );
                }
                solutionIsCleared_ = true;
            }

            std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getConcatenatedEquationsOfMotionResults(
                    const bool clearResults = false )
            {
                std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > concatenatedResults;
                if( !solutionIsCleared_ )
                {
                    for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                    {
                        if( clearResults )
                        {
                            concatenatedResults.push_back( std::move( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolution( ) ) );
                            singleArcResults_.at( i )->clearSolutionMaps( );
                        }
                        else
                        {
                            concatenatedResults.push_back( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolution( ) );
                        }
                    }
                    if( clearResults )
                    {
                        solutionIsCleared_ = true;
                    }
                }
                return concatenatedResults;
            }

            std::vector< std::map< TimeType, Eigen::VectorXd > > getConcatenatedDependentVariableResults( )
            {
                std::vector< std::map< TimeType, Eigen::VectorXd > > concatenatedResults;
                if( !solutionIsCleared_ )
                {
                    for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                    {
                        concatenatedResults.push_back( singleArcResults_.at( i )->getDependentVariableHistory( ) );
                    }
                }
                return concatenatedResults;
            }

            std::vector< std::map< TimeType, double > > getConcatenatedCumulativeComputationTimeHistory( )
            {
                std::vector< std::map< TimeType, double > > concatenatedResults;
                if( !solutionIsCleared_ )
                {
                    for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                    {
                        concatenatedResults.push_back( singleArcResults_.at( i )->getCumulativeComputationTimeHistory( ) );
                    }
                }
                return concatenatedResults;
            }

            std::vector< std::shared_ptr< PropagationTerminationDetails > > getConcatenatedTerminationReasons( )
            {
                std::vector< std::shared_ptr< PropagationTerminationDetails > > concatenatedResults;
                if( !solutionIsCleared_ )
                {
                    for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                    {
                        concatenatedResults.push_back( singleArcResults_.at( i )->getPropagationTerminationReason( ) );
                    }
                }
                return concatenatedResults;
            }

            std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > getMultiArcDependentVariablesInterface( )
            {
                return dependentVariableInterface_;
            }

            std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( )
            {
                return getMultiArcDependentVariablesInterface( );
            }

            void updateDependentVariableInterface( )
            {
                std::vector<std::shared_ptr<interpolators::OneDimensionalInterpolator<TimeType, Eigen::VectorXd> > > dependentVariablesInterpolators;
                for ( unsigned int i = 0; i < arcStartTimes_.size( ); i++ )
                {
                    if( singleArcResults_.at( i )->getDependentVariableHistory( ).size( ) > 0 )
                    {
                        std::shared_ptr<interpolators::LagrangeInterpolator<TimeType, Eigen::VectorXd> >
                            dependentVariablesInterpolator =
                            std::make_shared<interpolators::LagrangeInterpolator<TimeType, Eigen::VectorXd> >(
                                utilities::createVectorFromMapKeys<Eigen::VectorXd, TimeType>(
                                    singleArcResults_.at( i )->getDependentVariableHistory( )),
                                utilities::createVectorFromMapValues<Eigen::VectorXd, TimeType>(
                                    singleArcResults_.at( i )->getDependentVariableHistory( )), 8 );
                        dependentVariablesInterpolators.push_back( dependentVariablesInterpolator );
                    }
                    else
                    {
                        dependentVariablesInterpolators.push_back( nullptr );
                    }

                }

                // Update arc end times
                dependentVariableInterface_->updateDependentVariablesInterpolators(
                        dependentVariablesInterpolators, arcStartTimes_, arcEndTimes_ );
            }

        private:
            const std::vector< std::shared_ptr< SingleArcResults< StateScalarType, TimeType > > > singleArcResults_;

            bool propagationIsPerformed_;

            bool solutionIsCleared_;

            //! List of start times of each arc. NOTE:   This list is updated after every propagation.
            std::vector< double > arcStartTimes_;

            std::vector< double > arcEndTimes_;

            std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > dependentVariableInterface_;


        };


        template <template<class, class> class SingleArcResults, class StateScalarType, class TimeType>
        class HybridArcSimulationResults : public SimulationResults<StateScalarType, TimeType>
        {
        public:
            HybridArcSimulationResults(
                    const std::shared_ptr< SingleArcResults< StateScalarType, TimeType > > singleArcResults,
                    const std::shared_ptr< MultiArcSimulationResults< SingleArcResults, StateScalarType, TimeType > > multiArcResults ):
                    singleArcResults_( singleArcResults ), multiArcResults_( multiArcResults )
            {
                dependentVariableInterface_ = std::make_shared<HybridArcDependentVariablesInterface<TimeType> >(
                        singleArcResults_->getSingleArcDependentVariablesInterface( ),
                        multiArcResults_->getMultiArcDependentVariablesInterface( ) );
            }

            ~HybridArcSimulationResults( ) 
            {}

            std::shared_ptr< SingleArcResults< StateScalarType, TimeType > > getSingleArcResults( )
            {
                return singleArcResults_;
            }

            std::shared_ptr< MultiArcSimulationResults< SingleArcResults, StateScalarType, TimeType > > getMultiArcResults( )
            {
                return multiArcResults_;
            }

            bool integrationCompletedSuccessfully() const
            {
                return singleArcResults_->integrationCompletedSuccessfully( ) && multiArcResults_->integrationCompletedSuccessfully( );
            }

            std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getConcatenatedEquationsOfMotionResults( )
            {
                std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > concatenatedResults;
                if( singleArcResults_->getSolutionIsCleared( ) != multiArcResults_->getSolutionIsCleared( ) )
                {
                    throw std::runtime_error( "Error when getting concatenated hybrid-arc results, constituent results have inconsistent cleared state." );
                }

                if( singleArcResults_->getPropagationIsPerformed() != multiArcResults_->getPropagationIsPerformed( ) )
                {
                    throw std::runtime_error( "Error when getting concatenated hybrid-arc results, constituent results are inconsistently propagated." );
                }
                if( !singleArcResults_->getSolutionIsCleared( ) )
                {
                    concatenatedResults.push_back( singleArcResults_->getEquationsOfMotionNumericalSolution( ) );
                    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
                            multiArcResults = multiArcResults_->getConcatenatedEquationsOfMotionResults( );
                    concatenatedResults.insert( concatenatedResults.end( ), multiArcResults.begin( ), multiArcResults.end( ) );
                }
                return concatenatedResults;
            }
            
            std::vector< std::map< TimeType, Eigen::VectorXd > > getConcatenatedDependentVariableResults( )
            {
                std::vector< std::map< TimeType, Eigen::VectorXd  > > concatenatedResults;
                if( singleArcResults_->getSolutionIsCleared( ) != multiArcResults_->getSolutionIsCleared( ) )
                {
                    throw std::runtime_error( "Error when getting concatenated hybrid-arc results, constituent results have inconsistent cleared state." );
                }

                if( singleArcResults_->getPropagationIsPerformed() != multiArcResults_->getPropagationIsPerformed( ) )
                {
                    throw std::runtime_error( "Error when getting concatenated hybrid-arc results, constituent results are inconsistently propagated." );
                }
                if( !singleArcResults_->getSolutionIsCleared( ) )
                {
                    concatenatedResults.push_back( singleArcResults_->getDependentVariableHistory( ) );
                    std::vector< std::map< TimeType, Eigen::VectorXd > >
                            multiArcResults = multiArcResults_->getConcatenatedDependentVariableResults( );
                    concatenatedResults.insert( concatenatedResults.end( ), multiArcResults.begin( ), multiArcResults.end( ) );
                }
                return concatenatedResults;
            }

            std::vector< std::map< TimeType, double > > getConcatenatedCumulativeComputationTimeHistory( )
            {
                std::vector< std::map< TimeType, double > > concatenatedResults;
                if( singleArcResults_->getSolutionIsCleared( ) != multiArcResults_->getSolutionIsCleared( ) )
                {
                    throw std::runtime_error( "Error when getting concatenated hybrid-arc results, constituent results have inconsistent cleared state." );
                }

                if( singleArcResults_->getPropagationIsPerformed() != multiArcResults_->getPropagationIsPerformed( ) )
                {
                    throw std::runtime_error( "Error when getting concatenated hybrid-arc results, constituent results are inconsistently propagated." );
                }
                if( !singleArcResults_->getSolutionIsCleared( ) )
                {
                    concatenatedResults.push_back( singleArcResults_->getCumulativeComputationTimeHistory( ) );
                    std::vector< std::map< TimeType, double > >
                            multiArcResults = multiArcResults_->getConcatenatedCumulativeComputationTimeHistory( );
                    concatenatedResults.insert( concatenatedResults.end( ), multiArcResults.begin( ), multiArcResults.end( ) );
                }
                return concatenatedResults;
            }

            std::shared_ptr< HybridArcDependentVariablesInterface< TimeType > > getHybridArcDependentVariablesInterface( )
            {
                return dependentVariableInterface_;
            }

            std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( )
            {
                return getHybridArcDependentVariablesInterface( );
            }

            void updateDependentVariableInterface( )
            {
                singleArcResults_->updateDependentVariableInterface( );
                multiArcResults_->updateDependentVariableInterface( );
            }

        protected:
            std::shared_ptr< SingleArcResults< StateScalarType, TimeType > > singleArcResults_;

            std::shared_ptr< MultiArcSimulationResults< SingleArcResults, StateScalarType, TimeType > > multiArcResults_;

            std::shared_ptr< HybridArcDependentVariablesInterface< TimeType > > dependentVariableInterface_;
        };

    }

}
#endif //TUDAT_PROPAGATIONRESULTS_H
