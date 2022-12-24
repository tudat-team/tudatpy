//
// Created by dominic on 22-12-22.
//

#ifndef TUDAT_PROPAGATIONRESULTS_H
#define TUDAT_PROPAGATIONRESULTS_H

#include <map>
#include <string>

#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"

namespace tudat
{

    namespace propagators
    {
        template<typename StateScalarType = double, typename TimeType = double>
        class SimulationResults
        {
        public:
            SimulationResults() {}

            virtual ~SimulationResults() {}

        };


        template<typename StateScalarType, typename TimeType>
        class SingleArcDynamicsSimulator;

        template<typename StateScalarType, typename TimeType, int NumberOfStateColumns >
        class MultiArcSimulationResults;

        template<typename StateScalarType = double, typename TimeType = double, int NumberOfStateColumns = 1 >
        class SingleArcSimulationResults : public SimulationResults<StateScalarType, TimeType>
        {
        public:

            SingleArcSimulationResults(const std::map <std::pair<int, int>, std::string> &dependentVariableIds,
                                       const std::map <std::pair<int, int>, std::string> &stateIds,
                                       const std::shared_ptr <SingleArcPropagatorProcessingSettings> &outputSettings,
                                       const std::function< void ( std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&,
                                                                   const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& ) > rawSolutionConversionFunction ) :
                    SimulationResults<StateScalarType, TimeType>(),
                    dependentVariableIds_(dependentVariableIds),
                    stateIds_(stateIds),
                    outputSettings_(outputSettings),
                    rawSolutionConversionFunction_( rawSolutionConversionFunction ),
                    propagationIsPerformed_(false),
                    solutionIsCleared_( false ),
                    propagationTerminationReason_(
                            std::make_shared<PropagationTerminationDetails>(propagation_never_run)) {
            }

            void reset() {
                equationsOfMotionNumericalSolution_.clear();
                equationsOfMotionNumericalSolutionRaw_.clear();
                dependentVariableHistory_.clear();
                cumulativeComputationTimeHistory_.clear();
                cumulativeNumberOfFunctionEvaluations_.clear();
                propagationIsPerformed_ = false;
                solutionIsCleared_ = false;
                propagationTerminationReason_ = std::make_shared<PropagationTerminationDetails>(propagation_never_run);
            }

            void reset(
                    const std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, NumberOfStateColumns >>& equationsOfMotionNumericalSolutionRaw,
                    const std::map <TimeType, Eigen::VectorXd>& dependentVariableHistory,
                    const std::map<TimeType, double>& cumulativeComputationTimeHistory,
                    const std::map<TimeType, unsigned int>& cumulativeNumberOfFunctionEvaluations,
                    std::shared_ptr <PropagationTerminationDetails> propagationTerminationReason )
            {
                reset( );
                equationsOfMotionNumericalSolutionRaw_ = equationsOfMotionNumericalSolutionRaw;
                rawSolutionConversionFunction_( equationsOfMotionNumericalSolution_, equationsOfMotionNumericalSolutionRaw_ );
                dependentVariableHistory_ = dependentVariableHistory;
                cumulativeComputationTimeHistory_ = cumulativeComputationTimeHistory;
                cumulativeNumberOfFunctionEvaluations_ = cumulativeNumberOfFunctionEvaluations;
                propagationTerminationReason_ = propagationTerminationReason;
            }


            void clearSolutionMaps( )
            {
                equationsOfMotionNumericalSolution_.clear();
                equationsOfMotionNumericalSolutionRaw_.clear();
                dependentVariableHistory_.clear();
                cumulativeComputationTimeHistory_.clear();
                cumulativeNumberOfFunctionEvaluations_.clear();
                solutionIsCleared_ = true;
            }

            void finalizePropagation( const std::map<TimeType, unsigned int> cumulativeNumberOfFunctionEvaluations )
            {
                cumulativeNumberOfFunctionEvaluations_ = cumulativeNumberOfFunctionEvaluations;
                propagationIsPerformed_ = true;
            }

            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1>> &
            getEquationsOfMotionNumericalSolution() {
                return equationsOfMotionNumericalSolution_;
            }

            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, NumberOfStateColumns>> &
            getEquationsOfMotionNumericalSolutionRaw() {
                return equationsOfMotionNumericalSolutionRaw_;
            }

            std::map <TimeType, Eigen::VectorXd> &getDependentVariableHistory() {
                return dependentVariableHistory_;
            }

            std::map<TimeType, double> &getCumulativeComputationTimeHistory() {
                return cumulativeComputationTimeHistory_;
            }

            std::map<TimeType, unsigned int> &getCumulativeNumberOfFunctionEvaluations() {
                return cumulativeNumberOfFunctionEvaluations_;
            }

            std::shared_ptr <PropagationTerminationDetails> getPropagationTerminationReason() {
                return propagationTerminationReason_;
            }

            bool integrationCompletedSuccessfully() const {
                return (propagationTerminationReason_->getPropagationTerminationReason() ==
                        termination_condition_reached);
            }


            std::map <std::pair<int, int>, std::string> getDependentVariableId() {
                return dependentVariableIds_;
            }

            std::map <std::pair<int, int>, std::string> getStateIds() {
                return stateIds_;
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


        private:

            //! Map of state history of numerically integrated bodies.
            /*!
             *  Map of state history of numerically integrated bodies, i.e. the result of the numerical integration, transformed
             *  into the 'conventional form' (\sa SingleStateTypeDerivative::convertToOutputSolution). Key of map denotes time,
             *  values are concatenated vectors of integrated body states (order defined by propagatorSettings_).
             *  NOTE: this map is empty if clearNumericalSolutions_ is set to true.
             */
            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, NumberOfStateColumns > > equationsOfMotionNumericalSolution_;

            //! Map of state history of numerically integrated bodies.
            /*!
            *  Map of state history of numerically integrated bodies, i.e. the result of the numerical integration, in the
            *  original propagation coordinates. Key of map denotes time, values are concatenated vectors of integrated body
            * states (order defined by propagatorSettings_).
            *  NOTE: this map is empty if clearNumericalSolutions_ is set to true.
            */
            std::map <TimeType, Eigen::Matrix<StateScalarType, Eigen::Dynamic, NumberOfStateColumns>> equationsOfMotionNumericalSolutionRaw_;

            //! Map of dependent variable history that was saved during numerical propagation.
            std::map <TimeType, Eigen::VectorXd> dependentVariableHistory_;

            //! Map of cumulative computation time history that was saved during numerical propagation.
            std::map<TimeType, double> cumulativeComputationTimeHistory_;

            //! Map of cumulative number of function evaluations that was saved during numerical propagation.
            std::map<TimeType, unsigned int> cumulativeNumberOfFunctionEvaluations_;

            //! Map listing starting entry of dependent variables in output vector, along with associated ID.
            std::map <std::pair<int, int>, std::string> dependentVariableIds_;

            std::map <std::pair<int, int>, std::string> stateIds_;

            std::shared_ptr <SingleArcPropagatorProcessingSettings> outputSettings_;

            const std::function< void ( std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&,
                                        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& ) > rawSolutionConversionFunction_;

            bool propagationIsPerformed_;

            bool solutionIsCleared_;

            //! Event that triggered the termination of the propagation
            std::shared_ptr <PropagationTerminationDetails> propagationTerminationReason_;

            friend class SingleArcDynamicsSimulator<StateScalarType, TimeType>;

            friend class MultiArcSimulationResults<StateScalarType, TimeType, NumberOfStateColumns >;


        };


        template<typename StateScalarType = double, typename TimeType = double >
        class SingleArcVariationalSimulationResults
        {
        public:
            SingleArcVariationalSimulationResults( const std::shared_ptr <SingleArcSimulationResults<StateScalarType, TimeType>> singleArcDynamicsResults,
                                                   const int stateTransitionMatrixSize,
                                                   const int sensitivityMatrixSize ):
                    singleArcDynamicsResults_( singleArcDynamicsResults ), 
                    stateTransitionMatrixSize_( stateTransitionMatrixSize ),
                    sensitivityMatrixSize_( sensitivityMatrixSize ) { }


            void reset() {
                stateTransitionSolution_.clear( );
                sensitivitySolution_.clear( );
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
            
            void splitSolution(
                    const std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >>& fullSolution,
                    std::map <TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >>& equationsOfMotionNumericalSolutionRaw )
            {
                for( auto it : fullSolution )
                {
                    stateTransitionSolution_[ static_cast< double >( it.first ) ] = it.second.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).template cast< double >( );
                    sensitivitySolution_[ static_cast< double >( it.first ) ] = it.second.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ).template cast< double >( );
                    equationsOfMotionNumericalSolutionRaw[ static_cast< double >( it.first ) ] = it.second.block( 0, stateTransitionMatrixSize_ + sensitivityMatrixSize_, stateTransitionMatrixSize_, 1 ).template cast< double >( );
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

        protected:
            const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > singleArcDynamicsResults_;

            const int stateTransitionMatrixSize_;

            const int sensitivityMatrixSize_;

            std::map < TimeType, Eigen::MatrixXd > stateTransitionSolution_;

            std::map < TimeType, Eigen::MatrixXd > sensitivitySolution_;
        };

        template<typename StateScalarType = double, typename TimeType = double, int NumberOfStateColumns = 1 >
        class MultiArcSimulationResults : public SimulationResults<StateScalarType, TimeType> {

        public:
            MultiArcSimulationResults(
                    const std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > > singleArcResults ):
                    singleArcResults_( singleArcResults ), propagationIsPerformed_( false ), solutionIsCleared_( false ){ }

            ~MultiArcSimulationResults() {}

            bool getPropagationIsPerformed( )
            {
                return propagationIsPerformed_;
            }

            void restartPropagation( )
            {
                propagationIsPerformed_ = false;
                solutionIsCleared_ = false;
                arcStartTimes_.clear( );
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    singleArcResults_.at( i )->reset( );
                }
            }

            void setPropagationIsPerformed( )
            {
                propagationIsPerformed_ = true;
                if( arcStartTimes_.size( ) != 0 )
                {
                    throw std::runtime_error( "Error, arc start times not 0 when resetting" );
                }
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    if( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolutionRaw( ).size( ) == 0 )
                    {
                        throw std::runtime_error( "Error when setting multi-arc initial times from results; results of arc " +
                            std::to_string( i ) + " are empty." );
                    }
                    arcStartTimes_.push_back( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolutionRaw( ).begin( )->first );
                }
            }

            void setPropagationIsPerformed( const std::vector< double >& arcStartTimes )
            {
                propagationIsPerformed_ = true;
                arcStartTimes_ = arcStartTimes;
            }

            void manuallySetPropagationResults(
                    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > numericalMultiArcSolution )
            {
                for( unsigned int i = 0; i < numericalMultiArcSolution.size( ); i++ )
                {
                    singleArcResults_.at( i )->equationsOfMotionNumericalSolution_ = numericalMultiArcSolution.at( i );
                }
                propagationIsPerformed_ = true;
                for( unsigned int i = 0; i < singleArcResults_.size( ); i++ )
                {
                    arcStartTimes_.push_back( singleArcResults_.at( i )->getEquationsOfMotionNumericalSolution( ).begin( )->first );
                }
            }

            std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > > getSingleArcResults( )
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

        private:
            const std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > > singleArcResults_;

            bool propagationIsPerformed_;

            bool solutionIsCleared_;

            //! List of start times of each arc. NOTE: This list is updated after every propagation.
            std::vector< double > arcStartTimes_;

        };

//        template<typename StateScalarType = double, typename TimeType >
//        std::shared_ptr< MultiArcSimulationResults<StateScalarType, TimeType, Eigen::Dynamic > > createVariationalSimulationResults(
//                const std::shared_ptr< MultiArcSimulationResults<StateScalarType, TimeType, 1 > > simulationResults )
//        {
//            std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, 1 > > > singleArcResults =
//                    simulationResults->getSingleArcResults( );
//            std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, Eigen::Dynamic > > > singleArcVariationalResults;
//            for( unsigned int i = 0; i < singleArcResults.size( ); i++ )
//            {
//                singleArcVariationalResults.push_back( createVariationalSimulationResults( singleArcResults.at( i ) ) );
//            }
//
//            return std::make_shared< MultiArcSimulationResults<StateScalarType, TimeType, Eigen::Dynamic > >(
//                    singleArcVariationalResults );
//        }

        template<typename StateScalarType = double, typename TimeType = double, int NumberOfStateColumns = 1 >
        class HybridArcSimulationResults : public SimulationResults<StateScalarType, TimeType>
        {
        public:
            HybridArcSimulationResults(
                    const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > singleArcResults,
                    const std::shared_ptr< MultiArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > multiArcResults ):
                    singleArcResults_( singleArcResults ), multiArcResults_( multiArcResults ){ }

            ~HybridArcSimulationResults() {}

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
                std::vector< std::map< TimeType, Eigen::VectorXd > > concatenatedResults;
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
                    std::vector< std::map< TimeType, Eigen::VectorXd > >
                            multiArcResults = multiArcResults_->getConcatenatedEquationsOfMotionResults( );
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
            
        protected:
            std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > singleArcResults_;

            std::shared_ptr< MultiArcSimulationResults< StateScalarType, TimeType, NumberOfStateColumns > > multiArcResults_;
        };

    }

}
#endif //TUDAT_PROPAGATIONRESULTS_H
