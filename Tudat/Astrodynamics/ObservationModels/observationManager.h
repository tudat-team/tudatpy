#ifndef OBSERVATIONMANAGER_H
#define OBSERVATIONMANAGER_H

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <Eigen/Sparse>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observationSimulator.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"

namespace tudat
{

namespace observation_models
{


template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class ObservationManagerBase
{
public:
    ObservationManagerBase(
            const ObservableType observableType,
            const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
            const std::map< LinkEnds, boost::shared_ptr< observation_partials::PositionPartialScaling  > >& observationPartialScalers ):
        observableType_( observableType ), stateTransitionMatrixInterface_( stateTransitionMatrixInterface ),
        observationPartialScalers_( observationPartialScalers )
    {
        if( stateTransitionMatrixInterface_ != NULL )
        {
            stateTransitionMatrixSize_ = stateTransitionMatrixInterface_->getStateTransitionMatrixSize( );
        }
        else
        {
            stateTransitionMatrixSize_ = 0;
        }
    }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~ObservationManagerBase( ){ }

    virtual int getObservationSize( const LinkEnds& linkEnds ) = 0;


    //! Function to simulate observations between specified link ends and associated partials at set of observation times.
    /*!
     *  Function to simulate observations between specified link ends  and associated partials at set of observation times,
     *  used the sensitivity and state transition matrix interpolators set in this base class.
     *  \param recetionTimes Vector of times at which observations taked place (i.e. reception time)
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \return Pair of observable values and partial matrix
     */
    virtual std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, Eigen::MatrixXd >
    computeObservationsWithPartials( const std::vector< TimeType >& times,
                                     const LinkEnds linkEnds,
                                     const LinkEndType linkEndAssociatedWithTime ) = 0;

protected:
    //! Function to get the state transition and sensitivity matrix.
    /*!
     * Function to get the state transition matrix Phi and sensitivity matrix S at a given time as a singel matrix [Phi;S]
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime )
    {
        return stateTransitionMatrixInterface_->getFullCombinedStateTransitionAndSensitivityMatrix( evaluationTime );
    }

    //! Function to perform updates of dependent variables used by (subset of) observation partials.
    /*!
     *  Function to perform updates of dependent variables used by (subset of) observation partials, in order
     *  to decrease redundant computations (i.e. scaling of position partials for range partials).
     *  \param states States of link ends, order determined by derived class
     *  and calculatePartial( ) functions expected inputs.
     *  \param times Times at link ends (reception, transmission, reflection, etc. ), order determined by derived class
     *  and calculatePartial( ) functions expected inputs.
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     */
    virtual void updatePartials(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const LinkEnds& linkEnds,
            const LinkEndType linkEndAssociatedWithTime )
    {
        observationPartialScalers_.at( linkEnds )->update( states, times, linkEndAssociatedWithTime );
    }

    ObservableType observableType_;

    boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface_;


    std::map< LinkEnds, boost::shared_ptr< observation_partials::PositionPartialScaling  > > observationPartialScalers_;

    //! Size of (square) state transition matrix.
    /*!
     *  Size of (square) state transition matrix.
     */
    int stateTransitionMatrixSize_;

};

//! Class to manage simualtion of observables and associated partials.
/*!
 *  This class manages the simulation of observables and their partials. Separate observation models of a
 *  single kind (i.e. between different link ends) are all handled by a single derived class object of this type.
 * \tparam ObservationSize Number of values in single observation (i.e. 1 for range, 2 for VLBI)
 */
template< int ObservationSize, typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class ObservationManager: public ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType >
{
public:

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    using ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType >::stateTransitionMatrixSize_;
    using ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType >::updatePartials;
    using ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType >::stateTransitionMatrixInterface_;


    ObservationManager(
            const ObservableType observableType,
            const boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > >& observationSimulator,
            const std::map< LinkEnds, std::map< std::pair< int, int >, boost::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > > observationPartials,
            const std::map< LinkEnds, boost::shared_ptr< observation_partials::PositionPartialScaling  > > observationPartialScalers,
            const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface ):
        ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType >(
            observableType, stateTransitionMatrixInterface, observationPartialScalers ),
        observationSimulator_( observationSimulator ), observationPartials_( observationPartials ){ }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~ObservationManager( ){ }

    int getObservationSize( const LinkEnds& linkEnds )
    {
        return observationSimulator_->getObservationSize( linkEnds );
    }


    boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > getObservationModel(
            const LinkEnds linkEnds )
    {
       return observationSimulator_->getObservationModel( linkEnds );
    }

    boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > getObservationSimulator( )
    {
        return observationSimulator_;
    }

    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, Eigen::MatrixXd >
    computeObservationsWithPartials( const std::vector< TimeType >& times,
                                     const LinkEnds linkEnds,
                                     const LinkEndType linkEndAssociatedWithTime )
    {
        // Initialize return vectors.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > > observations;
        std::map< TimeType, Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > > observationMatrices;

        // Get observation model.
        boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > selectedObservationModel =
                observationSimulator_->getObservationModel( linkEnds );

        // Initialize vectors of states and times of link ends to be used in calculations.
        std::vector< basic_mathematics::Vector6d > vectorOfStates;
        std::vector< double > vectorOfTimes;

        Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation;

        for( unsigned int i = 0; i < times.size( ); i++ )
        {
            vectorOfTimes.clear( );
            vectorOfStates.clear( );
            //observationSimulator_->updateObservationSettings( times[ i ], linkEnds, linkEndAssociatedWithTime );

            currentObservation = selectedObservationModel->computeObservationsWithLinkEndData(
                        times[ i ], linkEndAssociatedWithTime, vectorOfTimes, vectorOfStates );

            observations[ times[ i ] ] = currentObservation;

            int currentObservationSize = currentObservation.rows( );

            observationMatrices[ times[ i ] ] = determineObservationPartialMatrix(
                        currentObservationSize, vectorOfStates, vectorOfTimes, linkEnds, linkEndAssociatedWithTime );

        }

        return std::make_pair( utilities::createConcatenatedEigenMatrixFromMapValues( observations ),
                               utilities::createConcatenatedEigenMatrixFromMapValues( observationMatrices ) );
    }

    std::map< LinkEnds, std::map< std::pair< int, int >, boost::shared_ptr<
    observation_partials::ObservationPartial< ObservationSize > > > > getObservationPartials( )
    {
        return observationPartials_;
    }

    std::map< std::pair< int, int >, boost::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > getObservationPartials(
            const LinkEnds& linkEnds )
    {
        return observationPartials_.at( linkEnds );
    }


protected:
    //! Function to calculate range partials at given states between link ends and reception and transmission time.
    /*!
     *  Function to calculate range partials at given states of link ends and reception and transmission times.
     *  \param states States of link ends, order determined by updatePartials( )
     *  and calculatePartial( ) functions expected inputs.
     *  \param times Times at link ends (reception, transmission, reflection, etc. ), order determined by updatePartials( )
     *  and calculatePartial( ) functions expected inputs.
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \return Matrix of partial derivative of observation w.r.t. parameter vector.
     */
    Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > determineObservationPartialMatrix(
            const int observationSize,
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const LinkEnds& linkEnds,
            const LinkEndType linkEndAssociatedWithTime )
    {
        // Initialize partial vector of observation w.r.t. all parameter.
        int fullParameterVector = stateTransitionMatrixInterface_->getFullParameterVectorSize( );
        Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > partialMatrix =
                Eigen::MatrixXd::Zero( observationSize, fullParameterVector );

        // Initialize list of [Phi;S] matrices at times required by calculation (key)
        std::map< double, Eigen::MatrixXd > combinedStateTransitionMatrices;

        // Perform updates of dependent variables used by (subset of) observation partials.
        updatePartials( states, times, linkEnds, linkEndAssociatedWithTime );

        currentLinkEndPartials = observationPartials_[ linkEnds ];

        // Iterate over all observation partials associated with given link ends.
        for( typename std::map< std::pair< int, int >, boost::shared_ptr<
             observation_partials::ObservationPartial< ObservationSize > > >::iterator
             partialIterator = currentLinkEndPartials.begin( );
             partialIterator != currentLinkEndPartials.end( ); partialIterator++ )
        {
            //std::cout<<"Current partial indices: "<<partialIterator->first.first<<" "<<partialIterator->first.second<<std::endl;
            // Get Observation partial start and size indices in parameter veector.
            std::pair< int, int > currentIndexInfo = partialIterator->first;

            // Calculate partials of observation w.r.t. parameters, with associated observation times (single partial
            // can consist of multiple partial matrices, associated at different times)
            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > singlePartialSet =
                    partialIterator->second->calculatePartial( states, times, linkEndAssociatedWithTime );

            // If start index is smaller than size of state transition,
            // current partial is w.r.t. to a body to be estimated current state.
            if( currentIndexInfo.first < stateTransitionMatrixSize_ )
            {
                for( unsigned int i = 0; i < singlePartialSet.size( ); i++ )
                {
                    // Evaluate [Phi;S] matrix at each time instant associated with partial, if not yet evaluated.
                    if( combinedStateTransitionMatrices.count( singlePartialSet[ i ].second ) == 0 )
                    {
                        combinedStateTransitionMatrices[ singlePartialSet[ i ].second ] =
                                this->getCombinedStateTransitionAndSensitivityMatrix( singlePartialSet[ i ].second );
                    }

                    // Add partial of observation h w.r.t. initial state x_{0} (dh/dx_{0}=dh/dx*dx/dx_{0})
                    partialMatrix += ( singlePartialSet[ i ].first ) *
                            combinedStateTransitionMatrices[ singlePartialSet[ i ].second ].block
                            ( currentIndexInfo.first, 0, currentIndexInfo.second, fullParameterVector );
                }
            }
            else
            {
                for( unsigned int i = 0; i < singlePartialSet.size( ); i++ )
                {
                    // Add direct partial of observation w.r.t. parameter.
                    partialMatrix.block( 0, currentIndexInfo.first, observationSize, currentIndexInfo.second ) +=
                            singlePartialSet[ i ].first;
                }
            }
        }

        return partialMatrix;
    }

    boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationSimulator_;

    //! Map of observation partials.
    /*!
     *  Map of observation partials, with link ends (receiver, transmitter, reflector, etc.) as key and
     *  map as value. This map has a pair of (start index/number of indices) in parameter vector as key
     *  and a pointer to an ObservationPartial object as value.
     */
    std::map< LinkEnds, std::map< std::pair< int, int >, boost::shared_ptr<
    observation_partials::ObservationPartial< ObservationSize > > > > observationPartials_;

    std::map< std::pair< int, int >, boost::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > currentLinkEndPartials;

};

}

}

#endif // OBSERVATIONMANAGER_H
