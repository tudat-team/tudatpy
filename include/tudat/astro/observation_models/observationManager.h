/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONMANAGER_H
#define TUDAT_OBSERVATIONMANAGER_H

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/propagators/stateTransitionMatrixInterface.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

namespace tudat
{

namespace observation_models
{

//! Virtual base class for the observation manager class,.
/*!
 *  Virtual base class for the observation manager class, which is used to compute observable values and partials of a
 *  single type of observable during estimation. This base class is sued for practical purposes, as  the derived class has a\
 *  template argument for the observable size, precluding the possibility of making a list of objects for all observation
 *  managers (e.g. one for each observable type)
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class ObservationManagerBase
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Type of observable for which the instance of this class will compute observations/observation
     * partials
     * \param stateTransitionMatrixInterface Object used to compute the state transition/sensitivity matrix at a given time
     * \param observationPartialScalers Map of objects (one per set of link ends) used to compute the scaling of
     * position partials that are used to compute the observation partials in the derived class
     */
    ObservationManagerBase(
            const ObservableType observableType,
            const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
            stateTransitionMatrixInterface,
            const std::map< LinkEnds, std::shared_ptr< observation_partials::PositionPartialScaling  > >&
            observationPartialScalers ):
        observableType_( observableType ), stateTransitionMatrixInterface_( stateTransitionMatrixInterface ),
        observationPartialScalers_( observationPartialScalers )
    {
        if( stateTransitionMatrixInterface_ != nullptr )
        {
            stateTransitionMatrixSize_ = stateTransitionMatrixInterface_->getStateTransitionMatrixSize( );
//            std::cout << "in observation manager, STM size: " << stateTransitionMatrixSize_ << "\n\n";
        }
        else
        {
            stateTransitionMatrixSize_ = 0;
        }
    }

    //! Virtual destructor
    virtual ~ObservationManagerBase( ){ }

    //! Pure virtual function to return the size of the observable for a given set of link ends
    /*!
     * Pure virtual function to return the size of the observable for a given set of link ends
     * \param linkEnds Link ends for which observable size is to be computed
     * \return Size of observable
     */
    virtual int getObservationSize( ) = 0;


    //! Function to simulate observations between specified link ends and associated partials at set of observation times.
    /*!
     *  Function(pure virtual)   to simulate observations between specified link ends  and associated partials at set of
     *  observation times, used the sensitivity and state transition matrix interpolators set in this base class.
     *  \param times Vector of times at which observations are performed
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param linkEndAssociatedWithTime Link end at which input times are valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \return Pair of observable values and partial matrix
     */
    virtual void computeObservationsWithPartials( const std::vector< TimeType >& times,
                                          const LinkEnds linkEnds,
                                          const LinkEndType linkEndAssociatedWithTime,
                                          const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings,
                                          Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector,
                                          Eigen::MatrixXd& partialsMatrix,
                                          const bool calculateObservations = true,
                                          const bool calculatePartials = true ) = 0;
    //! Function (ṕure virtual) to return the object used to simulate noise-free observations
    /*!
     * Function (ṕure virtual) to return the object used to simulate noise-free observations
     * \return Object used to simulate ideal observations
     */
    virtual std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > getObservationSimulator( ) = 0;


protected:

    //! Function to get the state transition and sensitivity matrix.
    /*!
     *  Function to get the state transition matrix Phi and sensitivity matrix S at a given time as a single matrix [Phi;S]
     *  \param evaluationTime Time at which matrices are to be evaluated
     *  \return Concatenated state transition and sensitivity matrices at given time.
     */
    Eigen::MatrixXd getCombinedStateTransitionAndSensitivityMatrix( const double evaluationTime,
                                                                    const std::vector< std::string >& arcDefiningBodies = std::vector< std::string >( ) )
    {
        return stateTransitionMatrixInterface_->getFullCombinedStateTransitionAndSensitivityMatrix( evaluationTime, true, arcDefiningBodies );
    }


    //! Type of observable for which the instance of this class will compute observations/observation partials
    ObservableType observableType_;

    //! Object used to compute the state transition/sensitivity matrix at a given time
    std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface_;

    //!  Map of objects (one per set of link ends) used to compute the scaling of position partials that are used to
    //! compute the observation partials in the derived class
    std::map< LinkEnds, std::shared_ptr< observation_partials::PositionPartialScaling  > > observationPartialScalers_;

    //! Size of (square) state transition matrix.
    /*!
     *  Size of (square) state transition matrix.
     */
    int stateTransitionMatrixSize_;

};

//! Class to manage simulation of observables and associated partials for a single type of observable.
/*!
 *  This class manages the simulation of observables and their partials, which are used during estimation.
 *  Separate observation models of a single kind (i.e. between different link ends) are all handled by a single object
 *  of this type.
 */
template< int ObservationSize, typename ObservationScalarType = double, typename TimeType = double >
class ObservationManager: public ObservationManagerBase< ObservationScalarType, TimeType >
{
public:

    //! Typedef for translational state type
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    // Using statements
    using ObservationManagerBase< ObservationScalarType, TimeType >::stateTransitionMatrixSize_;
    using ObservationManagerBase< ObservationScalarType, TimeType >::observationPartialScalers_;
    using ObservationManagerBase< ObservationScalarType, TimeType >::stateTransitionMatrixInterface_;

    //! Constructor
    /*!
     * Constructor
     * \param observableType Type of observable for which the instance of this class will compute observations/observation
     * partials
     * \param observationSimulator Object used to simulate ideal observations of the  observableType
     * \param observationPartials List of objects that compute the partial derivatives of the observations w.r.t. the
     * estimated parameters. Each partial (i.e. w.r.t. each parameter) has its own associated ObservationPartial object.
     * The LinkEnds key denotes the specific geomtry of the observable, while the pair< int, int > key denotes the start
     * index and size, respectively, of the current parameter partial in the estimated parameter vector.
     * \param observationPartialScalers Map of objects (one per set of link ends) used to compute the scaling of
     * position partials that are used to compute the observation partials in the derived class
     * \param stateTransitionMatrixInterface Object used to compute the state transition/sensitivity matrix at a given time
     */
    ObservationManager(
            const ObservableType observableType,
            const std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >&
            observationSimulator,
            const std::map< LinkEnds, std::map< std::pair< int, int >,
            std::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > > observationPartials,
            const std::map< LinkEnds, std::shared_ptr< observation_partials::PositionPartialScaling  > >
            observationPartialScalers,
            const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
            stateTransitionMatrixInterface,
            const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) ):
        ObservationManagerBase< ObservationScalarType, TimeType >(
            observableType, stateTransitionMatrixInterface, observationPartialScalers ),
        observationSimulator_( observationSimulator ), observationPartials_( observationPartials ), dependentVariablesInterface_( dependentVariablesInterface ){ }

    //! Virtual destructor
    virtual ~ObservationManager( ){ }

    //! Function to return the size of the observable for a given set of link ends
    /*!
     * Function to return the size of the observable for a given set of link ends
     * \param linkEnds Link ends for which observable size is to be computed
     * \return Size of observable
     */
    int getObservationSize( )
    {
        return observationSimulator_->getObservationSize( );
    }

    //! Function to return the observation model for a given set of link ends
    /*!
     *  Function to return the observation model for a given set of link ends
     *  \param linkEnds Link ends for which observation model is to be returned
     *  \return Observation model
     */
    std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > getObservationModel(
            const LinkEnds linkEnds )
    {
        return observationSimulator_->getObservationModel( linkEnds );
    }

    //! Function to return the object used to simulate noise-free observations
    /*!
     * Function to return the object used to simulate noise-free observations
     * \return Object used to simulate ideal observations
     */
    std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > getObservationSimulator( )
    {
        return observationSimulator_;
    }

    //! Function to simulate observations between specified link ends and associated partials at set of observation times.
    /*!
     *  Function to simulate observations between specified link ends  and associated partials at set of observation times,
     *  used the sensitivity and state transition matrix interpolators set in the base class.
     *  \param times Vector of times at which observations are performed
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param linkEndAssociatedWithTime Link end at which input times are valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \return Pair of observable values and partial matrix
     */
    void computeObservationsWithPartials( const std::vector< TimeType >& times,
                                     const LinkEnds linkEnds,
                                     const LinkEndType linkEndAssociatedWithTime,
                                     const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings,
                                     Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector,
                                     Eigen::MatrixXd& partialsMatrix,
                                     const bool calculateObservations = true,
                                     const bool calculatePartials = true )
    {
        // Initialize return vectors.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > > observations;
        std::map< TimeType, Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > > partialsMatrices;

        // Get observation model.
        std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > selectedObservationModel =
                observationSimulator_->getObservationModel( linkEnds );

        // Initialize vectors of states and times of link ends to be used in calculations.
        std::vector< Eigen::Vector6d > vectorOfStates;
        std::vector< double > vectorOfTimes;

        Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation;

        // Iterate over all observation times
        int currentObservationSize;
        for( unsigned int i = 0; i < times.size( ); i++ )
        {
            vectorOfTimes.clear( );
            vectorOfStates.clear( );

            // Compute observation
            currentObservation = selectedObservationModel->computeObservationsWithLinkEndData(
                        times[ i ], linkEndAssociatedWithTime, vectorOfTimes, vectorOfStates, ancilliarySettings );
            TimeType saveTime = times[ i ];
            while( observations.count( saveTime ) != 0 )
            {
                saveTime += std::numeric_limits< double >::epsilon( ) * 10.0 * times[ i ];
            }



            // Compute observation partial
            currentObservationSize = currentObservation.rows( );

            if( calculateObservations )
            {
                observations[ saveTime ] = currentObservation;
            }

            if( calculatePartials )
            {
                partialsMatrices[ saveTime ] = determineObservationPartialMatrix(
                    currentObservationSize, vectorOfStates, vectorOfTimes, linkEnds, currentObservation,
                    linkEndAssociatedWithTime, ancilliarySettings );
            }
        }

        if( calculateObservations )
        {
            observationsVector = utilities::createConcatenatedEigenMatrixFromMapValues<TimeType, ObservationScalarType, ObservationSize, 1>( observations );
        }

        if( calculatePartials )
        {
            partialsMatrix = utilities::createConcatenatedEigenMatrixFromMapValues<TimeType, double, ObservationSize, Eigen::Dynamic>( partialsMatrices );
        }
    }

    //! Function to return the full list of observation partial objects
    /*!
     * Function to return the full list of observation partial objects
     * \return Full list of observation partial objects
     */
    std::map< LinkEnds, std::map< std::pair< int, int >, std::shared_ptr<
    observation_partials::ObservationPartial< ObservationSize > > > > getObservationPartials( )
    {
        return observationPartials_;
    }

    //! Function to return the observation partial objects for a single set of link ends
    /*!
     * Function to return the observation partial objects for a single set of link ends
     * \return Observation partial objects for a single set of link ends
     */
    std::map< std::pair< int, int >, std::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > >
    getObservationPartials( const LinkEnds& linkEnds )
    {
        return observationPartials_.at( linkEnds );
    }


protected:

    //! Function to perform updates of dependent variables used by (subset of) observation partials.
    /*!
     *  Function to perform updates of dependent variables used by (subset of) observation partials, in order
     *  to decrease redundant computations (i.e. scaling of position partials for range partials).
     *  is kept constant (to input value)
     *  \param states List of times at each link end during observation
     *  \param times List of states at each link end during observation
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    virtual void updatePartials(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const LinkEnds& linkEnds,
            const LinkEndType linkEndAssociatedWithTime,
            const Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation)
    {
        if( observationPartialScalers_.at( linkEnds ) != nullptr )
        {
            observationPartialScalers_.at( linkEnds )->update( states, times, linkEndAssociatedWithTime,
                                                               currentObservation.template cast< double >( ) );
        }
    }

    //! Function to calculate range partials at given states between link ends and reception and transmission time.
    /*!
     *  Function to calculate range partials at given states of link ends and reception and transmission times.
     *  \param observationSize Size of single observation
     *  \param states States of link ends, order determined by updatePartials( )
     *  and calculatePartial( ) functions expected inputs.
     *  \param times Times at link ends (reception, transmission, reflection, etc. ), order determined by updatePartials( )
     *  and calculatePartial( ) functions expected inputs.
     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
     *  \param currentObservation Value of observation for which partials are to be computed
     *  \param linkEndAssociatedWithTime Reference link end for observations
     *  \return Matrix of partial derivative of observation w.r.t. parameter vector.
     */
    Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > determineObservationPartialMatrix(
            const int observationSize,
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const LinkEnds& linkEnds,
            const Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation,
            const LinkEndType linkEndAssociatedWithTime,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings )
    {
        // Initialize partial vector of observation w.r.t. all parameter.
        int fullParameterVector = stateTransitionMatrixInterface_->getFullParameterVectorSize( );

        Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > partialMatrix =
                Eigen::MatrixXd::Zero( observationSize, fullParameterVector );

        // Initialize list of [Phi;S] matrices at times required by calculation (key)
        std::map< double, Eigen::MatrixXd > combinedStateTransitionMatrices;

        // Perform updates of dependent variables used by (subset of) observation partials.
        updatePartials( states, times, linkEnds, linkEndAssociatedWithTime, currentObservation );

        currentLinkEndPartials = observationPartials_[ linkEnds ];

        // Get list of bodies involved in linkEnds
        std::vector< std::string > bodiesInLinkEnds;
        for ( auto itr : linkEnds )
        {
            if ( std::count( bodiesInLinkEnds.begin( ), bodiesInLinkEnds.end( ), itr.second.bodyName_  ) == 0 )
            {
                bodiesInLinkEnds.push_back( itr.second.bodyName_ );
            }
        }

        // Iterate over all observation partials associated with given link ends.
        for( typename std::map< std::pair< int, int >, std::shared_ptr<
             observation_partials::ObservationPartial< ObservationSize > > >::iterator
             partialIterator = currentLinkEndPartials.begin( );
             partialIterator != currentLinkEndPartials.end( ); partialIterator++ )
        {
            // Get Observation partial start and size indices in parameter vector.
            std::pair< int, int > currentIndexInfo = partialIterator->first;

            std::vector< std::string > bodiesOfInterestInLinkEnds;
            for ( unsigned int k = 0  ; k < bodiesInLinkEnds.size( ) ; k++ )
            {
                if ( partialIterator->second->getParameterIdentifier( ).second.first == bodiesInLinkEnds.at( k ) )
                {
                    bodiesOfInterestInLinkEnds.push_back( bodiesInLinkEnds.at( k ) );
                }
            }
            if ( bodiesOfInterestInLinkEnds.size( ) == 0 )
            {
                bodiesOfInterestInLinkEnds = bodiesInLinkEnds;
            }

            // Calculate partials of observation w.r.t. parameters, with associated observation times (single partial
            // can consist of multiple partial matrices, associated at different times)
            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > singlePartialSet =
                    partialIterator->second->calculatePartial(
                            states, times, linkEndAssociatedWithTime, ancilliarySettings, currentObservation.template cast< double >( ) );

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
                                this->getCombinedStateTransitionAndSensitivityMatrix( singlePartialSet[ i ].second, bodiesOfInterestInLinkEnds /*bodiesInLinkEnds*/ );
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
                    // Partial w.r.t. observation time property
                    if ( isParameterObservationLinkTimeProperty( partialIterator->second->getParameterIdentifier( ).first ) )
                    {
                        // Iterate (again) over all observation partials to retrieve those associated with given link ends states.
                        for( auto itr : currentLinkEndPartials )
                        {
                            // Get observation partial start and size indices in parameter vector.
                            std::pair< int, int > indexInfo = itr.first;

                            if( indexInfo.first < stateTransitionMatrixSize_ )
                            {
                                // Calculate partials of observation w.r.t. link end states, with associated observation times
                                // (single partial can consist of multiple partial matrices, associated at different times)
                                std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > linkEndStatePartialSet =
                                        itr.second->calculatePartial(
                                                states, times, linkEndAssociatedWithTime, ancilliarySettings,
                                                currentObservation.template cast< double >( ) );

                                for( unsigned int j = 0; j < linkEndStatePartialSet.size( ); j++ )
                                {
                                    std::string nameBody = itr.second->getParameterIdentifier( ).second.first;
                                    int indexLinkEndType = -1;
                                    for ( auto itrLinkEnds : linkEnds )
                                    {
                                        if ( itrLinkEnds.second.bodyName_ == nameBody )
                                        {
                                            indexLinkEndType =
                                                static_cast< int >( getLinkEndIndicesForLinkEndTypeAtObservable( this->observableType_, itrLinkEnds.first, linkEnds.size( ) ).at( 0 ) );
                                        }
                                    }

                                    if( indexLinkEndType < 0 )
                                    {
                                        throw std::runtime_error( "Error when retrieving partial for time bias; could not find associated link end index.");
                                    }
                                    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > totalAccelerationVariable
                                            = std::make_shared< propagators::SingleDependentVariableSaveSettings >( propagators::total_acceleration_dependent_variable, nameBody );

                                    if( dependentVariablesInterface_ == nullptr )
                                    {
                                        throw std::runtime_error( "Error, required dependent variable interfaces, but none found when computing partials." );
                                    }

                                    Eigen::VectorXd acceleration = dependentVariablesInterface_->getSingleDependentVariable(
                                        totalAccelerationVariable, times.at( indexLinkEndType ) );
                                    if( acceleration.rows( ) == 0 )
                                    {
                                        throw std::runtime_error( "Error when getting link time property partial; could not find acceleration of " + nameBody );
                                    }
                                    Eigen::Vector6d stateDerivativeVector = Eigen::Vector6d::Zero( );
                                    stateDerivativeVector.segment( 0, 3 ) = states.at( indexLinkEndType ).segment( 3, 3 );
                                    stateDerivativeVector.segment( 3, 3 ) = acceleration;

                                    Eigen::MatrixXd partialWrtStateDerivative = ( linkEndStatePartialSet[ j ].first ) * stateDerivativeVector;
                                    partialMatrix.block( 0, currentIndexInfo.first, observationSize, currentIndexInfo.second ) +=
                                            partialWrtStateDerivative * singlePartialSet[ i ].first;
                                }
                            }
                        }
                    }

                    // Partial w.r.t. observation bias
                    else
                    {
                        // Add direct partial of observation w.r.t. parameter.
                        partialMatrix.block( 0, currentIndexInfo.first, observationSize, currentIndexInfo.second ) +=
                                singlePartialSet[ i ].first;
                    }
                }
            }
        }
        return partialMatrix;
    }

    //! Object used to simulate ideal observations of the  observableType
    std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator_;

    //! Map of observation partials.
    /*!
     * List of objects that compute the partial derivatives of the observations w.r.t. the
     * estimated parameters. Each partial (i.e. w.r.t. each parameter) has its own associated ObservationPartial object.
     * The LinkEnds key denotes the specific geomtry of the observable, while the pair< int, int > key denotes the start
     * index and size, respectively, of the current parameter partial in the estimated parameter vector.
     */
    std::map< LinkEnds, std::map< std::pair< int, int >, std::shared_ptr<
    observation_partials::ObservationPartial< ObservationSize > > > > observationPartials_;

    //! Pre-declared map used in computation of partials.
    std::map< std::pair< int, int >, std::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > >
    currentLinkEndPartials;

    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface_;

};


extern template class ObservationManagerBase< double, double >;
extern template class ObservationManager< 1, double, double >;
extern template class ObservationManager< 2, double, double >;
extern template class ObservationManager< 3, double, double >;
extern template class ObservationManager< 6, double, double >;

}

}

#endif // TUDAT_OBSERVATIONMANAGER_H
