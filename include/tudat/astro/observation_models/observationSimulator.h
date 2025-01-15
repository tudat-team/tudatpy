/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONSIMULATOR_H
#define TUDAT_OBSERVATIONSIMULATOR_H

#include "tudat/basics/utilities.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/observationViabilityCalculator.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{

//! Virtual base class for the observation simulator class.
/*!
 *  Virtual base class for the observation simulator class, which is used to compute observable values of a
 *  single type of observable This base class is used for practical purposes, as the derived class has a
 *  template argument for the observable size, precluding the possibility of making a list of objects for all observation
 *  simulators (e.g. one for each observable type)
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class ObservationSimulatorBase
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Type of observable for which this object computes observations.
     */
    ObservationSimulatorBase(
            const ObservableType observableType ):
        observableType_( observableType ){ }

    //! Destructor
    virtual ~ObservationSimulatorBase( ){ }

    //! Function to get the type of observable for which this object computes observations
    /*!
     * Function to get the type of observable for which this object computes observations
     * \return Type of observable for which this object computes observations
     */
    ObservableType getObservableType( )
    {
        return observableType_;
    }

    //! Function to get the observation model for a given set of link ends
    /*!
     * Function to get the observation model for a given set of link ends
     * \return Observation model for a given set of link ends
     */
    virtual int getObservationSize( ) = 0;

    virtual void computeObservations( const std::vector< TimeType >& times,
                              const LinkEnds linkEnds,
                              const LinkEndType linkEndAssociatedWithTime,
                              const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings,
                              Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector ) = 0;

protected:

    //! Type of observable for which this object computes observations
    ObservableType observableType_;

};

//! Objects used to simulate a set of observations of a given kind
template< int ObservationSize, typename ObservationScalarType = double, typename TimeType = double >
class ObservationSimulator: public ObservationSimulatorBase< ObservationScalarType, TimeType >
{
public:

    using ObservationSimulatorBase< ObservationScalarType, TimeType >::observableType_;
    //! Constructor
    /*!
     * Constructor
     * \param observableType Type of observable for which this object computes observations
     * \param observationModels List of observation models of type observableType
     */
    ObservationSimulator(
            const ObservableType observableType,
            const std::map< LinkEnds, std::shared_ptr< ObservationModel< ObservationSize,
            ObservationScalarType, TimeType > > >& observationModels ):
        ObservationSimulatorBase< ObservationScalarType, TimeType >( observableType ), observationModels_( observationModels ){ }

    //! Virtual destructor
    virtual ~ObservationSimulator( ){ }

    //! Function to get the size of the observable for a given set of link ends
    /*!
     * Function to get the size of the observable for a given set of link ends
     * \return Size of the observable for a given set of link ends
     */
    int getObservationSize( )
    {
        return observationModels_.begin( )->second->getObservationSize( );
    }

    //! Function to get the observation model for a given set of link ends
    /*!
     * Function to get the observation model for a given set of link ends
     * \return Observation model for a given set of link ends
     */
    std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > getObservationModel(
            const LinkEnds linkEnds )
    {
        if( observationModels_.count( linkEnds ) == 0 )
        {
            throw std::runtime_error( "Error in observation manager when getting observation model, did not find model for given link ends: " + getLinkEndsString( linkEnds ) );
        }
        else
        {
            return observationModels_.at( linkEnds );
        }
    }

    //! Function to get the full list of observation models
    /*!
     * Function to get the full list of observation models
     * \return Full list of observation models
     */
    std::map< LinkEnds, std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
    getObservationModels( )
    {
        return observationModels_;
    }

    void computeObservations( const std::vector< TimeType >& times,
                                          const LinkEnds linkEnds,
                                          const LinkEndType linkEndAssociatedWithTime,
                                          const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings,
                                          Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector )
    {
        // Initialize return vectors.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > > observations;

        // Get observation model.
        std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > selectedObservationModel =
            this->getObservationModel( linkEnds );

        // Initialize vectors of states and times of link ends to be used in calculations.
        std::vector< Eigen::Vector6d > vectorOfStates;
        std::vector< double > vectorOfTimes;

        Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation;

        // Iterate over all observation times
//        int currentObservationSize;
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
//            currentObservationSize = currentObservation.rows( );
            observations[ saveTime ] = currentObservation;
        }

        observationsVector = utilities::createConcatenatedEigenMatrixFromMapValues<TimeType, ObservationScalarType, ObservationSize, 1>( observations );
    }

protected:

    //! List of observation models of type observableType
    std::map< LinkEnds, std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
    observationModels_;
};

template< int ObservationSize, typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >
getObservationSimulatorOfType(
        const std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators,
        const ObservableType observableType )
{
    std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator;
    for( unsigned int i = 0; i < observationSimulators.size( ); i++ )
    {
        if( observationSimulators.at( i )->getObservableType( ) == observableType )
        {
            if( observationSimulator != nullptr )
            {
                throw std::runtime_error( "Error when getting observation simulator of single type; multiple simulators detected" );
            }

            observationSimulator = std::dynamic_pointer_cast<  ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >(
                        observationSimulators.at( i ) );

            if( ObservationSize != observationSimulators.at( i )->getObservationSize( ) )
            {
                std::cout<<observableType<<" "<<ObservationSize<<" "<<observationSimulators.at( i )->getObservationSize( )<<std::endl;
                throw std::runtime_error( "Error when getting observation simulator of single type; sizes are incompatible" );
            }

            if( observationSimulator == nullptr )
            {
                throw std::runtime_error( "Error when getting observation simulator of single type; dynamic cast unsuccesful" );
            }
        }
    }
    if( observationSimulator == nullptr )
    {
        throw std::runtime_error( "Error when retrieving observation simulator from list for type " +
                                  std::to_string( observableType ) + ", no simualtor found" );
    }
    return observationSimulator;
}


}

}
#endif // TUDAT_OBSERVATIONSIMULATOR_H
