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


//    //! Function (pure virtual) to simulate a single observation between specified link ends.
//    /*!
//     *  Function (pure virtual) to simulate a single observation between specified link ends. Function is provided to have
//     *  an interface to simulate a single observation from an ObservationModel object (in ObservationSimulator derived class)
//     *  without the need for a ObservationSize template argument.
//     *  \param observationTime Time at which observable is to be evaluated.
//     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
//     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
//     *  is kept constant (to input value)
//     *  \param linkEndTimes List of times at each link end during observation.
//     *  \param linkEndStates List of states at each link end during observation.
//     *  \return Simulated observation at requested time and link ends.
//     */
//    virtual Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > simulateObservation(
//            const TimeType observationTime,
//            const LinkEnds linkEnds,
//            const LinkEndType linkEndAssociatedWithTime,
//            std::vector< double >& linkEndTimes,
//            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates ) = 0;

//    //! Function (pure virtual) to simulate observations between specified link ends.
//    /*!
//     *  Function (pure virtual) to simulate observations between specified link ends. Users can specify whether to check for
//     *  availability of link at given reception time.
//     *  \param observationTimes Vector of times at which observations taked place (i.e. reception time)
//     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
//     *  \param linkEndAssociatedWithTime Reference link end for observable
//     *  \param checkTimes Boolean denoting whether the observation times are to be checked for viability
//     *  \return Observations at given times (concatenated in an Eigen vector), with associated times
//     */
//    virtual std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
//    simulateObservations( const std::vector< TimeType >& observationTimes,
//                          const LinkEnds linkEnds,
//                          const LinkEndType linkEndAssociatedWithTime,
//                          const bool checkTimes = true ) = 0;

//    //! Function (pure virtual)  to simulate observations between specified link ends.
//    /*!
//     *  Function (pure virtual)  to simulate observations between specified link ends. Users can specify whether to check for
//     *  availability of  link at given reception time.
//     *  \param observationTimes Vector of times at which observations taked place (i.e. reception time)
//     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
//     *  \param linkEndAssociatedWithTime Reference link end for observable
//     *  \param checkTimes Boolean denoting whether the observation times are to be checked for viability
//     *  \return Observations at given times (concatenated in an Eigen vector), with associated times and reference link end.
//     */
//    virtual std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
//    simulateObservationsWithLinkEndId( const std::vector< TimeType >& observationTimes,
//                                       const LinkEnds linkEnds,
//                                       const LinkEndType linkEndAssociatedWithTime,
//                                       const bool checkTimes = true ) = 0;

//    //! Function to set observation viability calculators
//    /*!
//     * Function to set observation viability calculators, a different list must be provided for each set of link ends.
//     * \param viabilityCalculators List of observation viability calculators, sorted by LinkEnds.
//     */
//    void setViabilityCalculators(
//            const std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > >&
//                                  viabilityCalculators )
//    {
//        viabilityCalculators_ = viabilityCalculators;
//    }

//    //! Function to retrieve observation viability calculators for a single set of link ends
//    /*!
//     *  Function to retrieve observation viability calculators for a single set of link ends
//     *  \param linkEnds Link ends for which observation viability calculators are to be retrieved
//     *  \return Observation viability calculators for requestred single set of link ends (empty vector if no calculators for
//     *  requested link ends are given in this object)
//     */
//    std::vector< std::shared_ptr< ObservationViabilityCalculator > > getLinkViabilityCalculators(
//            const LinkEnds& linkEnds )
//    {
//        if( viabilityCalculators_.count( linkEnds ) > 0 )
//        {
//            return viabilityCalculators_.at( linkEnds );
//        }
//        else
//        {
//            return std::vector< std::shared_ptr< ObservationViabilityCalculator > >( );
//        }
//    }


protected:

    //! Type of observable for which this object computes observations
    ObservableType observableType_;

    //! List of observation viability calculators, sorted by LinkEnds.
    std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > > viabilityCalculators_;
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
            throw std::runtime_error(
                        "Error in observation manager when getting observation model, did not find model for given link ends " +
                        std::to_string( observationModels_.size( ) ) );
        }
        return observationModels_.at( linkEnds );
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


//    //! Function to simulate a single observation between specified link ends.
//    /*!
//     *  Function (to simulate a single observation between specified link ends.
//     *  \param observationTime Time at which observable is to be evaluated.
//     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
//     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
//     *  is kept constant (to input value)
//     *  \param linkEndTimes List of times at each link end during observation.
//     *  \param linkEndStates List of states at each link end during observation.
//     *  \return Simulated observation at requested time and link ends.
//     */
//    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > simulateObservation(
//            const TimeType observationTime,
//            const LinkEnds linkEnds,
//            const LinkEndType linkEndAssociatedWithTime,
//            std::vector< double >& linkEndTimes,
//            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
//    {
//        if( observationModels_.count( linkEnds ) == 0 )
//        {
//            throw std::runtime_error(
//                        "Error when simulating observtion, could not find observation model for given linke ends" );
//        }

//        return  observationModels_.at( linkEnds )->computeObservationsWithLinkEndData(
//                    observationTime, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
//    }

//    //! Function to simulate observations between specified link ends.
//    /*!
//     *  Function to simulate observations between specified link ends. Users can specify whether to check for availability of
//     *  link at given reception time.
//     *  \param observationTimes Vector of times at which observations taked place (i.e. reception time)
//     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
//     *  \param linkEndAssociatedWithTime Reference link end for observable
//     *  \param checkTimes Boolean denoting whether the observation times are to be checked for viability
//     *  \return Observations at given times (concatenated in an Eigen vector), with associated times
//     */
//    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
//    simulateObservations( const std::vector< TimeType >& observationTimes,
//                          const LinkEnds linkEnds,
//                          const LinkEndType linkEndAssociatedWithTime,
//                          const bool checkTimes = true )

//    {
//        if( observationModels_.count( linkEnds ) == 0 )
//        {
//            throw std::runtime_error(
//                        "Error when simulating observtions, could not find observation model for given linke ends" );
//        }

//        std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > selectedObservationModel =
//                observationModels_.at( linkEnds );
//        std::vector< std::shared_ptr< ObservationViabilityCalculator > > currentLinkViabilityCalculators;
//        if( checkTimes == true && this->viabilityCalculators_.count( linkEnds ) > 0 )
//        {
//            currentLinkViabilityCalculators = this->viabilityCalculators_.at( linkEnds );
//        }

//        return simulateObservationsWithCheck< ObservationSize, ObservationScalarType, TimeType >(
//                    observationTimes, selectedObservationModel, linkEndAssociatedWithTime, currentLinkViabilityCalculators );
//    }

//    //! Function to simulate observations between specified link ends.
//    /*!
//     *  Function to simulate observations between specified link ends. Users can specify whether to check for availability of
//     *  link at given reception time.
//     *  \param observationTimes Vector of times at which observations taked place (i.e. reception time)
//     *  \param linkEnds Set of stations, S/C etc. in link, with specifiers of type of link end.
//     *  \param linkEndAssociatedWithTime Reference link end for observable
//     *  \param checkTimes Boolean denoting whether the observation times are to be checked for viability
//     *  \return Observations at given times (concatenated in an Eigen vector), with associated times and reference link end.
//     */
//    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >
//    simulateObservationsWithLinkEndId( const std::vector< TimeType >& observationTimes,
//                                       const LinkEnds linkEnds,
//                                       const LinkEndType linkEndAssociatedWithTime,
//                                       const bool checkTimes = true )
//    {
//        std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > simulatedObservations =
//                simulateObservations( observationTimes, linkEnds, linkEndAssociatedWithTime, checkTimes );
//        std::vector< std::shared_ptr< ObservationViabilityCalculator > > currentLinkViabilityCalculators;
//        if( checkTimes == true && this->viabilityCalculators_.count( linkEnds ) > 0 )
//        {
//            currentLinkViabilityCalculators = this->viabilityCalculators_.at( linkEnds );
//        }

//        return std::make_pair( simulatedObservations.first, std::make_pair(
//                                   simulatedObservations.second, linkEndAssociatedWithTime ) );
//    }

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
    return observationSimulator;
}

//extern template class ObservationSimulatorBase< double, double >;
//extern template class ObservationSimulator< 1, double, double >;
//extern template class ObservationSimulator< 2, double, double >;
//extern template class ObservationSimulator< 3, double, double >;
//extern template class ObservationSimulator< 6, double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class ObservationSimulatorBase< double, Time >;
extern template class ObservationSimulatorBase< long double, double >;
extern template class ObservationSimulatorBase< long double, Time >;

extern template class ObservationSimulator< 1, double, Time >;
extern template class ObservationSimulator< 1, long double, double >;
extern template class ObservationSimulator< 1, long double, Time >;

extern template class ObservationSimulator< 2, double, Time >;
extern template class ObservationSimulator< 2, long double, double >;
extern template class ObservationSimulator< 2, long double, Time >;

extern template class ObservationSimulator< 3, double, Time >;
extern template class ObservationSimulator< 3, long double, double >;
extern template class ObservationSimulator< 3, long double, Time >;

extern template class ObservationSimulator< 6, double, Time >;
extern template class ObservationSimulator< 6, long double, double >;
extern template class ObservationSimulator< 6, long double, Time >;
#endif

}

}
#endif // TUDAT_OBSERVATIONSIMULATOR_H
