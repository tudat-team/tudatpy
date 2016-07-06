/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONMODEL_H
#define TUDAT_OBSERVATIONMODEL_H

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Astrodynamics/ObservationModels/observationBias.h"

namespace tudat
{

namespace observation_models
{

//! Base class for models of observables (i.e. range, range-rate, etc.).
/*!
 *  Base class for models of observables to be used in (for instance) orbit determination.
 *  Each type of observables (1-way range, 2-way range, Doppler, VLBI, etc.) is to have its own
 *  derived class capable of simulating observables of the given type using given link ends
 */
template< int ObservationSize = Eigen::Dynamic,
          typename ObservationScalarType = double,
          typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
class ObservationModel
{
public:

    //! Constructor
    /*!
     * Base class constructor. Implementation to be done in derived class.
     * \param observableType Type of observable, used for derived class type identification without
     * explicit casts.
     * \param observationBiasCalculator Object for calculating system-dependent errors in the
     * observable, i.e. deviations from the physically true observable (default none).
     */
    ObservationModel(
            const ObservableType observableType ,
            const boost::shared_ptr< ObservationBias< ObservationSize > > observationBiasCalculator = NULL ):
        observableType_( observableType ),
        observationBiasCalculator_( observationBiasCalculator )
    {
        if( observationBiasCalculator_ != NULL )
        {
            isBiasNull_ = 0;
            if( observationBiasCalculator_->getObservationSize( ) != ObservationSize )
            {
                throw std::runtime_error( "Error when making observation model, bias size is inconsistent" );
            }
        }
        else
        {
            isBiasNull_ = 1;
        }
    }

    //! Virtual destructor
    virtual ~ObservationModel( ) { }

    //! Function to return the type of observable.
    /*!
     *  Function to return the type of observable.
     *  \return Type of observable.
     */
    ObservableType getObservableType( )
    {
        return observableType_;
    }

    virtual Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeUnbiasedObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const
    {

        std::vector< TimeType > linkEndTimes;
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > linkEndStates;

        return this->computeUnbiasedObservationsWithLinkEndData( time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );

    }

    //! Function to compute observation at given time.
    /*!
     *  This function computes the observation at a given time and should be implemented in a
     *  derived class. The time argument can given at any of the link ends.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )
    {
        if( isBiasNull_ )
        {
            return computeObservations( time, linkEndAssociatedWithTime );
        }
        else
        {
            std::vector< TimeType > linkEndTimes;
            std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > linkEndStates;

            Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation =
                    computeUnbiasedObservationsWithLinkEndData(
                                            time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
            return currentObservation +
                    this->observationBiasCalculator_->getObservationBias( linkEndTimes ).template cast< ObservationScalarType >( );
        }
    }


    virtual Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeUnbiasedObservationsWithLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes,
                std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const = 0;

    //! Function to compute observation, as well as the state and time at each link end.
    /*!
     *  Function to compute observation, as well as the state and time at each link end. The
     *  times and states of the link ends are given in full precision (determined by class template
     *  arguments). These states and times are returned by reference. This function is to be
     *  implemented for each derived class.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservationsWithLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes,
                std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates )
    {

        if( isBiasNull_ )
        {
            return computeUnbiasedObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
        }
        else
        {
            Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation =
                    computeUnbiasedObservationsWithLinkEndData(
                                            time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
            return currentObservation +
                    this->observationBiasCalculator_->getObservationBias( linkEndTimes ).template cast< ObservationScalarType >( );
        }
    }

    //! Function to return the size of the observable
    /*!
     *  Function to return the size of the observable
     *  \return Size of the observable
     */
    int getObservationSize( )
    {
        return ObservationSize;
    }

protected:

    //! Type of observable, used for derived class type identification without explicit casts.
    ObservableType observableType_;

    //! Object for calculating system-dependent errors in the observable.
    /*!
     *  Object for calculating system-dependent errors in the observable, i.e. deviations from the
     *  physically true observable
     */
    boost::shared_ptr< ObservationBias< ObservationSize > > observationBiasCalculator_;

    bool isBiasNull_;

};


}

}
#endif // TUDAT_OBSERVATIONMODEL_H
