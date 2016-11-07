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

#include <iostream>
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
 *  Each type of observables (1-way range, 2-way range, Doppler, VLBI, etc.) has its own
 *  derived class capable of simulating observables of the given type using given link ends.
 *  The functions to be used for computing the observables can be called with/without deviations from ideal observable
 *  (see base class member functions). Corrections are computed from an observationBiasCalculator member object, which is
 *  empty by default. Also, the observable may be a with/without returning (by reference) the times and states
 *  at each of the link ends. Returning these times/states prevents recomputations of these quantities in later calculations.
 */
template< int ObservationSize = Eigen::Dynamic, typename ObservationScalarType = double, typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
class ObservationModel
{
public:

    //! Constructor
    /*!
     * Base class constructor.
     * \param observableType Type of observable, used for derived class type identification without
     * explicit casts.
     * \param observationBiasCalculator Object for calculating system-dependent errors in the
     * observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    ObservationModel(
            const ObservableType observableType ,
            const boost::shared_ptr< ObservationBias< ObservationSize > > observationBiasCalculator = NULL ):
        observableType_( observableType ),
        observationBiasCalculator_( observationBiasCalculator )
    {
        // Check if bias is empty
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

    //! Function to compute the observable without any corrections
    /*!
     * Function to compute the observable without any corrections, i.e. the ideal physical observable as computed
     *  from the defined link ends (in the derived class). Note that this observable does include e.g. light-time
     *  corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Ideal observable.
     */
    virtual Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeIdealObservationsWithLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes,
                std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) = 0;

    //! Function to compute full observation at given time.
    /*!
     *  Function to compute observation at given time (include any defined non-ideal corrections). The
     *  times and states of the link ends are given in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservationsWithLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes ,
                std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates )
    {
        // Check if any non-ideal models are set.
        if( isBiasNull_ )
        {
            return computeIdealObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
        }
        else
        {
            // Compute ideal observable
            Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation =
                    computeIdealObservationsWithLinkEndData(
                                            time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );

            // Add correction
            return currentObservation +
                    this->observationBiasCalculator_->getObservationBias( linkEndTimes, linkEndStates ).
                    template cast< ObservationScalarType >( );
        }
    }

    //! Function to compute the observable without any corrections.
    /*!
     * Function to compute the observable without any corrections, i.e. the ideal physical observable as computed
     * from the defined link ends (in the derived class). Note that this observable does include e.g. light-time
     * corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     * errors, such as biases or clock errors. This function may be redefined in derived class for improved efficiency.
     * \param time Time at which observable is to be evaluated.
     * \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     * is kept constant (to input value)
     * \return Ideal observable.
     */
    virtual Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )
    {
        // Compute ideal observable from derived class.
        return this->computeIdealObservationsWithLinkEndData(
                    time, linkEndAssociatedWithTime, this->linkEndTimes_, this->linkEndStates_ );
    }

    //! Function to compute full observation at given time.
    /*!
     *  Function to compute observation at given time (include any defined non-ideal corrections).
     * \param time Time at which observable is to be evaluated.
     * \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     * is kept constant (to input value)
     *  \return Calculated (non-ideal) observable value.
     */
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )
    {
        // Check if any non-ideal models are set.
        if( isBiasNull_ )
        {
            return computeIdealObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, linkEndTimes_, linkEndStates_ );
        }
        else
        {
            // Compute ideal observable
            Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > currentObservation =
                    computeIdealObservationsWithLinkEndData(
                                            time, linkEndAssociatedWithTime, linkEndTimes_, linkEndStates_ );

            // Add correction
            return currentObservation +
                    this->observationBiasCalculator_->getObservationBias( linkEndTimes_, linkEndStates_ ).
                    template cast< ObservationScalarType >( );
        }
    }

    ObservationScalarType computeObservationEntry(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            const int observationEntry )
    {
        if( observationEntry < ObservationSize )
        {
            return computeObservations( time, linkEndAssociatedWithTime )( observationEntry );
        }
        else
        {
            throw std::runtime_error( "Error, requesting out-of-bounds index for observation model" );
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

    //! Boolean set by constructor to denote whether observationBiasCalculator_ is NULL.
    bool isBiasNull_;


    //! Pre-define list of times used when calling function returning link-end states/times from interface function.
    std::vector< TimeType > linkEndTimes_;

    //! Pre-define list of states used when calling function returning link-end states/times from interface function.
    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > linkEndStates_;

};


} // namespace observation_models

} // namespace tudat
#endif // TUDAT_OBSERVATIONMODEL_H
