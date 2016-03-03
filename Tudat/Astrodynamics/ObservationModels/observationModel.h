#ifndef OBSERVATIONMODEL_H
#define OBSERVATIONMODEL_H

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
            const boost::shared_ptr< ObservationBiasInterface > observationBiasCalculator = NULL ):
        observableType_( observableType ),
        observationBiasCalculator_( observationBiasCalculator )
    {
        if( observationBiasCalculator_ != NULL )
        {
            if( observationBiasCalculator_->getObservationSize( ) != ObservationSize )
            {
                throw std::runtime_error( "Error when making observation model, bias size is inconsistent" );
            }
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

    //! Function to compute observation at given time.
    /*!
     *  This function computes the observation at a given time and should be implemented in a
     *  derived class. The time argument can given at any of the link ends.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \return Calculated observable value.
     */
    virtual Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const = 0;

    //! Function to compute observation, as well as the state and time at each link end.
    /*!
     *  Function to compute observation, as well as the state and time at each link end. The
     *  times and states of the link ends are given in double precision. They are returned by
     *  reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservationsAndLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< basic_mathematics::Vector6d >& linkEndStates ) const
    {
        // Create variables for call to full accuracy function.
        std::vector< TimeType > fullLinkEndTimes;
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > fullLinkEndStates;

        // Call function computing observation and full precion link end state/time.
        Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > observation =
                computeObservationsAndFullPrecisionLinkEndData(
                    time, linkEndAssociatedWithTime, fullLinkEndTimes, fullLinkEndStates );

        // Cast states and times to double precision
        linkEndTimes = utilities::staticCastSVectorToTVector< TimeType, double >( fullLinkEndTimes );
        linkEndStates = utilities::staticCastSEigenTypeVectorToTEigenTypeVector<
                Eigen::Matrix< StateScalarType, 6, 1 >, basic_mathematics::Vector6d, double >( fullLinkEndStates );

        // Return observable.
        return observation;
    }

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
    virtual Eigen::Matrix< ObservationScalarType, ObservationSize, 1 > computeObservationsAndFullPrecisionLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes,
                std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const = 0;

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
    boost::shared_ptr< ObservationBiasInterface > observationBiasCalculator_;

};

//! Template specialization for ObservationModel for observables of size 1.
template< typename ObservationScalarType, typename TimeType, typename StateScalarType >
class ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType >
{
public:

    //! Constructor
    /*!
     * Base class constructor for observtion model of size 1. Implementation to be done in derived
     * class.
     * \param observableType Type of observable, used for derived class type identification without
     * explicit casts.
     * \param observationBiasCalculator Object for calculating system-dependent errors in the
     * observable, i.e. deviations from the physically true observable
     */
    ObservationModel( const ObservableType observableType, const boost::shared_ptr< ObservationBiasInterface > observationBiasCalculator =
            boost::make_shared< ObservationBiasInterface >( 1 ) ):
        observableType_( observableType ),
        observationBiasCalculator_( observationBiasCalculator ){ }

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

    //! Function to compute observation at given time.
    /*!
     *  This function computes the observation at a given time and should be implemented in a
     *  derived class. The time argument can given at any of the link ends. This function
     *  returns a scalar, as opposed to the general vector return type of the unspecialized
     *  ObservationModel computeObservations function.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \return Calculated observable value.
     */
    virtual ObservationScalarType computeObservation(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const = 0 ;

    //! Function to compute observation at given time.
    /*!
     *  This function computes the observation at a given time. It uses the computeObservation
     *  function to calculate the observable. The time argument can given at any of the link ends.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const
    {
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( )<<
                 computeObservation( time, linkEndAssociatedWithTime ) ).finished( );
    }

    //! Function to compute observation, as well as the state and time at each link end.
    /*!
     *  Function to compute observation, as well as the state and time at each link end. The
     *  times and states of the link ends are given in full precision (determined by class template
     *  arguments). These states and times are returned by reference. This function is to be
     *  implemented for each derived class. This function returns a scalar, as opposed to the
     *  general vector return type of the unspecialized ObservationModel
     *  computeObservationsAndFullPrecisionLinkEndData function.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Calculated observable value.
     */
    virtual ObservationScalarType computeObservationAndFullPrecisionLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< TimeType >& linkEndTimes,
            std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const = 0;

    //! Function to compute observation, as well as the state and time at each link end.
    /*!
     *  Function to compute observation, as well as the state and time at each link end. The
     *  times and states of the link ends are given in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeObservationsAndFullPrecisionLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< TimeType >& linkEndTimes,
                    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const
    {
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( )<<
                 computeObservationAndFullPrecisionLinkEndData(
                     time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates ) ).finished( );
    }

    //! Function to compute observation, as well as the state and time at each link end.
    /*!
     *  Function to compute observation, as well as the state and time at each link end. The
     *  times and states of the link ends are given in double precision. They are returned by
     *  reference. This function returns a scalar, as opposed to the
     *  general vector return type of the unspecialized ObservationModel
     *  computeObservationsAndLinkEndData function.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Calculated observable value.
     */
    ObservationScalarType computeObservationAndLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< double >& linkEndTimes,
                std::vector< basic_mathematics::Vector6d >& linkEndStates ) const
    {
        std::vector< TimeType > fullLinkEndTimes;
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > fullLinkEndStates;
        ObservationScalarType observation = computeObservationAndFullPrecisionLinkEndData(
                    time, linkEndAssociatedWithTime, fullLinkEndTimes, fullLinkEndStates );
        linkEndTimes = utilities::staticCastSVectorToTVector< TimeType, double >( fullLinkEndTimes );
        linkEndStates = utilities::staticCastSEigenTypeVectorToTEigenTypeVector<
                Eigen::Matrix< StateScalarType, 6, 1 >, basic_mathematics::Vector6d, double >( fullLinkEndStates );
        return observation;
    }

    //! Function to compute observation, as well as the state and time at each link end.
    /*!
     *  Function to compute observation, as well as the state and time at each link end. The
     *  times and states of the link ends are given in double precision. They are returned by
     *  reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Calculated observable value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeObservationsAndLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< basic_mathematics::Vector6d >& linkEndStates )
    {
        Eigen::Matrix< ObservationScalarType, 1, 1 > observations;
        observations.x( ) = computeObservationAndLinkEndData( time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
        observations += this->observationBiasCalculator_->getObservationBias( linkEndTimes ).template cast< ObservationScalarType >( );

        return observations;
    }

    //! Function to return the size of the observable (=1).
    /*!
     *  Function to return the size of the observable (=1).
     *  \return Size of the observable (=1).
     */
    int getObservationSize( )
    {
        return 1;
    }
protected:


    //! Type of observable, used for derived class type identification without explicit casts.
    ObservableType observableType_;

    //! Object for calculating system-dependent errors in the observable.
    /*!
     *  Object for calculating system-dependent errors in the observable, i.e. deviations from the
     *  physically true observable
     */
    boost::shared_ptr< ObservationBiasInterface > observationBiasCalculator_;

};

}

}
#endif // OBSERVATIONMODEL_H
