/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONPARTIAL_H
#define TUDAT_OBSERVATIONPARTIAL_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

//! Base class for scaling position partials to observable partials.
/*!
 *  Base class for scaling position partials to observable partials. For observables computed from the three-dimensional
 *  positions of reference points, the partial derivatives are computed from the partials of these positions. Each
 *  observation model requires a specific derived class of PositionPartialScaling, that computes how position partials
 *  are converted to ('scaled to') observable partials
 */
class PositionPartialScaling
{
public:

    //! Destructor
    virtual ~PositionPartialScaling( ){ }

    //! Update the scaling object to the current times and states (pure virtual).
    /*!
     *  Update the scaling object to the current times and states (pure virtual).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     */
    virtual void update( const std::vector< basic_mathematics::Vector6d >& linkEndStates,
                         const std::vector< double >& times,
                         const observation_models::LinkEndType fixedLinkEnd ) = 0;
};

//! Base class of partial derivative of an observable w.r.t. an estimated parameter.
/*!
 *  Base class of partial derivative of an observable w.r.t. an estimated parameter.
 *  Calculated partial function must be defined in each derived class for calculation of specific derivative.
 */
template< int ObservationSize >
class ObservationPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Parameter id of for specifc parameter of which partial is computed by object.
     */
    ObservationPartial( const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        parameterIdentifier_( parameterIdentifier ){ }

    //! Virtual destructor.
    /*!
     *  Virtual destrcutor.
     */
    virtual ~ObservationPartial( ) { }

    //! Pure virtual function to calculate the obsevration partial(s) at required time(s) and state(s)
    /*!
     *  Pure virtual function to calculate the obsevration partial(s) at required time(s) and state(s). States and times
     *  are typically obtained from evaluation of associated observation model.
     *  Derived class functions implement this for a specific observable.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver ) = 0;

    //! Function to get parameter id of for specifc parameter of which partial is computed by object.
    /*!
     * Function to get parameter id of for specifc parameter of which partial is computed by object.
     * \return Parameter id of for specifc parameter of which partial is computed by object.
     */
    estimatable_parameters::EstimatebleParameterIdentifier getParameterIdentifier( )
    {
        return parameterIdentifier_;
    }


protected:

    //! Parameter id of for specifc parameter of which partial is computed by object.
    estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier_;


};


//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 1.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 2.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 2 > > > SingleLinkObservationTwoPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 3.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 3 > > > SingleLinkObservationThreePartialList;



}

}

#endif // OBSERVATIONPARTIAL_H
