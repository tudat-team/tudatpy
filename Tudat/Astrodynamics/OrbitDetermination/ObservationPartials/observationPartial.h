#ifndef OBSERVATIONPARTIAL_H
#define OBSERVATIONPARTIAL_H

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

class PositionPartialScaling
{
public:

    //! Destructor
    virtual ~PositionPartialScaling( ){ }

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

    ObservationPartial( const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        parameterIdentifier_( parameterIdentifier ){ }

    //! Virtual destructor.
    /*!
     *  Virtual destrcutor.
     */
    virtual ~ObservationPartial( ) { }

    //! Pure virtual function to calculate the obsevration partial(s) at required time(s)
    /*!
     *  Pure virtual function to calculate the obsevration partial(s) at required time(s). Derived class functions implement this
     *  for a specific observable and parameter.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end times.
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver ) = 0;

    estimatable_parameters::EstimatebleParameterIdentifier getParameterIdentifier( )
    {
        return parameterIdentifier_;
    }


protected:
    estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier_;


};


//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 1. The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 2. The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 2 > > > SingleLinkObservationTwoPartialList;

typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 3 > > > SingleLinkObservationThreePartialList;



}

}

#endif // OBSERVATIONPARTIAL_H
