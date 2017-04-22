/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/observationBiasParameter.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

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
    virtual void update( const std::vector< Eigen::Vector6d >& linkEndStates,
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
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
            Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) ) = 0;

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

template< int ObservationSize >
class ObservationPartialWrtConstantAdditiveBias: public ObservationPartial< ObservationSize >
{
public:
    ObservationPartialWrtConstantAdditiveBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::constant_additive_observation_bias, linkEnds.begin( )->second ) ),
        observableType_( observableType ), linkEnds_( linkEnds )
    {
        constantPartial_ = Eigen::Matrix< double, ObservationSize, 1 >::Constant( 1.0 );
    }

    ~ObservationPartialWrtConstantAdditiveBias( ){ }

    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
            Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        return { std::make_pair( constantPartial_, times.at( 0 ) ) };
    }

private:
    observation_models::ObservableType observableType_;

    observation_models::LinkEnds linkEnds_;

    Eigen::Matrix< double, ObservationSize, 1 > constantPartial_;

};

template< int ObservationSize >
class ObservationPartialWrtConstantMultiplicativeBias: public ObservationPartial< ObservationSize >
{
public:
    ObservationPartialWrtConstantMultiplicativeBias( const observation_models::ObservableType observableType,
                                                     const observation_models::LinkEnds& linkEnds ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::constant_additive_observation_bias, linkEnds.begin( )->second ) ),
        observableType_( observableType ), linkEnds_( linkEnds )
    {
        constantPartial_ = Eigen::Matrix< double, ObservationSize, 1 >::Constant( 1.0 );
    }

    ~ObservationPartialWrtConstantMultiplicativeBias( ){ }

    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
            Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        return { std::make_pair( currentObservation, times.at( 0 ) ) };
    }

private:
    observation_models::ObservableType observableType_;

    observation_models::LinkEnds linkEnds_;

    Eigen::Matrix< double, ObservationSize, 1 > constantPartial_;

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



template< int ObservationSize >
boost::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtLinkProperty(
        const observation_models::LinkEnds& linkEnds,
        const observation_models::ObservableType observableType,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate )
{
    boost::shared_ptr< ObservationPartial< ObservationSize > > observationPartial;
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::constant_additive_observation_bias:
    {
        boost::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > constantBias =
                boost::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                    parameterToEstimate );
        if( constantBias == NULL )
        {
            throw std::runtime_error( "Error when making partial w.r.t. observation bias, type is inconsistent" );
        }
        else
        {
            if( linkEnds == constantBias->getLinkEnds( ) && observableType == constantBias->getObservableType( ) )
            {
                observationPartial = boost::make_shared< ObservationPartialWrtConstantAdditiveBias< ObservationSize > >(
                            observableType, linkEnds );
            }
        }
        break;
    }
    case estimatable_parameters::constant_relative_observation_bias:
    {
        boost::shared_ptr< estimatable_parameters::ConstantRelativeObservationBiasParameter > constantBias =
                boost::dynamic_pointer_cast< estimatable_parameters::ConstantRelativeObservationBiasParameter >(
                    parameterToEstimate );
        if( constantBias == NULL )
        {
            throw std::runtime_error( "Error when making partial w.r.t. observation bias, type is inconsistent" );
        }
        else
        {
            if( linkEnds == constantBias->getLinkEnds( ) && observableType == constantBias->getObservableType( ) )
            {
                observationPartial = boost::make_shared< ObservationPartialWrtConstantMultiplicativeBias< ObservationSize > >(
                            observableType, linkEnds );
            }
        }
        break;
    }
    default:
        break;
    }
    return observationPartial;
}

}

}

#endif // OBSERVATIONPARTIAL_H
