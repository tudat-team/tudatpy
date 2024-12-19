/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/math/basic/mathematicalConstants.h"

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
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    virtual void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                         const std::vector< double >& times,
                         const observation_models::LinkEndType fixedLinkEnd,
                         const Eigen::VectorXd currentObservation ) = 0;


};

template< int ObservationSize >
class DirectPositionPartialScaling: public PositionPartialScaling
{
public:

    DirectPositionPartialScaling( const observation_models::ObservableType observableType ):
        observableType_( observableType )
    { }

    virtual ~DirectPositionPartialScaling( ){ }

    virtual observation_models::LinkEndType getCurrentLinkEndType( ) = 0;

    virtual Eigen::Matrix< double, ObservationSize, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType ) = 0;

    virtual Eigen::Matrix< double, ObservationSize, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return Eigen::MatrixXd::Zero( ObservationSize, 3 );
    }

    virtual Eigen::Matrix< double, ObservationSize, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        throw std::runtime_error( "Error, fixed-time position partial not yet implemented for observable " + std::to_string( observableType_ ) );
    }

    virtual Eigen::Matrix< double, ObservationSize, 3 > getFixedTimeVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return Eigen::MatrixXd::Zero( ObservationSize, 3 );
    }


    virtual Eigen::Matrix< double, ObservationSize, 1 > getLightTimePartialScalingFactor( ) = 0;

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    virtual bool isVelocityScalingNonZero( )
    {
        return false;
    }

    virtual std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > getLinkIndependentPartials(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType )
    {
        if( useLinkIndependentPartials( ) == false )
        {
            throw std::runtime_error( "Error when getting link independent partials of observable type " +
                                      observation_models::getObservableName( observableType_) +
                                      ", no such model is defined." );
        }
        else
        {
            throw std::runtime_error( "Error when getting link independent partials of observable type " +
                                      observation_models::getObservableName( observableType_) +
                                      ", derived class model is not implemented." );
        }
    }

    virtual Eigen::Matrix< double, ObservationSize, 3 > getLightTimeGradientPartialScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        throw std::runtime_error( "Error when getting light time gradient partials of observable type " +
                                  observation_models::getObservableName( observableType_) +
                                  ", derived class model is not implemented." );
    }


    virtual bool useLinkIndependentPartials( )
    {
        return false;
    }

    virtual bool useLightTimeGradientPartials( )
    {
        return false;
    }

protected:

    observation_models::ObservableType observableType_;

    bool doesVelocityScalingFactorExist_;
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

    //! Pure virtual function to calculate the observation partial(s) at required time(s) and state(s)
    /*!
     *  Pure virtual function to calculate the observation partial(s) at required time(s) and state(s). States and times
     *  are typically obtained from evaluation of associated observation model.
     *  Derived class functions implement this for a specific observable.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes).
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) ) = 0;

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


//! Class for computing the derivative of any observable w.r.t. a constant absolute observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a constant absolute observation bias. Note that this partial is
 *  distinct from most other ObservationPartial partial derived classes, as its implementation is based on the parameter
 *  (constant observation bias), not the type of observable: the implementation is identical for each observable.
 */
template< int ObservationSize >
class ObservationPartialWrtConstantAbsoluteBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     */
    ObservationPartialWrtConstantAbsoluteBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::constant_additive_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds )
    {
        // Compute partial (vector of ObservationSize with 1.0 entries).
        constantPartial_ = Eigen::Matrix< double, ObservationSize, ObservationSize >::Identity( );
    }

    //! Destructor
    ~ObservationPartialWrtConstantAbsoluteBias( ){ }

    //! Function to calculate the observation partial w.r.t. constant absolute bias
    /*!
     *  Function to calculate the observation partial w.r.t. constant absolute bias. Note that output is independent of input.
     *  Associated time defined at times[ 0 ].
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        return { std::make_pair( constantPartial_, times.at( 0 ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Observation partial: constant for all conditions.
    Eigen::Matrix< double, ObservationSize, ObservationSize > constantPartial_;

};


//! Class for computing the derivative of any observable w.r.t. an arc-wise constant absolute observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a n arc-wiseconstant absolute observation bias. Note that this
 *  partial is distinct from most other ObservationPartial partial derived classes, as its implementation is based on the
 *  parameter (arc-wise constant observation bias), not the type of observable: the implementation is identical for each
 *  observable.
 */
template< int ObservationSize >
class ObservationPartialWrtArcWiseAbsoluteBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param arcLookupScheme Object used to determine the index from observationBiases_ to be used, based on the current time.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param numberOfArcs Number of arcs for which biases are defined
     */
    ObservationPartialWrtArcWiseAbsoluteBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds,
                                               const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                               const int linkEndIndex,
                                               const int numberOfArcs ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::arcwise_constant_additive_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds ), arcLookupScheme_( arcLookupScheme ),
        linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs )
    {
        constantPartial_ = Eigen::Matrix< double, ObservationSize, ObservationSize >::Identity( );
        totalPartial_ = Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >::Zero(
                    ObservationSize, ObservationSize * numberOfArcs_ );
    }

    //! Destructor
    ~ObservationPartialWrtArcWiseAbsoluteBias( ){ }

    //! Function to calculate the observation partial w.r.t. arc-wise constant absolute bias
    /*!
     *  Function to calculate the observation partial w.r.t. arc-wise constant absolute bias. Note that output is independent of
     *  input. Associated time defined by linkEndIndex_.
     *  \param states Link end states (unused).
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
//        std::cout << "TEST CALCULATE PARTIALS W.R.T. ARC WISE CONSTANT BIAS" << "\n\n";
        int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );


        totalPartial_.setZero( );
        if( arcLookupScheme_->getMinimumValue( ) <= times.at( linkEndIndex_ ) )
        {
            int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );
            totalPartial_.block( 0, currentIndex * ObservationSize, ObservationSize, ObservationSize ) = constantPartial_;
        }

        return { std::make_pair( totalPartial_, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

       //! Observation partial for singe arc: constant for all conditions.
    Eigen::Matrix< double, ObservationSize, ObservationSize > constantPartial_;

    //! Pre-allocated partial vector
    Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > totalPartial_;

};

//! Class for computing the derivative of any observable w.r.t. a constant relative observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a constant relative observation bias. Note that this partial is
 *  distinct from most other ObservationPartial partial derived classes, as its implementation is based on the parameter
 *  (constant observation bias), not the type of observable: the implementation is identical for each observable.
 */
template< int ObservationSize >
class ObservationPartialWrtConstantRelativeBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     */
    ObservationPartialWrtConstantRelativeBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::constant_additive_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds )
    {  }

    //! Destructor
    ~ObservationPartialWrtConstantRelativeBias( ){ }

    //! Function to calculate the observation partial w.r.t. constant relative bias
    /*!
     *  Function to calculate the observation partial w.r.t. constant bias. Note that output is independent of input times and
     *  states. Associated time defined at times[ 0 ].
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed.
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        return { std::make_pair( currentObservation, times.at( 0 ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;
};

//! Class for computing the derivative of any observable w.r.t. an arc-wise constant relative observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a n arc-wiseconstant relative observation bias. Note that this
 *  partial is distinct from most other ObservationPartial partial derived classes, as its implementation is based on the
 *  parameter (arc-wise constant observation bias), not the type of observable: the implementation is identical for each
 *  observable.
 */
template< int ObservationSize >
class ObservationPartialWrtArcWiseRelativeBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param arcLookupScheme Object used to determine the index from observationBiases_ to be used, based on the current time.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param numberOfArcs Number of arcs for which biases are defined
     */
    ObservationPartialWrtArcWiseRelativeBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds,
                                               const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                               const int linkEndIndex,
                                               const int numberOfArcs ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::arcwise_constant_relative_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds ), arcLookupScheme_( arcLookupScheme ),
        linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs )
    {
       totalPartial_ = Eigen::VectorXd::Zero( ObservationSize * numberOfArcs_ );
    }

    //! Destructor
    ~ObservationPartialWrtArcWiseRelativeBias( ){ }

    //! Function to calculate the observation partial w.r.t. arc-wise constant relative bias
    /*!
     *  Function to calculate the observation partial w.r.t. arc-wise constant relative bias. Note that output is independent of
     *  input. Associated time defined by linkEndIndex_.
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        totalPartial_.setZero( );

        if( arcLookupScheme_->getMinimumValue( ) <= times.at( linkEndIndex_ ) )
        {
            int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );
            totalPartial_.segment( currentIndex * ObservationSize, ObservationSize ) = currentObservation;
        }
        return { std::make_pair( totalPartial_, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

    //! Pre-allocated partial vector
    Eigen::VectorXd totalPartial_;

};



//extern template class ObservationPartial< 1 >;
//extern template class ObservationPartial< 2 >;
//extern template class ObservationPartial< 3 >;

//extern template class ObservationPartialWrtConstantAbsoluteBias< 1 >;
//extern template class ObservationPartialWrtConstantAbsoluteBias< 2 >;
//extern template class ObservationPartialWrtConstantAbsoluteBias< 3 >;

//extern template class ObservationPartialWrtArcWiseAbsoluteBias< 1 >;
//extern template class ObservationPartialWrtArcWiseAbsoluteBias< 2 >;
//extern template class ObservationPartialWrtArcWiseAbsoluteBias< 3 >;

//extern template class ObservationPartialWrtConstantRelativeBias< 1 >;
//extern template class ObservationPartialWrtConstantRelativeBias< 2 >;
//extern template class ObservationPartialWrtConstantRelativeBias< 3 >;

//extern template class ObservationPartialWrtArcWiseRelativeBias< 1 >;
//extern template class ObservationPartialWrtArcWiseRelativeBias< 2 >;
//extern template class ObservationPartialWrtArcWiseRelativeBias< 3 >;




//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 1.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 2.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 2 > > > SingleLinkObservationTwoPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 3.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > > SingleLinkObservationThreePartialList;


}

}

#endif // OBSERVATIONPARTIAL_H
