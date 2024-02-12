/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONBIAS_H
#define TUDAT_OBSERVATIONBIAS_H

#include <vector>
#include <iostream>

#include <memory>

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/math/interpolators/lookupScheme.h"


namespace tudat
{

namespace observation_models
{

//! Enum listing types of observation biases that are availabe
enum ObservationBiasTypes
{
    multiple_observation_biases,
    constant_absolute_bias,
    constant_relative_bias,
    arc_wise_constant_absolute_bias,
    arc_wise_constant_relative_bias,
    constant_time_drift_bias,
    arc_wise_time_drift_bias,
    constant_time_bias,
    arc_wise_time_bias
};

//! Base class (non-functional) for describing observation biases
/*!
 * Base class (non-functional) for describing observation biases. In this context, an observation bias denotes any deviation
 * from the ideal observable between two reference points.
 */
template< int ObservationSize = 1 >
class ObservationBias
{
public:

    //! Constructor
    ObservationBias( const bool hasTimeBias = false ): hasTimeBias_( hasTimeBias ){ }

    //! Destructor
    virtual ~ObservationBias( ){ }

    //! Pure virtual function to retrieve the observation bias.
    /*!
     * Pure virtual function to retrieve the observation bias as a function of the observation times and states (which are
     * typically computed by the ObservationModel)
     * \param linkEndTimes List of times at each link end during observation.
     * \param linkEndStates List of states at each link end during observation.
     * \param currentObservableValue  Unbiased value of the observable (default NAN for compatibility purposes with original
     * version of code).
     * \return Observation bias at given times and states.
     */
    virtual Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) ) = 0;

    virtual double getTimeBias(
        const double nominalObservationTime,
        const LinkEndType referenceLinkEnd )
    {
        if( hasTimeBias_ )
        {
            throw std::runtime_error( "Error when getting time bias, value should be non-zero, but none is implemented." );
        }
        return 0.0;
    }

    //! Function to return the size of the associated observation
    /*!
     * Function to return the size of the associated observation
     * \return Size of the associated observation
     */
    int getObservationSize( )
    {
        return ObservationSize;
    }

    bool getHasTimeBias( )
    {
        return hasTimeBias_;
    }

protected:

    bool hasTimeBias_;
};

//! Class for a constant absolute observation bias of a given size
/*!
 *  Class for a constant absolute observation bias of a given size. For unbiases observation h and bias A, the biased observation
 *  is computed as h + A
 */
template< int ObservationSize = 1 >
class ConstantObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observationBias Constant (entry-wise) observation bias.
     */
    ConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > observationBias ):
        observationBias_( observationBias ){ }

    //! Destructor
    ~ConstantObservationBias( ){ }

    //! Function to retrieve the constant observation bias.
    /*!
     * Function to retrieve the constant observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return observationBias_;
    }


    //! Function retrieve the constant (entry-wise) absolute observation bias.
    /*!
     * Function retrieve the constant (entry-wise) absolute observation bias.
     * \return The constant (entry-wise) absolute observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return observationBias_;
    }

    //! Function to reset the constant (entry-wise) absolute observation bias.
    /*!
     * Function to reset the constant (entry-wise) absolute observation bias.
     * \param observationBias The new constant (entry-wise) absolute observation bias.
     */
    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& observationBias )
    {
        observationBias_ = observationBias;
    }

    //! Function retrieve the constant (entry-wise) absolute observation bias as a variable-size vector.
    /*!
     * Function retrieve the constant (entry-wise) absolute observation bias as a variable-size vector
     * \return The constant (entry-wise) absolute observation bias.
     */
    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return observationBias_;
    }

    //! Function to reset the constant (entry-wise) absolute observation bias with variable-size input.
    /*!
     * Function to reset the constant (entry-wise) absolute observation bias with variable-size input. Input VectorXd size
     * must match ObservationSize class template parameter.
     * \param observationBias The new constant (entry-wise) absolute observation bias.
     */
    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& observationBias )
    {
        if( observationBias.rows( ) == ObservationSize )
        {
            observationBias_ = observationBias;
        }
        else
        {
            throw std::runtime_error( "Error when resetting constant bias, size is inconsistent" );
        }
    }


private:

    //! Constant (entry-wise) observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > observationBias_;

};

//! Class for an arc-wise constant absolute observation bias of a given size
/*!
 *  Class for an  arc-wise constant absolute observation bias of a given size. For unbiases observation h and bias A,
 *  the biased observation is computed as h + A. The bias A is provided per arc, with the arc start times provided to the
 *  class constructor.
 */
template< int ObservationSize = 1 >
class ConstantArcWiseObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases Absolute biases, constant per arc
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used
     * in getObservationBias function.
     */
    ConstantArcWiseObservationBias(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& observationBiases,
            const int linkEndIndexForTime ):
        arcStartTimes_( arcStartTimes ), observationBiases_( observationBiases ), linkEndIndexForTime_( linkEndIndexForTime )
    {
        resetArcStartTimes( arcStartTimes_, observationBiases_ );
    }

    //! Destructor
    ~ConstantArcWiseObservationBias( ){ }

    //! Function to retrieve the observation bias, determining the current arc from linkEndTimes.
    /*!
     * Function to retrieve the constant observation bias, determining the current arc from linkEndTimes.
     * \param linkEndTimes List of times at each link end during observation
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return observationBiases_.at(
                    lookupScheme_->findNearestLowerNeighbour( linkEndTimes.at( linkEndIndexForTime_ ) ) );
    }

    //! Function retrieve the list of arc-wise constant absolute observation biases as a variable-size vector.
    /*!
     * Function retrieve the  list of absolute observation biases as a variable-size vector, conatenated first by
     * entry, and then by arc.
     * \return List of absolute observation biases as a variable-size vector, conatenated first by
     * entry, and then by arc.
     */
    std::vector< Eigen::VectorXd > getTemplateFreeConstantObservationBias( )
    {
        std::vector< Eigen::VectorXd > templateFreeObservationBiases;
        for( unsigned int i = 0; i < observationBiases_.size( ); i++ )
        {
            templateFreeObservationBiases.push_back( observationBiases_.at( i ) );
        }
        return templateFreeObservationBiases;
    }

    //! Function reset the list of absolute observation biases
    /*!
     *  Function reset the list of absolute observation biases
     *  \param observationBiases The new list of absolute observation biases as a variable-size vector, with the bias for arc i
     *  in index i of the input vector
     */
    void resetConstantObservationBiasTemplateFree( const std::vector< Eigen::VectorXd >& observationBiases )
    {
        if( observationBiases_.size( ) == observationBiases.size( ) )
        {
            for( unsigned int i = 0; i < observationBiases.size( ); i++ )
            {
                if( ! ( observationBiases.at( i ).rows( ) == ObservationSize ) )
                {
                    throw std::runtime_error( "Error when resetting arc-wise constant bias, single entry size is inconsistent" );
                }
                else
                {
                    observationBiases_[ i ] = observationBiases.at( i );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error when resetting arc-wise constant bias, size is inconsistent" );
        }
    }

    //! Function to reset the arc start times for the bias (and associated biases)
    /*!
     * Function to reset the arc start times for the bias (and associated biases)
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases Absolute biases, constant per arc
     */
    void resetArcStartTimes( const std::vector< double >& arcStartTimes,
                             const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& observationBiases )
    {
        arcStartTimes_ = arcStartTimes;
        observationBiases_ = observationBiases;

        if( arcStartTimes_.size( ) != observationBiases_.size( ) )
        {
            throw std::runtime_error( "Error when creating constant arc-wise biases, input is inconsistent" );
        }

        // Create current arc lookup scheme
        std::vector< double > lookupSchemeTimes = arcStartTimes_;
        lookupSchemeTimes.push_back( std::numeric_limits< double >::max( ) );
        lookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    lookupSchemeTimes );

    }
    //! Function to retrieve start times for arcs in which biases (observationBiases) are used
    /*!
     * Function to retrieve start times for arcs in which biases (observationBiases) are used
     * \return Start times for arcs in which biases (observationBiases) are used
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Function to retrieve link end index from which the 'current time' is determined
    /*!
     * Function to retrieve link end index from which the 'current time' is determined
     * \return Link end index from which the 'current time' is determined
     */
    int getLinkEndIndexForTime( )
    {
        return linkEndIndexForTime_;
    }

    //! Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
    /*!
     * Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
     * \return Object used to determine the index from observationBiases_ to be used, based on the current time.
     */
    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

private:

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! Absolute biases, constant per arc
    std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Class for a constant relative observation bias of a given size
/*!
 *  Class for a constant relative observation bias of a given size. For unbiases observation h and bias A, the biased observation
 *  is computed as h .* A, where .* is the component-wise multiplication.
 */
template< int ObservationSize = 1 >
class ConstantRelativeObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param relativeObservationBias Constant (entry-wise) observation bias.
     */
    ConstantRelativeObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > relativeObservationBias ):
        relativeObservationBias_( relativeObservationBias ){ }

    //! Destructor
    ~ConstantRelativeObservationBias( ){ }

    //! Function to retrieve the constant relative observation bias.
    /*!
     * Function to retrieve the constant relative observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue Unbiased value of the observable
     * \return Relative observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue )
    {
        return relativeObservationBias_.cwiseProduct( currentObservableValue );
    }

    //! Function retrieve the constant (entry-wise) relative observation bias.
    /*!
     * Function retrieve the constant (entry-wise) relative observation bias.
     * \return The constant (entry-wise) relative observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return relativeObservationBias_;
    }

    //! Function to reset the constant (entry-wise) relative observation bias.
    /*!
     * Function to reset the constant (entry-wise) relative observation bias.
     * \param relativeObservationBias The new constant (entry-wise) relative observation bias.
     */
    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& relativeObservationBias )
    {
        relativeObservationBias_ = relativeObservationBias;
    }

    //! Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector.
    /*!
     * Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector
     * \return The constant (entry-wise) relative observation bias.
     */
    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return relativeObservationBias_;
    }

    //! Function to reset the constant (entry-wise) relative observation bias with variable-size input.
    /*!
     * Function to reset the constant (entry-wise) relative observation bias with variable-size input. Input VectorXd size
     * must match ObservationSize class template parameter.
     * \param relativeObservationBias The new constant (entry-wise) relative observation bias.
     */
    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& relativeObservationBias )
    {
        if( relativeObservationBias.rows( ) == ObservationSize )
        {
            relativeObservationBias_ = relativeObservationBias;
        }
        else
        {
            throw std::runtime_error( "Error when resetting constant bias, size is inconsistent" );
        }
    }

private:

    //! Constant (entry-wise) relative observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > relativeObservationBias_;

};


//! Class for an arc-wise constant relative observation bias of a given size
/*!
 *  Class for an  arc-wise constant relative observation bias of a given size. For unbiases observation h and bias A,
 *  the biased observation is computed as h .* A, where .* is the component-wise multiplication.The bias A is provided per arc,
 *  with the arc start times provided to the class constructor.
 */
template< int ObservationSize = 1 >
class ConstantRelativeArcWiseObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases Relative biases, constant per arc
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used
     * in getObservationBias function.
     */
    ConstantRelativeArcWiseObservationBias(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& observationBiases,
            const int linkEndIndexForTime ):
        arcStartTimes_( arcStartTimes ), observationBiases_( observationBiases ), linkEndIndexForTime_( linkEndIndexForTime )
    {
        if( arcStartTimes_.size( ) != observationBiases_.size( ) )
        {
            throw std::runtime_error( "Error when creating constant arc-wise biases, input is inconsistent" );
        }

        // Create current arc lookup scheme
        std::vector< double > lookupSchemeTimes = arcStartTimes_;
        lookupSchemeTimes.push_back( std::numeric_limits< double >::max( ) );
        lookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    lookupSchemeTimes );
    }

    //! Destructor
    ~ConstantRelativeArcWiseObservationBias( ){ }

    //! Function to retrieve the constant observation bias, determining the current arc from linkEndTimes.
    /*!
     * Function to retrieve the constant observation bias, determining the current arc from linkEndTimes.
     * \param linkEndTimes List of times at each link end during observation
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return observationBiases_.at( lookupScheme_->findNearestLowerNeighbour( linkEndTimes.at( linkEndIndexForTime_ ) ) ).
                cwiseProduct( currentObservableValue );
    }

    //! Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector.
    /*!
         * Function retrieve the constant (entry-wise) relative observation bias as a variable-size vector
         * \return The constant (entry-wise) relative observation bias.
         */
    std::vector< Eigen::VectorXd > getTemplateFreeConstantObservationBias( )
    {
        std::vector< Eigen::VectorXd > templateFreeObservationBiases;
        for( unsigned int i = 0; i < observationBiases_.size( ); i++ )
        {
            templateFreeObservationBiases.push_back( observationBiases_.at( i ) );
        }
        return templateFreeObservationBiases;
    }

    //! Function to reset the constant (entry-wise) relative observation bias with variable-size input.
    /*!
     *  Function to reset the constant (entry-wise) relative observation bias with variable-size input. Input VectorXd size
     *  must match ObservationSize class template parameter.
     *  \param observationBiases The new constant arc-wise list of (entry-wise) relative observation bias, with the bias for arc i
     *  in index i of the input vector
     */
    void resetConstantObservationBiasTemplateFree( const std::vector< Eigen::VectorXd >& observationBiases )
    {
        if( observationBiases_.size( ) == observationBiases.size( ) )
        {
            for( unsigned int i = 0; i < observationBiases.size( ); i++ )
            {
                if( ! ( observationBiases.at( i ).rows( ) == ObservationSize ) )
                {
                    throw std::runtime_error( "Error when resetting arc-wise constant bias, single entry size is inconsistent" );
                }
                else
                {
                    observationBiases_[ i ] = observationBiases.at( i );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error when resetting arc-wise constant bias, size is inconsistent" );
        }
    }

    //! Function to retrieve start times for arcs in which biases (observationBiases) are used
    /*!
     * Function to retrieve start times for arcs in which biases (observationBiases) are used
     * \return Start times for arcs in which biases (observationBiases) are used
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Function to retrieve link end index from which the 'current time' is determined
    /*!
     * Function to retrieve link end index from which the 'current time' is determined
     * \return Link end index from which the 'current time' is determined
     */
    int getLinkEndIndexForTime( )
    {
        return linkEndIndexForTime_;
    }

    //! Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
    /*!
     * Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
     * \return Object used to determine the index from observationBiases_ to be used, based on the current time.
     */
    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

private:

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! Absolute biases, constant per arc
    std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Class for combining multiple observation bias models into a single bias value
/*!
 *  Class for combining multiple observation bias models into a single bias value. This class computes a list of biases,
 *  all based on the nominal, unbiased, observation and sums them up to form the total observation bias.
 */
template< int ObservationSize = 1 >
class MultiTypeObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param biasList List of bias objects that are to be combined.
     */
    MultiTypeObservationBias( const std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > biasList ):
        biasList_( biasList )
    {
        for( unsigned int i = 0; i < biasList_.size( ); i++ )
        {
            if( biasList_.at( i )->getHasTimeBias( ) )
            {
                this->hasTimeBias_ = true;
            }
        }
    }

    //! Destructor
    ~MultiTypeObservationBias( ){ }

    //! Function to retrieve the total observation bias.
    /*!
     * Function to retrieve the total observation bias.
     * \param linkEndTimes List of times at each link end during observation.
     * \param linkEndStates List of states at each link end during observation.
     * \param currentObservableValue  Unbiased value of the observable.
     * \return Total observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue )
    {
        Eigen::Matrix< double, ObservationSize, 1 > totalBias = Eigen::Matrix< double, ObservationSize, 1 >::Zero( );
        for( unsigned int i = 0; i < biasList_.size( ); i++ )
        {
            totalBias += biasList_.at( i )->getObservationBias( linkEndTimes, linkEndStates, currentObservableValue );
        }
        return totalBias;
    }

    double getTimeBias(
        const double nominalObservationTime,
        const LinkEndType referenceLinkEnd )
    {
        double timeBias = 0.0;
        if( this->hasTimeBias_ )
        {
            for ( unsigned int i = 0; i < biasList_.size( ); i++ )
            {
                if ( biasList_.at( i )->getHasTimeBias( ) )
                {
                    timeBias += biasList_.at( i )->getTimeBias( nominalObservationTime, referenceLinkEnd );
                }
            }
        }
        return timeBias;
    }

    //! Function to retrieve the list of bias objects that are to be combined.
    /*!
     * Function to retrieve the list of bias objects that are to be combined.
     * \return The list of bias objects that are to be combined.
     */
    std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > getBiasList( )
    {
        return biasList_;
    }


private:


    //! List of bias objects that are to be combined.
    std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > biasList_;

};

//! Class for a constant time drift observation bias of a given size
/*!
 *  Class for a constant time drift observation bias of a given size. For unbiases observation h and time drift c, the biased observation
 *  is computed as h(t) + c*t
 */
template< int ObservationSize = 1 >
class ConstantTimeDriftBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param timeDriftBias Constant (entry-wise) time drift bias.
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined
     * \param referenceEpoch Reference epoch at which the time drift is initialised
     */
    ConstantTimeDriftBias( const Eigen::Matrix< double, ObservationSize, 1 > timeDriftBias,
                                 const int linkEndIndexForTime,
                                 const double referenceEpoch ):
            timeDriftBias_( timeDriftBias ), linkEndIndexForTime_( linkEndIndexForTime ), referenceEpoch_( referenceEpoch ){ }

    //! Destructor
    ~ConstantTimeDriftBias( ){ }

    //! Function to retrieve the constant time drift bias.
    /*!
     * Function to retrieve the constant time drift bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant time drift bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return timeDriftBias_ * ( linkEndTimes.at( linkEndIndexForTime_ ) - referenceEpoch_ );
    }


    //! Function retrieve the constant (entry-wise) time drift bias.
    /*!
     * Function retrieve the constant (entry-wise) time drift bias.
     * \return The constant (entry-wise) time drift bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return timeDriftBias_;
    }

    //! Function to reset the constant (entry-wise) time drift bias.
    /*!
     * Function to reset the constant (entry-wise) time drift bias.
     * \param observationTimeBias The new constant (entry-wise) time drift bias.
     */
    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& timeDriftBias )
    {
        timeDriftBias_ = timeDriftBias;
    }

    //! Function retrieve the constant (entry-wise) time drift bias as a variable-size vector.
    /*!
     * Function retrieve the constant (entry-wise) time drift bias as a variable-size vector
     * \return The constant (entry-wise) time drift bias.
     */
    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return timeDriftBias_;
    }

    //! Function to reset the constant (entry-wise) time drift bias with variable-size input.
    /*!
     * Function to reset the constant (entry-wise) time drift bias with variable-size input. Input VectorXd size
     * must match ObservationSize class template parameter.
     * \param observationTimeBias The new constant (entry-wise) time drift bias.
     */
    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& timeDriftBias )
    {
        if( timeDriftBias.rows( ) == ObservationSize )
        {
            timeDriftBias_ = timeDriftBias;
        }
        else
        {
            throw std::runtime_error( "Error when resetting time drift bias, size is inconsistent" );
        }
    }

private:

    //! Constant (entry-wise) time drift bias.
    Eigen::Matrix< double, ObservationSize, 1 > timeDriftBias_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Reference epoch at which the time drift is initialised.
    double referenceEpoch_;

};


//! Class for an arc-wise time drift bias of a given size
/*!
*  Class for an arc-wise time drift bias of a given size. For unbiases observation h and time drift c, the biased observation
*  is computed as h(t) + c*t, where .* is the component-wise multiplication.The bias c is provided per arc,
*  with the arc start times provided to the class constructor.
*/
template< int ObservationSize = 1 >
class ArcWiseTimeDriftBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param arcStartTimes Start times for arcs in which biases (timeDriftBiases) are used
     * \param timeDriftBiases Time drift biases, constant per arc
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used
     * in getObservationBias function.
     * \param referenceEpochs Reference epochs (per arc) at which the time drift is initialised
     */
    ArcWiseTimeDriftBias(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::Matrix< double, ObservationSize, 1 > >& timeDriftBiases,
            const int linkEndIndexForTime,
            const std::vector< double > referenceEpochs ):
            arcStartTimes_( arcStartTimes ), timeDriftBiases_( timeDriftBiases ), linkEndIndexForTime_( linkEndIndexForTime ), referenceEpochs_( referenceEpochs )
    {
        if( arcStartTimes_.size( ) != timeDriftBiases_.size( ) )
        {
            throw std::runtime_error( "Error when creating arc-wise time drift biases, input is inconsistent" );
        }

        // Create current arc lookup scheme
        std::vector< double > lookupSchemeTimes = arcStartTimes_;
        lookupSchemeTimes.push_back( std::numeric_limits< double >::max( ) );
        lookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                lookupSchemeTimes );
    }

    //! Destructor
    ~ArcWiseTimeDriftBias( ){ }

    //! Function to retrieve the arc-wise time drift bias, determining the current arc from linkEndTimes.
    /*!
     * Function to retrieve the arc-wise time drift bias, determining the current arc from linkEndTimes.
     * \param linkEndTimes List of times at each link end during observation
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return arc-wise time drift bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return timeDriftBiases_.at( lookupScheme_->findNearestLowerNeighbour( linkEndTimes.at( linkEndIndexForTime_ ) ) ) *
               ( linkEndTimes.at( linkEndIndexForTime_ )
                 - referenceEpochs_.at( lookupScheme_->findNearestLowerNeighbour( linkEndTimes.at( linkEndIndexForTime_ ) ) ) );
    }

    //! Function retrieve the constant (entry-wise) time drift bias as a variable-size vector.
    /*!
         * Function retrieve the constant (entry-wise) time drift bias as a variable-size vector
         * \return The constant (entry-wise) time drift bias.
         */
    std::vector< Eigen::VectorXd > getTemplateFreeConstantObservationBias( )
    {
        std::vector< Eigen::VectorXd > templateFreeObservationBiases;
        for( unsigned int i = 0; i < timeDriftBiases_.size( ); i++ )
        {
            templateFreeObservationBiases.push_back( timeDriftBiases_.at( i ) );
        }
        return templateFreeObservationBiases;
    }

    //! Function to reset the constant (entry-wise) time drift bias with variable-size input.
    /*!
     *  Function to reset the constant (entry-wise) time drift bias with variable-size input. Input VectorXd size
     *  must match ObservationSize class template parameter.
     *  \param observationTimeBiases The new constant arc-wise list of (entry-wise) time drift bias, with the bias for arc i
     *  in index i of the input vector
     */
    void resetConstantObservationBiasTemplateFree( const std::vector< Eigen::VectorXd >& timeDriftBiases )
    {
        if( timeDriftBiases_.size( ) == timeDriftBiases.size( ) )
        {
            for( unsigned int i = 0; i < timeDriftBiases.size( ); i++ )
            {
                if( ! ( timeDriftBiases.at( i ).rows( ) == ObservationSize ) )
                {
                    throw std::runtime_error( "Error when resetting arc-wise time drift bias, single entry size is inconsistent" );
                }
                else
                {
                    timeDriftBiases_[ i ] = timeDriftBiases.at( i );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error when resetting arc-wise time drift bias, size is inconsistent" );
        }
    }

    //! Function to retrieve start times for arcs in which biases (observationBiases) are used
    /*!
     * Function to retrieve start times for arcs in which biases (observationBiases) are used
     * \return Start times for arcs in which biases (observationBiases) are used
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Function to retrieve link end index from which the 'current time' is determined
    /*!
     * Function to retrieve link end index from which the 'current time' is determined
     * \return Link end index from which the 'current time' is determined
     */
    int getLinkEndIndexForTime( )
    {
        return linkEndIndexForTime_;
    }

    //! Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
    /*!
     * Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
     * \return Object used to determine the index from observationBiases_ to be used, based on the current time.
     */
    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

private:

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! Time drift biases, constant per arc
    std::vector< Eigen::Matrix< double, ObservationSize, 1 > > timeDriftBiases_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;

    //! Reference epochs at which the time drift is initialised, for each arc.
    std::vector< double > referenceEpochs_;
};


//! Class for a constant time observation bias of a given size
/*!
 *  Class for a constant time observation bias of a given size. For unbiases observation h and time bias c, the biased observation
 *  is computed as h(t+c)
 */
template< int ObservationSize = 1 >
class ConstantTimeBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param timeBias Constant (entry-wise) time bias.
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined
     */
    ConstantTimeBias( const double timeBias,
                      const int linkEndIndexForTime ): ObservationBias< ObservationSize >( true ),
            timeBias_( timeBias ), linkEndIndexForTime_( linkEndIndexForTime ){ }

    //! Destructor
    ~ConstantTimeBias( ){ }

    //! Function to retrieve the constant time drift bias.
    /*!
     * Function to retrieve the constant time drift bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return Constant time drift bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return Eigen::Matrix< double, ObservationSize, 1 >::Zero( );
    }

    double getTimeBias(
        const double nominalObservationTime,
        const LinkEndType referenceLinkEnd )
    {
        return timeBias_;
    }


    //! Function to reset the constant (entry-wise) time bias.
    /*!
     * Function to reset the constant (entry-wise) time bias.
     * \param observationTimeBias The new constant (entry-wise) time bias.
     */
    void resetConstantObservationBias( const double timeBias )
    {
        timeBias_ = timeBias;
    }

    //! Function retrieve the constant (entry-wise) time bias as a variable-size vector.
    /*!
     * Function retrieve the constant (entry-wise) time bias as a variable-size vector
     * \return The constant (entry-wise) time bias.
     */
    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        Eigen::VectorXd bias = Eigen::VectorXd( 1 );
        bias[ 0 ] = timeBias_;
        return bias;
    }

    //! Function to reset the constant (entry-wise) time bias with variable-size input.
    /*!
     * Function to reset the constant (entry-wise) time bias with variable-size input. Input VectorXd size
     * must match ObservationSize class template parameter.
     * \param observationTimeBias The new constant (entry-wise) time bias.
     */
    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& timeBias )
    {
        timeBias_ = timeBias[ 0 ];
    }

    int getLinkEndIndexForTime( ) const
    {
        return linkEndIndexForTime_;
    }

private:

    //! Constant (entry-wise) time bias.
    double timeBias_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

};


//! Class for an arc-wise time bias of a given size
/*!
*  Class for an arc-wise time bias of a given size. For unbiases observation h and time bias c, the biased observation
*  is computed as h(t+c).The bias c is provided per arc,
*  with the arc start times provided to the class constructor.
*/
template< int ObservationSize = 1 >
class ArcWiseTimeBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param arcStartTimes Start times for arcs in which biases (timeDriftBiases) are used
     * \param timeBiases Time biases, constant per arc
     * \param linkEndIndexForTime Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used
     * in getObservationBias function.
     */
    ArcWiseTimeBias(
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& timeBiases,
            const int linkEndIndexForTime ): ObservationBias< ObservationSize >( true ),
            arcStartTimes_( arcStartTimes ), timeBiases_( timeBiases ), linkEndIndexForTime_( linkEndIndexForTime )
    {
        if( arcStartTimes_.size( ) != timeBiases_.size( ) )
        {
            throw std::runtime_error( "Error when creating arc-wise time biases, input is inconsistent" );
        }

        // Create current arc lookup scheme
        std::vector< double > lookupSchemeTimes = arcStartTimes_;
        lookupSchemeTimes.push_back( std::numeric_limits< double >::max( ) );
        lookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                lookupSchemeTimes );
    }

    //! Destructor
    ~ArcWiseTimeBias( ){ }

    //! Function to retrieve the arc-wise time bias, determining the current arc from linkEndTimes.
    /*!
     * Function to retrieve the arc-wise time bias, determining the current arc from linkEndTimes.
     * \param linkEndTimes List of times at each link end during observation
     * \param linkEndStates List of states at each link end during observation (unused).
     * \param currentObservableValue  Unbiased value of the observable (unused and default NAN).
     * \return arc-wise time bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) )
    {
        return Eigen::Matrix< double, ObservationSize, 1 >::Zero( );
    }

    double getTimeBias(
        const double nominalObservationTime,
        const LinkEndType referenceLinkEnd )
    {
        return timeBiases_.at( lookupScheme_->findNearestLowerNeighbour( nominalObservationTime ) );
    }

    //! Function retrieve the constant (entry-wise) time bias as a variable-size vector.
    /*!
         * Function retrieve the constant (entry-wise) time bias as a variable-size vector
         * \return The constant (entry-wise) time bias.
         */
    std::vector< Eigen::VectorXd > getTemplateFreeConstantObservationBias( )
    {
        std::vector< Eigen::VectorXd > templateFreeObservationBiases;
        for( unsigned int i = 0; i < timeBiases_.size( ); i++ )
        {
            Eigen::VectorXd biases = ( Eigen::Vector1d( ) << timeBiases_.at( i ) ).finished( );
            templateFreeObservationBiases.push_back( biases );
        }
        return templateFreeObservationBiases;
    }

    //! Function to reset the constant (entry-wise) time bias with variable-size input.
    /*!
     *  Function to reset the constant (entry-wise) time bias with variable-size input. Input VectorXd size
     *  must match ObservationSize class template parameter.
     *  \param timeBiases The new constant arc-wise list of (entry-wise) time bias, with the bias for arc i
     *  in index i of the input vector
     */
    void resetConstantObservationBiasTemplateFree( const std::vector< Eigen::VectorXd >& timeBiases )
    {
        if( timeBiases_.size( ) == timeBiases.size( ) )
        {
            for( unsigned int i = 0; i < timeBiases.size( ); i++ )
            {
                if( ( timeBiases.at( i ).rows( ) != 1 ) )
                {
                    throw std::runtime_error( "Error when resetting arc-wise time bias, single entry size is inconsistent" );
                }
                else
                {
                    timeBiases_[ i ] = timeBiases.at( i )[ 0 ];
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error when resetting arc-wise time bias, size is inconsistent" );
        }
    }

    //! Function to retrieve start times for arcs in which biases (observationBiases) are used
    /*!
     * Function to retrieve start times for arcs in which biases (observationBiases) are used
     * \return Start times for arcs in which biases (observationBiases) are used
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Function to retrieve link end index from which the 'current time' is determined
    /*!
     * Function to retrieve link end index from which the 'current time' is determined
     * \return Link end index from which the 'current time' is determined
     */
    int getLinkEndIndexForTime( )
    {
        return linkEndIndexForTime_;
    }

    //! Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
    /*!
     * Function to retrieve object used to determine the index from observationBiases_ to be used, based on the current time.
     * \return Object used to determine the index from observationBiases_ to be used, based on the current time.
     */
    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

private:

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! Time biases, constant per arc
    std::vector< double > timeBiases_;

    //! Link end index from which the 'current time' is determined (e.g. entry from linkEndTimes used in getObservationBias
    //! function.
    int linkEndIndexForTime_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;

};


//! Function to retrieve the type of an observation bias
/*!
 *  Function to retrieve the type of an observation bias.
 *  \param biasObject Bias for which the type is to be determined
 *  \return Bias type for biasObject.
 */
template< int ObservationSize >
ObservationBiasTypes getObservationBiasType(
        const std::shared_ptr< ObservationBias< ObservationSize > > biasObject )
{
    ObservationBiasTypes biasType;

    // Check available bias types
    if( std::dynamic_pointer_cast< ConstantObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = constant_absolute_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantArcWiseObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = arc_wise_constant_absolute_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantRelativeObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = constant_relative_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantRelativeArcWiseObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = arc_wise_constant_relative_bias;
    }
    else if( std::dynamic_pointer_cast< MultiTypeObservationBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = multiple_observation_biases;
    }
    else if( std::dynamic_pointer_cast< ConstantTimeDriftBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = constant_time_drift_bias;
    }
    else if( std::dynamic_pointer_cast< ArcWiseTimeDriftBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = arc_wise_time_drift_bias;
    }
    else if( std::dynamic_pointer_cast< ConstantTimeBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = constant_time_bias;
    }
    else if( std::dynamic_pointer_cast< ArcWiseTimeBias< ObservationSize > >( biasObject ) != nullptr )
    {
        biasType = arc_wise_time_bias;
    }
    else
    {
        std::string errorMessage = "Error, did not recognize observation bias when retrieving bias type";
        throw std::runtime_error( errorMessage );
    }
    return biasType;
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONBIAS_H
