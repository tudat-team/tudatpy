/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_models
{

//! Enum listing types of observation biases that are availabe
enum ObservationBiasTypes
{
    multiple_observation_biases,
    constant_additive_bias,
    constant_multiplicative_bias
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
    ObservationBias( ){ }

    //! Destructor
    virtual ~ObservationBias( ){ }

    //! Pure virtual function to retrieve the observation bias.
    /*!
     * Pure virtual function to retrieve the observation bias as a function of the observation times and states (which are
     * typically computed by the ObservationModel)
     * \param linkEndTimes List of times at each link end during observation.
     * \param linkEndStates List of states at each link end during observation.
     * \return Observation bias at given times and states.
     */
    virtual Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue =
            ( Eigen::Matrix< double, ObservationSize, 1 >( ) << TUDAT_NAN ).finished( ) ) = 0;

    //! Function to return the size of the associated observation
    /*!
     * Function to return the size of the associated observation
     * \return Size of the associated observation
     */
    int getObservationSize( )
    {
        return ObservationSize;
    }
};

//! Class for a constant observation bias of a given size
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


    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return observationBias_;
    }

    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& observationBias )
    {
        observationBias_ = observationBias;
    }

    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return observationBias_;
    }

    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& observationBias )
    {
        if( observationBias.rows( ) == ObservationSize )
        {
            observationBias_ = observationBias;
        }
        else
        {
            throw std::runtime_error( "Error whem resetting constant bias, size is inconsistent" );
        }
    }


private:

    //! Constant (entry-wise) observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > observationBias_;

};

//! Class for a constant relative observation bias of a given size
template< int ObservationSize = 1 >
class ConstantRelativeObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observationBias Constant (entry-wise) observation bias.
     */
    ConstantRelativeObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > relativeObservationBias ):
        relativeObservationBias_( relativeObservationBias ){ }

    //! Destructor
    ~ConstantRelativeObservationBias( ){ }

    //! Function to retrieve the constant observation bias.
    /*!
     * Function to retrieve the constant observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \return Relative observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservableValue )
    {
        return relativeObservationBias_.cwiseProduct( currentObservableValue );
    }

    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return relativeObservationBias_;
    }

    void resetConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 >& relativeObservationBias )
    {
        relativeObservationBias_ = relativeObservationBias;
    }


    Eigen::VectorXd getTemplateFreeConstantObservationBias( )
    {
        return relativeObservationBias_;
    }

    void resetConstantObservationBiasTemplateFree( const Eigen::VectorXd& relativeObservationBias )
    {
        if( relativeObservationBias.rows( ) == ObservationSize )
        {
            relativeObservationBias_ = relativeObservationBias;
        }
        else
        {
            throw std::runtime_error( "Error whem resetting constant bias, size is inconsistent" );
        }
    }

private:

    //! Constant (entry-wise) observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > relativeObservationBias_;

};

template< int ObservationSize = 1 >
class MultiTypeObservationBias: public ObservationBias< ObservationSize >
{
public:


    MultiTypeObservationBias( const std::vector< boost::shared_ptr< ObservationBias< ObservationSize > > > biasList ):
        biasList_( biasList ){ }

    //! Destructor
    ~MultiTypeObservationBias( ){ }

    //! Function to retrieve the constant observation bias.
    /*!
     * Function to retrieve the constant observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \return Constant observation bias.
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

    std::vector< boost::shared_ptr< ObservationBias< ObservationSize > > > getBiasList( )
    {
        return biasList_;
    }


private:

    std::vector< boost::shared_ptr< ObservationBias< ObservationSize > > > biasList_;

};

template< int ObservationSize >
ObservationBiasTypes getObservationBiasType(
        const boost::shared_ptr< ObservationBias< ObservationSize > > biasObject )
{
    ObservationBiasTypes biasType;
    if( boost::dynamic_pointer_cast< ConstantObservationBias< ObservationSize > >( biasObject ) != NULL )
    {
        biasType = constant_additive_bias;
    }
    else if( boost::dynamic_pointer_cast< ConstantRelativeObservationBias< ObservationSize > >( biasObject ) != NULL )
    {
        biasType = constant_multiplicative_bias;
    }
    else if( boost::dynamic_pointer_cast< MultiTypeObservationBias< ObservationSize > >( biasObject ) != NULL )
    {
        biasType = multiple_observation_biases;
    }
    else
    {
        std::string errorMessage = "Error, did not recognize observation bias when retrieveing bias type";
        throw std::runtime_error( errorMessage );
    }
    return biasType;
 }

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONMODEL_H
