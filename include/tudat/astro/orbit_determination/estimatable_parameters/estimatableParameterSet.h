/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ESTIMATABLEPARAMETERSET_H
#define TUDAT_ESTIMATABLEPARAMETERSET_H

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <map>



#include <memory>
#include <Eigen/Geometry>

#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace estimatable_parameters
{

//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//extern template class EstimatableParameter< Eigen::Matrix< long double, Eigen::Dynamic, 1 > >;
//#endif

//! Container class for all parameters that are to be estimated.
/*!
 *  Container class for all parameters that are to be estimated. Class is templated with the scalar type used for the
 *  estimation of any initial dynamical states that may be included
 */
template< typename InitialStateParameterType = double >
class EstimatableParameterSet
{
public:

    //! Constructor of parameter set.
    /*!
     *  Constructor of parameter set.
     *  \param estimatedDoubleParameters List of double parameters that are estimated.
     *  \param estimatedVectorParameters List of vector parameters that are estimated.
     *  \param estimateInitialStateParameters List of initial dynamical states that are to be estimated.
     */
    EstimatableParameterSet(
            const std::vector< std::shared_ptr< EstimatableParameter< double > > >& estimatedDoubleParameters,
            const std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& estimatedVectorParameters,
            const std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix
            < InitialStateParameterType, Eigen::Dynamic, 1 > > > >& estimateInitialStateParameters =
            ( std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix
              < InitialStateParameterType, Eigen::Dynamic, 1 > > > >( ) ) ):
        estimatedDoubleParameters_( estimatedDoubleParameters ), estimatedVectorParameters_( estimatedVectorParameters ),
        totalConstraintSize_( 0 )
    {
        // Initialize total number of parameters to 0.
        estimatedParameterSetSize_ = 0;
        initialDynamicalStateParameterSize_ = 0;
        initialDynamicalSingleArcStateParameterSize_ = 0;
        initialDynamicalMultiArcStateParameterSize_ = 0;

        // Iterate over all double parameters and add to parameter size.
        for( unsigned int i = 0; i < estimateInitialStateParameters.size( ); i++ )
        {
            if( isDynamicalParameterSingleArc( estimateInitialStateParameters[ i ] ) )
            {
                estimateSingleArcInitialStateParameters_.push_back( estimateInitialStateParameters[ i ] );
            }
            else
            {
                estimateMultiArcInitialStateParameters_.push_back( estimateInitialStateParameters[ i ] );
            }
        }

        estimateInitialStateParameters_ = estimateSingleArcInitialStateParameters_;
        estimateInitialStateParameters_.insert(
                    estimateInitialStateParameters_.end( ), estimateMultiArcInitialStateParameters_.begin( ),
                    estimateMultiArcInitialStateParameters_.end( ) );

        for( unsigned int i = 0; i < estimateSingleArcInitialStateParameters_.size( ); i++ )
        {
            initialStateParameters_[ estimatedParameterSetSize_ ] = estimateSingleArcInitialStateParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimateInitialStateParameters_[ i ]->getParameterSize( ) ) );
            totalConstraintSize_ += estimateInitialStateParameters_[ i ]->getConstraintSize( );

            initialDynamicalSingleArcStateParameterSize_ += estimateSingleArcInitialStateParameters_[ i ]->getParameterSize( );
            initialSingleArcStateParameters_[ estimatedParameterSetSize_ ] = estimateSingleArcInitialStateParameters_[ i ];

            initialDynamicalStateParameterSize_ += estimateSingleArcInitialStateParameters_[ i ]->getParameterSize( );
            estimatedParameterSetSize_ += estimateSingleArcInitialStateParameters_[ i ]->getParameterSize( );
        }


        for( unsigned int i = 0; i < estimateMultiArcInitialStateParameters_.size( ); i++ )
        {
            initialStateParameters_[ estimatedParameterSetSize_ ] = estimateMultiArcInitialStateParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( ) ) );

            initialDynamicalMultiArcStateParameterSize_ += estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( );
            initialMultiArcStateParameters_[ estimatedParameterSetSize_ ] = estimateMultiArcInitialStateParameters_[ i ];

            initialDynamicalStateParameterSize_ += estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( );
            estimatedParameterSetSize_ += estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( );
        }


        // Iterate over all double parameters and add to parameter size and set indices in parameterIndices_
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            doubleParameters_[ estimatedParameterSetSize_ ] = estimatedDoubleParameters_[ i ];
            totalConstraintSize_ += estimatedDoubleParameters_[ i ]->getConstraintSize( );

            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_, 1 ) );
            estimatedParameterSetSize_++;
        }

        // Iterate over all vector parameter, add to total number of parameters and set indices in parameterIndices_
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            vectorParameters_[ estimatedParameterSetSize_ ] = estimatedVectorParameters_[ i ];
            totalConstraintSize_ += estimatedVectorParameters_[ i ]->getConstraintSize( );

            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimatedVectorParameters_[ i ]->getParameterSize( ) ) );
            estimatedParameterSetSize_ += estimatedVectorParameters_[ i ]->getParameterSize( );
        }

        totalParameterSetSize_ = estimatedParameterSetSize_;

        // Process multi-arc parameters
        multiArcInitialStateParametersPerArc_ = getMultiArcDynamicalStateToEstimatePerArc(
                estimateMultiArcInitialStateParameters_, bodiesToEstimatePerArc_, multiArcStateParametersSizePerArc_ );
    }

    //! Function to return the total number of parameter values (including consider parameters)
    /*!
     *  Function to return the total number of parameter values (including consider parameters)
     *  \return Size of parameter vector (including consider parameters)
     */
    int getParameterSetSize( )
    {
        return totalParameterSetSize_;
    }

    //! Function to return the total number of parameter values (excluding consider parameters).
    /*!
     *  Function to return the total number of parameter values (excluding consider parameters)
     *  \return Size of parameter vector (excluding consider parameters)
     */
    int getEstimatedParameterSetSize( )
    {
        return estimatedParameterSetSize_;
    }

    //! Function to return the total number of initial state values that are estimated.
    /*!
     *  Function to return the total number of initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalStateParameterSize( )
    {
        return initialDynamicalStateParameterSize_;
    }

    //! Function to return the total number of single-arc initial state values that are estimated.
    /*!
     *  Function to return the total number of single-arc initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalSingleArcStateParameterSize( )
    {
        return initialDynamicalSingleArcStateParameterSize_;
    }

    //! Function to return the total number of multi-arc initial state values that are estimated.
    /*!
     *  Function to return the total number of multi-arc initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalMultiArcStateParameterSize( )
    {
        return initialDynamicalMultiArcStateParameterSize_;
    }

    //! Function to return the total number of non-dynamical state parameters that are estimated.
    /*!
     * Function to return the total number of non-dynamical state parameters that are estimated.
     * \return Total number of estimated parameters that are NOT initial dynamical state parameters
     */
    int getNonDynamicalStateParameterSize( )
    {
        return estimatedParameterSetSize_ - initialDynamicalStateParameterSize_;
    }

    //! Function that returns a vector containing all current parameter values
    /*!
     *  Function that returns a vector containing all current parameter values. The total vector starts with the initial
     *  state parameters, followed by the double and vector parameters, respectively.
     *  Initial state, double and vector parameter values are concatenated in the order in which they are set in the
     *  estimateInitialStateParameters_, doubleParameters_ and vectorParameters_ members.
     *  \return Vector containing all parameter values
     */
    template< typename ParameterScalar >
    Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 > getFullParameterValues( )
    {
        Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >  parameterValues =
                Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >::Zero( totalParameterSetSize_ );

        int currentStartIndex = 0;

        // Retrieve initial state parameter values.
        for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
        {
            parameterValues.segment( currentStartIndex, estimateInitialStateParameters_[ i ]->getParameterSize( ) ) =
                    estimateInitialStateParameters_[ i ]->getParameterValue( ).template cast< ParameterScalar >( );
            currentStartIndex += estimateInitialStateParameters_[ i ]->getParameterSize( );
        }

        // Retrieve double parameter values.
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            parameterValues( currentStartIndex ) = static_cast< ParameterScalar >(
                        estimatedDoubleParameters_[ i ]->getParameterValue( ) );
            currentStartIndex++;
        }

        // Retrieve vector parameter values.
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            parameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ) =
                    estimatedVectorParameters_[ i ]->getParameterValue( ).template cast< ParameterScalar >( );
            currentStartIndex += estimatedVectorParameters_[ i ]->getParameterSize( );
        }

        return parameterValues;
    }

    //! Function to reset all parameter values.
    /*!
     *  Function to reset all parameter values.
     *  \param newParameterValues New parameter values. Order of values in vector must be same order as return vector of getFullParameterValues
     */
    template< typename ParameterScalar >
    void resetParameterValues( const Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >& newParameterValues )
    {
        // Check input consistency
        if( newParameterValues.rows( ) != totalParameterSetSize_ )
        {
            throw std::runtime_error( "Error when resetting parameters of parameter set, given vector has size " +
                                      std::to_string( newParameterValues.rows( ) ) +
                                      ", while internal size is " + std::to_string( totalParameterSetSize_ ) );
        }
        else
        {
            int currentStartIndex = 0;

//            std::cout << "before reset initial state parameters" << "\n\n";
            for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
            {
                estimateInitialStateParameters_[ i ]->setParameterValue(
                            newParameterValues.segment( currentStartIndex, estimateInitialStateParameters_[ i ]->getParameterSize( ) ).
                            template cast< InitialStateParameterType >( ) );
                currentStartIndex += estimateInitialStateParameters_[ i ]->getParameterSize( );
            }
//            std::cout << "after reset initial state parameters" << "\n\n";

            // Set double parameter values.
//            std::cout << "before reset double parameters" << "\n\n";
            for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
            {
                estimatedDoubleParameters_[ i ]->setParameterValue( static_cast< double >( newParameterValues( currentStartIndex ) ) );
                currentStartIndex++;
            }
//            std::cout << "after reset double parameters" << "\n\n";

            // Set vector parameter values.
//            std::cout << "before reset vector parameters" << "\n\n";
//            std::cout << "new parameters values: " << newParameterValues.transpose( ) << "\n\n";
            for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
            {
//                std::cout << "current start index: " << currentStartIndex << "\n\n";
//                std::cout << "current parameter size: " << estimatedVectorParameters_[ i ]->getParameterSize( ) << "\n\n";
//                std::cout << "newParameterValues current parameter: "
//                          << newParameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ).transpose( ) << "\n\n";
                estimatedVectorParameters_[ i ]->setParameterValue(
                            newParameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ).
                            template cast< double >( ) );

                currentStartIndex += estimatedVectorParameters_[ i ]->getParameterSize( );
            }
//            std::cout << "after reset vector parameters" << "\n\n";
        }
    }

    //! Function to retrieve double parameter objects.
    /*!
     *  Function to retrieve double parameter objects.
     *  \return Map containing all double parameter objects, with map key start index of parameter in total vector.
     */
    std::map< int, std::shared_ptr< EstimatableParameter< double > > > getDoubleParameters( )
    {
        return doubleParameters_;
    }

    //! Function to retrieve vector parameter objects.
    /*!
     *  Function to retrieve vector parameter objects.
     *  \return Map containing all vector parameter objects, with map key start index of parameter in total vector.
     */
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getVectorParameters( )
    {
        return vectorParameters_;
    }

    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getInitialStateParameters( )
    {
        return initialStateParameters_;
    }

    //! Function to retrieve all single-arc initial state parameter objects.
    /*!
     *  Function to retrieve all single-arc initial state parameter objects.
     *  \return Map containing all single-arc initial state parameter objects, with map key start index of parameter in total
     *  vector.
     */
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getInitialSingleArcStateParameters( )
    {
        return initialSingleArcStateParameters_;
    }

    //! Function to retrieve all multi-arc initial state parameter objects.
    /*!
     *  Function to retrieve all multi-arc initial state parameter objects.
     *  \return Map containing all multi-arc initial state parameter objects, with map key start index of parameter in total
     * vector.
     */
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getInitialMultiArcStateParameters( )
    {
        return initialMultiArcStateParameters_;
    }

    std::vector< std::shared_ptr< EstimatableParameter< double > > > getEstimatedDoubleParameters( )
    {
        return estimatedDoubleParameters_;
    }

    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getEstimatedVectorParameters( )
    {
        return estimatedVectorParameters_;
    }

    //! Function to retrieve the start index and size of (a) parameters(s) with a given identifier
    /*!
     * Function to retrieve the start index and size of (a) parameters(s) with a given identifier
     * \param requiredParameterId Parameter identifier that is to be searched in full list of patameters
     * \return List of start indices and sizes of parameters corresponding to requiredParameterId
     */
    std::vector< std::pair< int, int > > getIndicesForParameterType(
            const EstimatebleParameterIdentifier requiredParameterId )
    {
        std::vector< std::pair< int, int > > typeIndices;

        for( auto parameterIterator : initialSingleArcStateParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        for( auto parameterIterator : initialMultiArcStateParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        for( auto parameterIterator : doubleParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        for( auto parameterIterator : vectorParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        return typeIndices;
    }

    //! Function to get list of initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedInitialStateParameters( )
    {
        return estimateInitialStateParameters_;
    }

    //! Function to get list of single-arc initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of single-arc initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedSingleArcInitialStateParameters( )
    {
        return estimateSingleArcInitialStateParameters_;
    }

    //! Function to get list of multi-arc initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of multi-arc initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedMultiArcInitialStateParameters( )
    {
        return estimateMultiArcInitialStateParameters_;
    }

    //! Function to retrieve list of start indices and sizes (map keys) of estimated parameters.
    /*!
     *  Function to retrieve list of start indices and sizes (map keys) of estimated parameters.
     *  \return List of start indices and sizes (map keys) of estimated parameters.
     */
    std::vector< std::pair< int, int > > getParametersIndices( )
    {
        return parameterIndices_;
    }

    //! Function to retrieve total multiplier and right-hand side for parameter estimation linear constraint
    /*!
     * Function to retrieve total multiplier and right-hand side for parameter estimation linear constraint
     * \param constraintStateMultiplier Multiplier for parameter linear constraint
     * \param constraintRightHandSide Right-hand side for parameter linear constraint
     */
    void getConstraints( Eigen::MatrixXd& constraintStateMultiplier, Eigen::VectorXd& constraintRightHandSide )
    {
        // Resize constraint elements
        constraintStateMultiplier.setZero( totalConstraintSize_, estimatedParameterSetSize_ );
        constraintRightHandSide.setZero( totalConstraintSize_, 1 );

        // Iterate over all state parameters
        int currentConstraintRow = 0;
        int currentConstraintSize = 0;
        for( auto parameterIterator = initialStateParameters_.begin( ); parameterIterator != initialStateParameters_.end( );
             parameterIterator++ )
        {
            // Add constraint if of non-zero size
            currentConstraintSize = parameterIterator->second->getConstraintSize( );
            if( currentConstraintSize > 0 )
            {
                constraintStateMultiplier.block(
                            currentConstraintRow, parameterIterator->first, currentConstraintSize,
                            parameterIterator->second->getParameterSize( )  ) =
                        parameterIterator->second->getConstraintStateMultipler( );
                constraintRightHandSide.segment( currentConstraintRow, currentConstraintSize ) =
                        parameterIterator->second->getConstraintRightHandSide( );

                currentConstraintRow += currentConstraintSize;
            }

        }

        // Iterate over all double parameters
        for( auto parameterIterator = doubleParameters_.begin( ); parameterIterator != doubleParameters_.end( );
             parameterIterator++ )
        {
            // Add constraint if of non-zero size
            currentConstraintSize = parameterIterator->second->getConstraintSize( );
            if( currentConstraintSize > 0 )
            {
                constraintStateMultiplier.block(
                            currentConstraintRow, parameterIterator->first, currentConstraintSize,
                            parameterIterator->second->getParameterSize( )  ) =
                        parameterIterator->second->getConstraintStateMultipler( );
                constraintRightHandSide.segment( currentConstraintRow, currentConstraintSize ) =
                        parameterIterator->second->getConstraintRightHandSide( );

                currentConstraintRow += currentConstraintSize;
            }
        }

        // Iterate over all vector parameters
        for( auto parameterIterator = vectorParameters_.begin( ); parameterIterator != vectorParameters_.end( );
             parameterIterator++ )
        {
            // Add constraint if of non-zero size
            currentConstraintSize = parameterIterator->second->getConstraintSize( );
            if( currentConstraintSize > 0 )
            {
                constraintStateMultiplier.block(
                            currentConstraintRow, parameterIterator->first, currentConstraintSize,
                            parameterIterator->second->getParameterSize( )  ) =
                        parameterIterator->second->getConstraintStateMultipler( );
                constraintRightHandSide.segment( currentConstraintRow, currentConstraintSize ) =
                        parameterIterator->second->getConstraintRightHandSide( );

                currentConstraintRow += currentConstraintSize;
            }
        }
    }

    //! Total size of linear constraint that is to be applied during estimation
    /*!
     * Total size of linear constraint that is to be applied during estimation
     * \return Size of linear constraint that is to be applied during estimation
     */
    int getConstraintSize( )
    {
        return totalConstraintSize_;
    }

    //! Function to return a map with the names of the arc-wise bodies to be estimated for each arc, and the arc starting times as keys.
    /*!
     * Function to return a map with the names of the arc-wise bodies to be estimated for each arc, and the arc starting times as keys.
     * \return Map with names of the arc-wise bodies to be estimated, for each arc.
     */
    std::map< double, std::vector< std::string > > getBodiesToEstimatePerArc( )
    {
        return bodiesToEstimatePerArc_;
    }

    //! Function to return the list of multi-arc dynamical state parameters sizes (arc-wise).
    /*!
     * Function to return the list of multi-arc dynamical state parameters sizes (arc-wise).
     * \return Vector with multi-arc state parameters sizes.
     */
    std::vector< int > getMultiArcStateParametersSizePerArc( )
    {
        return multiArcStateParametersSizePerArc_;
    }

    //! Function to return a map with the multi-arc state parameters to estimate, ordered per arc (with the map keys containing the arc starting times)
    /*!
     * Function to return a map with the multi-arc state parameters to estimate, ordered per arc (with the map keys containing the arc starting times)
     * \return Map with multi-arc state parameters to be estimated, for each arc (key = arc starting time).
     */
    std::map< double, std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > >
            getMultiArcInitialStateParametersPerArc( )
    {
        return multiArcInitialStateParametersPerArc_;
    }

    //! Function to return the number of arcs detected from the multi-arc parameters to be estimated.
    /*!
     * Function to return the number of arcs detected from the multi-arc parameters to b estimated.
     * \return Number of arcs detected from the multi-arc parameters to estimated.
     */
     int getNumberArcsFromMultiArcInitialStateParameters( )
    {
         return bodiesToEstimatePerArc_.size( );
    }

    //! Function to return the arc starting times detected from the multi-arc parameters to be estimated.
    /*!
     * Function to return the arc starting times detected from the multi-arc parameters to be estimated.
     * \return Vector containing the arc starting times detected from the multi-arc dynamical parameters.
     */
     std::vector< double > getArcStartingTimes( )
    {
         return utilities::createVectorFromMapKeys( bodiesToEstimatePerArc_ );
    }

protected:

    //! Total size of all initial dynamical states that are to be estimated.
    int initialDynamicalStateParameterSize_;

    //! Total size of all initial single-arc dynamical states that are to be estimated.
    int initialDynamicalSingleArcStateParameterSize_;

    //! Total size of all initial multi-arc dynamical states that are to be estimated.
    int initialDynamicalMultiArcStateParameterSize_;

    //! Total number of parameter values (including currently non yet implemented consider parameters).
    int totalParameterSetSize_;

    //! Total number of estimated parameter values (excluding currently non yet implemented consider parameters).
    int estimatedParameterSetSize_;

    //! List of start indices and sizes (map keys) of estimated parameters.
    /*!
     * List of start indices and sizes (map keys) of estimated parameters, in order of vector
     * estimateInitialStateParameters_, followed by estimatedDoubleParameters_, followed by estimatedVectorParameters_.
     */
    std::vector< std::pair< int, int > > parameterIndices_;

    //! List of double parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< double > > > estimatedDoubleParameters_;

    //! List of vector parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedVectorParameters_;

    //! List of initial dynamical states that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateInitialStateParameters_;

    //! List of initial single-arc dynamical states that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateSingleArcInitialStateParameters_;

    //! List of initial multi-arc dynamical states that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateMultiArcInitialStateParameters_;

    //! Map of double parameters that are to be estimated, with start index in total parameter vector as key.
    std::map< int, std::shared_ptr< EstimatableParameter< double > > > doubleParameters_;

    //! Map of vector parameters that are to be estimated, with start index in total parameter vector as key.
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters_;

    //! Map of initial dynamical states that are to be estimated, with start index in total parameter vector as key.
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialStateParameters_;

    //! Size of linear constraint that is to be applied during estimation
    int totalConstraintSize_;

    //! Map containing all single-arc initial state parameter objects, with map key start index of parameter in total vector.
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialSingleArcStateParameters_;

    //! Map containing all multi-arc initial state parameter objects, with map key start index of parameter in total vector.
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialMultiArcStateParameters_;

    //! Map containing the names of the multi-arc bodies to be estimated for each arc, with map key being the arc starting time.
    std::map< double, std::vector< std::string > > bodiesToEstimatePerArc_;

    //! List of multi-arc dynamical state parameters sizes (arc-wise).
    std::vector< int > multiArcStateParametersSizePerArc_;

    //! List of initial multi-arc dynamical states that are to be estimated, for each arc.
    std::map< double, std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > >
    multiArcInitialStateParametersPerArc_;

};

//! Function to create a subset of all estimated parameters, with either only single-arc or multi-arc initial state parameter
/*!
 *  Function to create a subset of all estimated parameters, with either only single-arc or multi-arc initial state parameter.
 *  All non-dynamical state parameters are copied from the input to the output
 *  \param parametersToEstimate Total set of parameters
 *  \param getSingleArcParameters Boolean denoting whether single- or multi-arc initial state parameters are to be kept
 *  \return Subset of all estimated parameters, with either only single-arc or multi-arc initial state parameter
 */
template< typename StateScalarType >
std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > createEstimatableParameterSetArcSubSet(
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const bool getSingleArcParameters)
{
    if( getSingleArcParameters )
    {
        return std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >(
                    parametersToEstimate->getEstimatedDoubleParameters( ),
                    parametersToEstimate->getEstimatedVectorParameters( ),
                    parametersToEstimate->getEstimatedSingleArcInitialStateParameters( ) );
    }
    else
    {
        return std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >(
                    parametersToEstimate->getEstimatedDoubleParameters( ),
                    parametersToEstimate->getEstimatedVectorParameters( ),
                    parametersToEstimate->getEstimatedMultiArcInitialStateParameters( )  );
    }
}

template< typename InitialStateParameterType >
void printEstimatableParameterEntries(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::map< int, std::shared_ptr<
            EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialStateParameters =
            estimatableParameters->getInitialStateParameters( );
    std::map< int, std::shared_ptr<
            EstimatableParameter< double > > > doubleParameters = estimatableParameters->getDoubleParameters( );
    std::map< int, std::shared_ptr<
            EstimatableParameter< Eigen::VectorXd > > > vectorParameters = estimatableParameters->getVectorParameters( );

    std::cout << "Parameter start index, Parameter definition" << std::endl;
    for( typename  std::map< int, std::shared_ptr<  EstimatableParameter< Eigen::Matrix<
         InitialStateParameterType, Eigen::Dynamic, 1 > > > >::const_iterator parameterIterator = initialStateParameters.begin( );
         parameterIterator != initialStateParameters.end( ); parameterIterator++ )
    {
        std::cout << parameterIterator->first << ", " << parameterIterator->second->getParameterDescription( ) << std::endl;
    }

    for( typename  std::map< int, std::shared_ptr<  EstimatableParameter< double > > >::const_iterator
         parameterIterator = doubleParameters.begin( );
         parameterIterator != doubleParameters.end( ); parameterIterator++ )
    {
        std::cout << parameterIterator->first << ", " << parameterIterator->second->getParameterDescription( ) << std::endl;
    }

    for( typename  std::map< int, std::shared_ptr<  EstimatableParameter< Eigen::VectorXd > > >::const_iterator
         parameterIterator = vectorParameters.begin( );
         parameterIterator != vectorParameters.end( ); parameterIterator++ )
    {
        std::cout << parameterIterator->first << ", " << parameterIterator->second->getParameterDescription( ) << std::endl;
    }
    std::cout << std::endl;
}

//! Function to get the list of names of bodies for which initial translational dynamical state is estimated.
/*!
 *  Function to get the list of names of bodies for which initial translational dynamical state is estimated.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return List of names of bodies for which initial state is estimated.
 */
template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithTranslationalStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

//! Function to get the list of names of bodies for which initial rotational dynamical state is estimated.
/*!
 *  Function to get the list of names of bodies for which initial rotational dynamical state is estimated.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return List of names of bodies for which initial rotational state is estimated.
 */
template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithRotationalStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithMassStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_mass_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

//! Function to retrieve the list of bodies for which the translational state is estimated in a multi-arc fashion
/*!
 * Function to retrieve the list of bodies for which the translational state is estimated in a multi-arc fashion
 * \param estimatableParameters Full set of estimated parameters
 * \return List of parameters (with body names as keys) used for the multi-arc estimation of initial translational state
 */
template< typename InitialStateParameterType >
std::map< std::string, std::shared_ptr< EstimatableParameter<
Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  > > >
getListOfBodiesWithTranslationalMultiArcStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::map< std::string, std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  > > > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state )
        {
            bodiesToEstimate[ initialDynamicalParameters.at( i )->getParameterName( ).second.first ] =
                    initialDynamicalParameters.at( i );
        }
    }

    return bodiesToEstimate;
}

//! Function to retrieve the list of bodies for which an initial dynamical state is to be estimated
/*!
 * Function to retrieve the list of bodies for which an initial dynamical state is to be estimated
 * \param estimatableParameters Full set of estimated parameters
 * \return List of bodies for which an initial dynamical state is estimated.
 */
template< typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, std::vector< std::string > > getListOfBodiesToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::map< propagators::IntegratedStateType, std::vector< std::string > > bodiesToEstimate;

    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )  ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state ) )
        {
            bodiesToEstimate[ propagators::translational_state ].push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
        else if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state ) )
         {
             bodiesToEstimate[ propagators::rotational_state ].push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
         }
    }

    return bodiesToEstimate;
}

//! Function to retrieve the list of translational state parameters from full parameter list
/*!
 * Function to retrieve the list of translational state parameters from full parameter list
 * \param estimatableParameters Full set of estimated parameters
 * \return List of translational state parameters
 */
template< typename InitialStateParameterType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getListOfTranslationalStateParametersToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > translationalStateParameters;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )  ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state ) )
        {
            translationalStateParameters.push_back(  initialDynamicalParameters.at( i ) );
        }
    }

    return translationalStateParameters;
}

//! Function to retrieve the list of rotational state parameters from full parameter list
/*!
 * Function to retrieve the list of rotational state parameters from full parameter list
 * \param estimatableParameters Full set of estimated parameters
 * \return List of rotational state parameters
 */
template< typename InitialStateParameterType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getListOfRotationalStateParametersToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > rotationalStateParameters;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state )   )
        {
            rotationalStateParameters.push_back(  initialDynamicalParameters.at( i ) );
        }
    }

    return rotationalStateParameters;
}

template< typename InitialStateParameterType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getListOfMassStateParametersToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > massStateParameters;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_mass_state )   )
        {
            massStateParameters.push_back(  initialDynamicalParameters.at( i ) );
        }
    }

    return massStateParameters;
}

//! Function to get the complete list of initial dynamical states that are to be estimated, sorted by dynamics type.
/*!
 *  Function to get the complete list of initial dynamical states that are to be estimated, sorted by dynamics type.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return Map containing dynamics type (key) and vector of pairs: list of bodies (first in pair) with reference point
 *  identifier (second in pair; empty if not relevant) for which given dynamics type is estimated.
 */
template< typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
getListOfInitialDynamicalStateParametersEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > initialDynamicalStateParametersEstimate;
    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state ) ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state ) )
        {
            initialDynamicalStateParametersEstimate[ propagators::translational_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }
        else if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state ) )
        {
            initialDynamicalStateParametersEstimate[ propagators::rotational_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }
        else if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_mass_state ) )
        {
            initialDynamicalStateParametersEstimate[ propagators::body_mass_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }
    }

    return initialDynamicalStateParametersEstimate;
}

//! Function to get initial state vector of estimated dynamical states.
/*!
 *  Function to get initial state vector of estimated dynamical states (i.e. presently estimated state at propagation
 *  start time.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return State vector of estimated dynamics at propagation start time.
 */
template< typename InitialStateParameterType = double >
Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getInitialStateVectorOfBodiesToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Initialize state vector.
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateVector =
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >::Zero(
                estimatableParameters->getInitialDynamicalStateParameterSize( ), 1 );

    int vectorSize = 0;
    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( isParameterDynamicalPropertyInitialState( initialDynamicalParameters.at( i )->getParameterName( ).first ) )
        {
            int currentParameterSize = initialDynamicalParameters.at( i )->getParameterSize( );
            initialStateVector.block( vectorSize, 0, currentParameterSize, 1 ) = initialDynamicalParameters.at( i )->getParameterValue( );

            vectorSize += currentParameterSize;
        }
    }

    return initialStateVector.block( 0, 0, vectorSize, 1 );
}

//! Function to get the complete list of multi-arc dynamical states that are to be estimated, sorted by arc start times.
/*!
 *  Function to get the complete list of multi-arc dynamical states that are to be estimated, sorted by arc start times.
 *  \param initialMultiArcStateParameters Object containing all parameters that are to be estimated.
 *  \param bodiesToEstimatePerArc List of arc starting times (keys) and bodies to be estimated for each arc (to be returned by reference).
 *  \param multiArcStateParametersSizePerArc List of multi-arc state parameters to be estimated for each each (to be returned by reference).
 *  \return Map containing arc starting times (key) and vector of multi-arc dynamical state parameters to estimate
 */
template< typename InitialStateParameterType >
std::map< double, std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > >
        getMultiArcDynamicalStateToEstimatePerArc(
                const std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >&
                        initialMultiArcStateParameters,
                std::map< double, std::vector< std::string > >& bodiesToEstimatePerArc,
                std::vector< int >& multiArcStateParametersSizePerArc )
{
    bodiesToEstimatePerArc.clear( );
    multiArcStateParametersSizePerArc.clear( );
    std::map< double, std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > >
            multiArcDynamicalStateParametersPerArc;
    std::map< double, int > arcWiseParameterSize;

    for ( unsigned int j = 0 ; j < initialMultiArcStateParameters.size( ) ; j++ )
    {
        switch ( initialMultiArcStateParameters.at( j )->getParameterName( ).first )
        {
            case arc_wise_initial_body_state:
            {
                std::shared_ptr< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > > arcWiseTranslationalStateParameter =
                        std::dynamic_pointer_cast< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > >( initialMultiArcStateParameters.at( j ) );
                std::vector< double > currentArcTimes = arcWiseTranslationalStateParameter->getArcStartTimes( );

                for ( unsigned int i = 0 ; i < currentArcTimes.size( ) ; i++ )
                {
                    // Multi-arc state parameter limited to current arc.
                    std::vector< double > currentArcStartTime = { currentArcTimes.at( i ) };
                    std::vector< std::string > currentArcCentralBody = { arcWiseTranslationalStateParameter->getCentralBodies( ).at( i ) };
                    std::shared_ptr< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > > currentArcTranslationalStateParameter =
                            std::make_shared< ArcWiseInitialTranslationalStateParameter< InitialStateParameterType > >(
                                    arcWiseTranslationalStateParameter->getParameterName( ).second.first,
                                    currentArcStartTime, arcWiseTranslationalStateParameter->getParameterValue( ).segment( i * 6, 6 ),
                                    currentArcCentralBody, arcWiseTranslationalStateParameter->getFrameOrientation( ) );

                    // Check whether the current arc was already detected for other multi-arc initial state parameters
                    // and add the current parameter to the per-arc list.
                    bool alreadyDetectedArc = false;
                    for ( auto itr: multiArcDynamicalStateParametersPerArc )
                    {
                        if( std::fabs( currentArcTimes.at( i ) - itr.first ) <
                            std::max( 4.0 * itr.first * std::numeric_limits< double >::epsilon( ), 1.0E-12 ) )
                        {
                            bodiesToEstimatePerArc[ itr.first ].push_back( arcWiseTranslationalStateParameter->getParameterName( ).second.first );
                            itr.second.push_back( currentArcTranslationalStateParameter );
                            arcWiseParameterSize[ itr.first ] += 6;
                            alreadyDetectedArc = true;
                        }
                    }
                    // Create a new entry in the per-arc list of  multi-arc initial state parameters if the current arc is detected
                    // for the first time.
                    if ( !alreadyDetectedArc )
                    {
                        bodiesToEstimatePerArc[ currentArcTimes.at( i ) ] = { arcWiseTranslationalStateParameter->getParameterName( ).second.first };
                        multiArcDynamicalStateParametersPerArc[ currentArcTimes.at( i ) ] = { currentArcTranslationalStateParameter };
                        arcWiseParameterSize[ currentArcTimes.at( i ) ] = 6;
                    }
                }

                // Compute size of multi-arc state parameters to be estimated for each arc
                multiArcStateParametersSizePerArc = utilities::createVectorFromMapValues( arcWiseParameterSize );
                break;
            }
            default:
                throw std::runtime_error( "Multi-arc dynamical parameter not recognised (only multi-arc translational initial state currently "
                                          " implemented.)" );
        }
    }

    return multiArcDynamicalStateParametersPerArc;
}

//! Function to retrieve the size of the estimatable parameter set.
/*!
 *  Function to retrieve the size of the estimatable parameter set.
 *  \param estimatableParameterSet Set of estimatable parameters.
 *  \return Size of parameter set.
 */
template< typename InitialStateParameterType = double >
int getSingleArcParameterSetSize(
        std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameterSet,
        const int arcIndex = 0 )
{
    // Check the consistency of the arc index
    int detectedNumberArcs = estimatableParameterSet->getNumberArcsFromMultiArcInitialStateParameters( );
    if ( ( ( arcIndex < 0 ) || ( arcIndex + 1 > detectedNumberArcs ) ) && ( detectedNumberArcs > 0 ) )
    {
//        std::cout << "detected number arcs: " << detectedNumberArcs << "\n\n";
//        std::cout << "arc index: " << arcIndex << "\n\n";
        throw std::runtime_error( "Error when getting single arc parameter size for a given arc, the required arc index is unconsistent with "
                                  "the detected number of arcs." );
    }


    // Check consistency between number of detected arcs and size of the vector containing the multi-arc dynamical state sizes.
    if ( static_cast< unsigned int >( detectedNumberArcs ) != estimatableParameterSet->getMultiArcStateParametersSizePerArc( ).size( ) )
    {
        throw std::runtime_error(  "Error when getting single arc parameter size for a given arc, inconsistency between the detected number of arcs and "
                                   "the vector giving the multi-arc dynamical state parameter sizes per arc." );
    }

    // Compute parameter size for the given arc.
    if ( detectedNumberArcs > 0 )
    {
        return estimatableParameterSet->getInitialDynamicalSingleArcStateParameterSize( ) +
               estimatableParameterSet->getNonDynamicalStateParameterSize( ) +
               estimatableParameterSet->getMultiArcStateParametersSizePerArc( ).at(arcIndex);
    }
    else
    {
        return estimatableParameterSet->getInitialDynamicalSingleArcStateParameterSize( ) +
               estimatableParameterSet->getNonDynamicalStateParameterSize( );
    }
}

//! Function to retrieve the size of the dynamical state.
/*!
 *  Function to retrieve the size of the dynamical state.
 *  \param estimatableParameterSet Set of estimatable parameters.
 *  \return Size of the initial dynamical state.
 */
template< typename InitialStateParameterType = double >
int getSingleArcInitialDynamicalStateParameterSetSize(
        std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameterSet,
        const int arcIndex = 0 )
{
    return getSingleArcParameterSetSize( estimatableParameterSet, arcIndex ) -
            estimatableParameterSet->getNonDynamicalStateParameterSize( );
}

//! Function to get arc start times from list of estimated parameters
/*!
 *  Function to get arc start times from list of estimated parameters. Function throws an error if multiple arc-wise
 *  estimations are found, but arc times are not compatible
 *  \param estimatableParameters List of estimated parameters
 *  \param throwErrorOnSingleArcDynamics Boolean denoting whether to throw an exception if single arc dynamics are used (default true)
 *  \return Start times for estimation arcs
 */
template< typename InitialStateParameterType >
std::vector< double > getMultiArcStateEstimationArcStartTimes(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters,
        const bool throwErrorOnSingleArcDynamics = true )

{
    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
                                  Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
    estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required and check if multi-arc dynamics are detected.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first != arc_wise_initial_body_state )
        {
            if( throwErrorOnSingleArcDynamics )
            {
                throw std::runtime_error( "Error when getting arc times from estimated parameters, single arc dynamics found" );
            }
        }
    }

    return estimatableParameters->getArcStartingTimes( );
}

///! Retrieve parameters to be estimated for each arc (arc-wise parameters might differ from one arc to another).
/*!
 * Retrieve parameters to be estimated for each arc (arc-wise parameters might differ from one arc to another).
 * \param parametersToEstimate Pointer for estimated parameters, provided as input of the whole multi-arc variational equations solver.
 * \param arcWiseParametersToEstimate Vector containing the estimated parameters for each arc (returned by reference).
 * \param estimatedBodiesPerArc list of bodies to be estimated, for each arc.
 */
template< typename StateScalarType = double >
void getParametersToEstimatePerArcTest(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > >& arcWiseParametersToEstimate,
        const std::vector< double >& arcStartTimes,
        std::map< int, std::vector< std::string > >& estimatedBodiesPerArc,
        std::map< int, std::map< std::string, int > >& arcIndicesPerBody )
{

    arcWiseParametersToEstimate.clear( );

    // Get list of objets and associated bodies to estimate initial arc-wise translational states
    typedef std::map< std::string, std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > ArcWiseParameterList;
    ArcWiseParameterList estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalMultiArcStateToEstimate(
            parametersToEstimate );

    int numberEstimatedBodies = estimatedBodies.size( );

    // Check that the arc starting times are provided in correct order.
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) - 1 ; i++ )
    {
        if ( ( arcStartTimes[ i + 1 ] - arcStartTimes[ i ] ) < 0.0 )
        {
            throw std::runtime_error( "Error, arc start times not provided in increasing order." );
        }
    }

    // Initialising vector keeping track of whether each arc is associated with at least one body whose multi-arc state is to be estimated.
    std::vector< bool > detectedEstimatedStatesPerArc;
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        detectedEstimatedStatesPerArc.push_back( false );
    }

    estimatedBodiesPerArc.clear( );
    arcIndicesPerBody.clear( );
    std::vector< int > counterStateIndicesPerBody;
    for ( int i = 0 ; i < numberEstimatedBodies ; i++ )
    {
        counterStateIndicesPerBody.push_back( 0 );
    }

    // Iterate over all parameters and check consistency
    unsigned int counterEstimatedBody = 0;
    for( typename ArcWiseParameterList::const_iterator parameterIterator = estimatedBodies.begin( ); parameterIterator !=
                                                                                                     estimatedBodies.end( ); parameterIterator++ )
    {
        // Get arc start times of current parameter
        std::vector< double > parameterArcStartTimes =
                std::dynamic_pointer_cast< estimatable_parameters::
                ArcWiseInitialTranslationalStateParameter< StateScalarType > >(
                        parameterIterator->second )->getArcStartTimes( );

        // Check that each arc has at least one body whose state is to be estimated.
        for ( unsigned int i = 0 ; i < parameterArcStartTimes.size( ) ; i++ )
        {
//            bool detectedArc = false;
            int indexDetectedArc = 0;
            for ( unsigned int j = indexDetectedArc ; j < arcStartTimes.size( ) ; j++ )
            {
                if( std::fabs( arcStartTimes.at( j ) - parameterArcStartTimes.at( i ) ) <
                    std::max( 4.0 * parameterArcStartTimes.at( i ) * std::numeric_limits< double >::epsilon( ), 1.0E-12 ) )
                {
//                    detectedArc = true;
                    indexDetectedArc = j;
                    detectedEstimatedStatesPerArc[ j ] = true;

                    estimatedBodiesPerArc[ indexDetectedArc ].push_back( parameterIterator->first );
                    arcIndicesPerBody[ indexDetectedArc ][ parameterIterator->first ] = counterStateIndicesPerBody[ counterEstimatedBody ];
                    counterStateIndicesPerBody[ counterEstimatedBody ] += 1;
                }
            }
        }

        counterEstimatedBody += 1;
    }

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
            initialStatesParameters = parametersToEstimate->getEstimatedInitialStateParameters( );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
            parametersToEstimate->getEstimatedDoubleParameters( );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parametersToEstimate->getEstimatedVectorParameters( );

    for ( unsigned int arc = 0 ; arc < estimatedBodiesPerArc.size( ) ; arc++ )
    {
        std::vector< std::string > arcWiseBodiesToEstimate = estimatedBodiesPerArc.at( arc );

        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
                arcWiseStatesParameters;

        for ( unsigned int j = 0 ; j < initialStatesParameters.size( ) ; j++ )
        {
            for ( unsigned int body = 0 ; body < arcWiseBodiesToEstimate.size( ) ; body++ )
            {
                if ( arcWiseBodiesToEstimate[ body ] == initialStatesParameters[ j ]->getParameterName( ).second.first )
                {
                    //arcWiseStatesParameters.push_back( initialStatesParameters[ j ] );

                    std::shared_ptr< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< StateScalarType > > currentArcInitialStateParameter =
                            std::dynamic_pointer_cast< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< StateScalarType > >
                                    ( initialStatesParameters[ j ] );

                    if ( currentArcInitialStateParameter == nullptr )
                    {
                        throw std::runtime_error( "Error, initial state parameter type for multi-arc is not ArcWiseInitialTranslationalStateParameter." );
                    }
                    else
                    {
                        int arcIndexForBody = arcIndicesPerBody.at( arc ).at( arcWiseBodiesToEstimate[ body ] );
                        std::vector< double > currentArcStartTime = { currentArcInitialStateParameter->getArcStartTimes( )[ arcIndexForBody ] };
                        std::vector< std::string > currentArcCentralBody = { currentArcInitialStateParameter->getCentralBodies( )[ arcIndexForBody ] };

                        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentArcInitialState =
                                currentArcInitialStateParameter->getParameterValue( ).block( 6 * arcIndexForBody, 0, 6, 1 );

                        std::shared_ptr< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< StateScalarType > > arcWiseInitialStateParameter
                                = std::make_shared< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< StateScalarType > >
                                        ( arcWiseBodiesToEstimate[ body ], currentArcStartTime,
                                         currentArcInitialState, currentArcCentralBody,
                                         currentArcInitialStateParameter->getFrameOrientation( ) );

                        arcWiseStatesParameters.push_back( arcWiseInitialStateParameter );

//                        arcWiseStatesParameters.push_back( initialStatesParameters[ j ] );

                    }
                }
            }
        }

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > arcWiseEstimatableParamatersSet =
                std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >
                        ( doubleParameters, vectorParameters, arcWiseStatesParameters );

        arcWiseParametersToEstimate.push_back( arcWiseEstimatableParamatersSet );

    }
}


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_ESTIMATABLEPARAMETERSET_H
