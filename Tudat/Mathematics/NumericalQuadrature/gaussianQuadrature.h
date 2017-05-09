/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GAUSSIAN_QUADRATURE_H
#define TUDAT_GAUSSIAN_QUADRATURE_H

#include <vector>
#include <map>

#include <boost/function.hpp>
// #include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{

namespace numerical_quadrature
{


//! Read Gaussian nodes from text file
template< typename IndependentVariableType >
static void readGaussianQuadratureNodes(
        std::map< unsigned int, Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1> >& gaussQuadratureNodes )
{
   gaussQuadratureNodes =
           utilities::convertSTLVectorMapToEigenVectorMap< unsigned int, double >(
               input_output::readStlVectorMapFromFile< unsigned int, IndependentVariableType >(
                input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianNodes.txt" ) );
}

//! Read Gaussian weight factors from text file
template< typename IndependentVariableType >
static void readGaussianQuadratureWeights(
        std::map< unsigned int, Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1> >& gaussQuadratureWeights )
{
    gaussQuadratureWeights = utilities::convertSTLVectorMapToEigenVectorMap< unsigned int, double >(
                input_output::readStlVectorMapFromFile< unsigned int, IndependentVariableType >(
                input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianWeights.txt" ) );
}


//! Gaussian numerical quadrature wrapper class.
/*!
 * Numerical method that uses the Gaussian nodes and weight factors to compute definite integrals of a function.
 * The Gaussian nodes and weight factors are not calculated, but read from text files. The number of nodes (or
 * weight factors) has to be at least n = 2. The current text files contain tabulated values up to n = 64.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class GaussianQuadrature : public NumericalQuadrature< IndependentVariableType , DependentVariableType >
{
public:

    //! Empty constructor.
    GaussianQuadrature( ) { }


    //! Constructor.
    /*!
     * Constructor
     * \param integrand Function to be integrated numerically.
     * \param lowerLimit Lower limit for the integral.
     * \param upperLimit Upper limit for the integral.
     * \param numberOfNodes Number of nodes (i.e. nodes) at which the integrand will be evaluated.
     * Must be an integer value between 2 and 64.
     */
    GaussianQuadrature( const boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                        const IndependentVariableType lowerLimit, const IndependentVariableType upperLimit,
                        const unsigned int numberOfNodes ):
        integrand_ ( integrand ), lowerLimit_( lowerLimit ), upperLimit_ ( upperLimit ),
        numberOfNodes_( numberOfNodes ), quadratureHasBeenPerformed_( false )
    {
        readGaussianQuadratureNodes( uniqueNodes_ );
        readGaussianQuadratureWeights( uniqueWeights_ );
    }


    //! Reset the current Gaussian quadrature.
    /*!
     * The nodes and weights are not read/computed again if they had already been used previously.
     * \param integrand Function to be integrated numerically.
     * \param lowerLimit Lower limit for the integral.
     * \param upperLimit Upper limit for the integral.
     * \param numberOfNodes Number of nodes (i.e. nodes) at which the integrand will be evaluated.
     * Must be an integer value between 2 and 64.
     */
    void reset( const boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                const IndependentVariableType lowerLimit, const IndependentVariableType upperLimit,
                const unsigned int numberOfNodes )
    {
        integrand_ = integrand;
        lowerLimit_ = lowerLimit;
        upperLimit_ = upperLimit;
        numberOfNodes_ = numberOfNodes;
        quadratureHasBeenPerformed_ = false;
    }


    //! Function to return computed value of the quadrature.
    /*!
     *  Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     *  \return Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     */
    DependentVariableType getQuadrature( )
    {
        if ( ! quadratureHasBeenPerformed_ )
        {
            if ( integrand_.empty() )
            {
                throw std::runtime_error(
                            "The integrand for the Gaussian quadrature has not been set." );
            }

            if ( lowerLimit_ > upperLimit_ )
            {
                throw std::runtime_error(
                            "The lower limit for the Gaussian quadrature is larger than the upper limit." );
            }

            if ( numberOfNodes_ < 2 || numberOfNodes_ > 64 )
            {
                throw std::runtime_error(
                            "The number of nodes for the Gaussian quadrature must be between 2 and 64." );
            }

            performQuadrature();
            quadratureHasBeenPerformed_ = true;
        }

        return quadratureResult_;
    }


    typedef Eigen::Array< DependentVariableType, Eigen::Dynamic, 1 > DependentVariableArray;
    typedef Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1 > IndependentVariableArray;

    //! Get all the nodes (i.e. n nodes for nth order) from uniqueNodes_
    IndependentVariableArray getNodes( const unsigned int n )
    {
        if ( nodes_.count( n ) == 0 )
        {
            IndependentVariableArray newNodes( n );

            // Include node 0.0 if n is odd
            unsigned int i = 0;
            if ( n % 2 == 1 )
            {
                newNodes.row( i++ ) = 0.0;
            }

            // Include ± nodes
            IndependentVariableArray uniqueNodes_ = getUniqueNodes( n );
            for ( unsigned int j = 0; j < uniqueNodes_.size(); j++ )
            {
                newNodes.row( i++ ) = -uniqueNodes_[ j ];
                newNodes.row( i++ ) =  uniqueNodes_[ j ];
            }

            nodes_[ n ] = newNodes;
        }

        return nodes_.at( n );
    }

    //! Get all the weight factors (i.e. n weight factors for nth order) from uniqueWeights_
    IndependentVariableArray getWeights( const unsigned int n )
    {
        if ( weights_.count( n ) == 0 )
        {
            IndependentVariableArray newWeights( n );

            IndependentVariableArray orderNWeights = getUniqueWeights( n );
            // Include non-repeated weight factor if n is odd
            unsigned int i = 0;
            unsigned int j = 0;
            if ( n % 2 == 1 )
            {
                newWeights.row( i++ ) = orderNWeights[ j++ ];
            }

            // Include repeated weight factors
            for ( ; j < orderNWeights.size( ); j++ )
            {
                newWeights.row( i++ ) = orderNWeights[ j ];
                newWeights.row( i++ ) = orderNWeights[ j ];
            }

            weights_[ n ] = newWeights;
        }

        return weights_.at( n );
    }


protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     * Function that is called to perform the numerical quadrature. Sets the result in the quadratureResult local
     * variable.
     */
    void performQuadrature( )
    {
        // Determine the values of the auxiliary independent variable (nodes)
        const IndependentVariableArray nodes = getNodes( numberOfNodes_ );

        // Determine the values of the weight factors
        const IndependentVariableArray weights = getWeights( numberOfNodes_ );

        // Change of variable -> from range [-1, 1] to range [lowerLimit, upperLimit]
        const IndependentVariableArray independentVariables =
                0.5 * ( ( upperLimit_ - lowerLimit_ ) * nodes + upperLimit_ + lowerLimit_ );

        // Determine the value of the dependent variable
        DependentVariableArray weighedIntegrands( numberOfNodes_ );
        for ( unsigned int i = 0; i < numberOfNodes_; i++ )
        {
            weighedIntegrands( i ) = weights( i ) * integrand_( independentVariables( i ) );
        }

        quadratureResult_ = 0.5 * ( upperLimit_ - lowerLimit_ ) * weighedIntegrands.sum( );
    }


private:

    //! Map containing the nodes read from the text file (currently up to `n = 64`).
    //! The following relation holds: `size( uniqueNodes_[n] ) = floor( n / 2 )`
    //! For the actual nodes, the following must hold: `size( nodes[n] ) = n`
    //! The actual nodes are generated from `uniqueNodes_` by `getNodes()`
    std::map< unsigned int, IndependentVariableArray > uniqueNodes_;
    std::map< unsigned int, IndependentVariableArray > nodes_;

    //! Map containing the weight factors read from the text file (currently up to `n = 64`).
    //! The following relation holds: `size( uniqueWeights_[n] ) = ceil( n / 2 )`
    //! For the actual weight factors, the following must hold: `size( uniqueWeights_[n] ) = n`
    //! The actual weight factors are generated from `uniqueWeights_` by `getWeights()`
    std::map< unsigned int, IndependentVariableArray > uniqueWeights_;
    std::map< unsigned int, IndependentVariableArray > weights_;

    //! Function returning the integrand.
    boost::function< DependentVariableType( IndependentVariableType ) > integrand_;

    //! Lower limit for the integral.
    IndependentVariableType lowerLimit_;

    //! Upper limit for the integral.
    IndependentVariableType upperLimit_;

    //! Number of nodes.
    unsigned int numberOfNodes_;

    //! Whether quadratureResult has been set for the current integrand, lowerLimit, upperLimit and numberOfNodes.
    bool quadratureHasBeenPerformed_;

    //! Computed value of the quadrature, as computed by last call to performQuadrature.
    DependentVariableType quadratureResult_;

    //! Get the unique nodes for a specified order `n`.
    /*!
     * \param n The number of nodes or weight factors.
     * \return `uniqueNodes_[n]`, after reading the text file with the tabulated nodes if necessary.
     */
    IndependentVariableArray getUniqueNodes( const unsigned int n )
    {
        if ( uniqueNodes_.count( n ) == 0 )
        {
            std::string errorMessage = "Error in Gaussian quadrature, nodes not available for n=" +
                    boost::lexical_cast< std::string >( n );
            throw std::runtime_error( errorMessage );
        }
        return uniqueNodes_.at( n );
    }

    //! Get the unique weight factors for a specified order `n`.
    /*!
     * \param n The number of nodes or weight factors.
     * \return `uniqueWeights_[n]`, after reading the text file with the tabulated nodes if necessary.
     */
    IndependentVariableArray getUniqueWeights( const unsigned int n )
    {
        if ( uniqueWeights_.count( n ) == 0 )
        {
            std::string errorMessage = "Error in Gaussian quadrature, weights not available for n=" +
                    boost::lexical_cast< std::string >( n );
            throw std::runtime_error( errorMessage );
        }
        return uniqueWeights_.at( n );
    }
};

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_GAUSSIAN_QUADRATURE_H
