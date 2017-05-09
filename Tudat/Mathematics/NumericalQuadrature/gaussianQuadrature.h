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

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{

namespace numerical_quadrature
{


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
    GaussianQuadrature( boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                        IndependentVariableType lowerLimit, IndependentVariableType upperLimit,
                        const unsigned int numberOfNodes )
    {
        reset( integrand, lowerLimit, upperLimit, numberOfNodes );
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
    void reset( boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                IndependentVariableType lowerLimit, IndependentVariableType upperLimit,
                const unsigned int numberOfNodes )
    {
        this->integrand            = integrand;
        this->lowerLimit           = lowerLimit;
        this->upperLimit           = upperLimit;
        this->numberOfNodes        = numberOfNodes;
        quadratureHasBeenPerformed = false;
    }


    //! Function to return computed value of the quadrature.
    /*!
     *  Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     *  \return Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     */
    DependentVariableType getQuadrature( )
    {
        if ( ! quadratureHasBeenPerformed ) {
            if ( integrand.empty() ) {
                throw std::runtime_error(
                            "The integrand for the Gaussian quadrature has not been set." );
            }

            if ( lowerLimit > upperLimit ) {
                throw std::runtime_error(
                            "The lower limit for the Gaussian quadrature is larger than the upper limit." );
            }

            if ( numberOfNodes < 2 || numberOfNodes > 64 ) {
                throw std::runtime_error(
                            "The number of nodes for the Gaussian quadrature must be between 2 and 64." );
            }

            performQuadrature();
            quadratureHasBeenPerformed = true;
        }

        return quadratureResult;
    }


    typedef Eigen::Array<   DependentVariableType, 1, Eigen::Dynamic >   DependentVariableArray;
    typedef Eigen::Array< IndependentVariableType, 1, Eigen::Dynamic > IndependentVariableArray;

    //! Get all the nodes (i.e. n nodes for nth order) from uniqueNodes
    IndependentVariableArray getNodes( const unsigned int n )
    {
        if ( nodes.count( n ) == 0 ) {
            IndependentVariableArray newNodes( n );

            // Include node 0.0 if n is odd
            unsigned int i = 0;
            if ( n % 2 == 1 ) {
                newNodes.col( i++ ) = 0.0;
            }

            // Include ± nodes
            IndependentVariableArray uniqueNodes = getUniqueNodes( n );
            for ( unsigned int j = 0; j < uniqueNodes.size(); j++ ) {
                newNodes.col( i++ ) = -uniqueNodes[ j ];
                newNodes.col( i++ ) =  uniqueNodes[ j ];
            }

            nodes[ n ] = newNodes;
        }

        return nodes.at( n );
    }

    //! Get all the weight factors (i.e. n weight factors for nth order) from uniqueWeights
    IndependentVariableArray getWeights( const unsigned int n )
    {
        if ( weights.count( n ) == 0 ) {
            IndependentVariableArray newWeights( n );

            IndependentVariableArray uniqueWeights = getUniqueWeights( n );

            // Include non-repeated weight factor if n is odd
            unsigned int i = 0;
            unsigned int j = 0;
            if ( n % 2 == 1 ) {
                newWeights.col( i++ ) = uniqueWeights[ j++ ];
            }

            // Include repeated weight factors
            for ( ; j < uniqueWeights.size(); j++ ) {
                newWeights.col( i++ ) = uniqueWeights[ j ];
                newWeights.col( i++ ) = uniqueWeights[ j ];
            }

            weights[ n ] = newWeights;
        }

        return weights.at( n );
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
        const IndependentVariableArray nodes = getNodes( numberOfNodes );

        // Determine the values of the weight factors
        const IndependentVariableArray weights = getWeights( numberOfNodes );

        // Change of variable -> from range [-1, 1] to range [lowerLimit, upperLimit]
        const IndependentVariableArray independentVariables =
                0.5 * ( ( upperLimit - lowerLimit ) * nodes + upperLimit + lowerLimit );

        // Determine the value of the dependent variable
        DependentVariableArray weighedIntegrands( numberOfNodes );
        for ( unsigned int i = 0; i < numberOfNodes; i++ ) {
            weighedIntegrands( i ) = weights( i ) * integrand( independentVariables( i ) );
        }

        quadratureResult = 0.5 * ( upperLimit - lowerLimit ) * weighedIntegrands.sum();
    }


private:

    //! Map containing the nodes read from the text file (currently up to `n = 64`).
    //! The following relation holds: `size( uniqueNodes[n] ) = floor( n / 2 )`
    //! For the actual nodes, the following must hold: `size( nodes[n] ) = n`
    //! The actual nodes are generated from `uniqueNodes` by `getNodes()`
    std::map< unsigned int, IndependentVariableArray > uniqueNodes;
    std::map< unsigned int, IndependentVariableArray > nodes;

    //! Map containing the weight factors read from the text file (currently up to `n = 64`).
    //! The following relation holds: `size( uniqueWeights[n] ) = ceil( n / 2 )`
    //! For the actual weight factors, the following must hold: `size( uniqueWeights[n] ) = n`
    //! The actual weight factors are generated from `uniqueWeights` by `getWeights()`
    std::map< unsigned int, IndependentVariableArray > uniqueWeights;
    std::map< unsigned int, IndependentVariableArray > weights;

    //! Function returning the integrand.
    boost::function< DependentVariableType( IndependentVariableType ) > integrand;

    //! Lower limit for the integral.
    IndependentVariableType lowerLimit = 0;

    //! Upper limit for the integral.
    IndependentVariableType upperLimit = 0;

    //! Number of nodes.
    unsigned int numberOfNodes = 0;

    //! Whether quadratureResult has been set for the current integrand, lowerLimit, upperLimit and numberOfNodes.
    bool quadratureHasBeenPerformed = false;

    //! Computed value of the quadrature, as computed by last call to performQuadrature.
    DependentVariableType quadratureResult;

    //! Get the unique nodes for a specified order `n`.
    /*!
     * \param n The number of nodes or weight factors.
     * \return `uniqueNodes[n]`, after reading the text file with the tabulated nodes if necessary.
     */
    IndependentVariableArray getUniqueNodes( const unsigned int n )
    {
        if ( uniqueNodes.count( n ) == 0 )
        {
            readNodes();  // this reads all the nodes from n=2 to n=64
        }
        return uniqueNodes.at( n );
    }

    //! Get the unique weight factors for a specified order `n`.
    /*!
     * \param n The number of nodes or weight factors.
     * \return `uniqueWeights[n]`, after reading the text file with the tabulated nodes if necessary.
     */
    IndependentVariableArray getUniqueWeights( const unsigned int n )
    {
        if ( uniqueWeights.count( n ) == 0 )
        {
            readWeights();  // this reads all the weights from n=2 to n=64
        }
        return uniqueWeights.at( n );
    }


    //! Transform from map of std::vector (output of text file reader) to map of Eigen::Array
    std::map< unsigned int, IndependentVariableArray > vector2eigen(
            std::map< unsigned int, std::vector< IndependentVariableType > > map )
    {
        std::map< unsigned int, IndependentVariableArray > eigenMap;
        for ( auto ent: map )
        {
            IndependentVariableArray array( ent.second.size() );
            for ( unsigned int i = 0; i < array.cols(); i++ )
            {
                array.col( i ) = ent.second.at( i );
            }
            eigenMap[ ent.first ] = array;
        }
        return eigenMap;
    }

    //! Read Gaussian nodes from text file
    void readNodes()
    {
        auto nodes = input_output::readMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianNodes.txt" );
        uniqueNodes = vector2eigen( nodes );
    }

    //! Read Gaussian weight factors from text file
    void readWeights()
    {
        auto weights = input_output::readMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianWeights.txt" );
        uniqueWeights = vector2eigen( weights );
    }


};

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_GAUSSIAN_QUADRATURE_H
