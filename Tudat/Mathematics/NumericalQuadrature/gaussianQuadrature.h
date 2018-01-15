/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

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
/*!
 *  Read Gaussian nodes from text file, file name is hard-coded into this function, read nodes are returned by reference
 *  \param gaussQuadratureNodes Gauss quadrature nodes that are read from the file. Map key denotes the order, map value the list
 *  of nodes for that order.
 */
template< typename IndependentVariableType >
void readGaussianQuadratureNodes(
        std::map< unsigned int, Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1> >& gaussQuadratureNodes )
{
    gaussQuadratureNodes =
            utilities::convertSTLVectorMapToEigenVectorMap< unsigned int, IndependentVariableType >(
                input_output::readStlVectorMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianNodes.txt" ) );
}

//! Read Gaussian weight factors from text file
/*!
*  Read Gaussian weight factors from text file, file name is hard-coded into this function, read nodes are returned by reference
*  \param gaussQuadratureWeights Gauss quadrature weights that are read from the file.  Map key denotes the order, map value the
 *  list of weights for that order.
*/
template< typename IndependentVariableType >
void readGaussianQuadratureWeights(
        std::map< unsigned int, Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1> >& gaussQuadratureWeights )
{
    gaussQuadratureWeights = utilities::convertSTLVectorMapToEigenVectorMap< unsigned int, IndependentVariableType >(
                input_output::readStlVectorMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianWeights.txt" ) );
}

//! Container object for Gauss quadrature nodes and weights (templated by data variable type, e.g. float, double, long double)
template< typename IndependentVariableType >
struct GaussQuadratureNodesAndWeights
{
    //! Typedef for vector of IndependentVariableType scalar type
    typedef Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1 > IndependentVariableArray;

    //! COnstructor, reads nodes and weights from file
    GaussQuadratureNodesAndWeights( )
    {
        readGaussianQuadratureNodes< IndependentVariableType >( uniqueNodes_ );
        readGaussianQuadratureWeights< IndependentVariableType >( uniqueWeights_ );
    }

    //! Get the unique nodes for a specified order `n`.
    /*!
     * \param numberOfNodes The number of nodes or weight factors.
     * \return `uniqueNodes_[n]`, after reading the text file with the tabulated nodes if necessary.
     */
    IndependentVariableArray getUniqueNodes( const unsigned int numberOfNodes )
    {
        if ( uniqueNodes_.count( numberOfNodes ) == 0 )
        {
            std::string errorMessage = "Error in Gaussian quadrature, nodes not available for n=" +
                    std::to_string( numberOfNodes );
            throw std::runtime_error( errorMessage );
        }
        return uniqueNodes_.at( numberOfNodes );
    }

    //! Get the unique weight factors for a specified order.
    /*!
     * Get the unique weight factors for a specified order.
     * \param order The number of nodes or weight factors.
     * \return `uniqueWeights_ at entry order`, after reading the text file with the tabulated nodes if necessary.
     */
    IndependentVariableArray getUniqueWeights( const unsigned int order )
    {
        if ( uniqueWeights_.count( order ) == 0 )
        {
            std::string errorMessage = "Error in Gaussian quadrature, weights not available for n=" +
                    std::to_string( order );
            throw std::runtime_error( errorMessage );
        }
        return uniqueWeights_.at( order );
    }

    //! Get all the nodes at given order from uniqueNodes_
    /*!
    * Get all the nodes at given order from uniqueNodes_
    * \param order The number of nodes or weight factors.
    * \return `uniqueWeights_ at entry order`, after reading the text file with the tabulated nodes if necessary.
    */
    IndependentVariableArray getNodes( const unsigned int order )
    {
        if ( nodes_.count( order ) == 0 )
        {
            IndependentVariableArray newNodes( order );

            // Include node 0.0 if order is odd
            unsigned int i = 0;
            if ( order % 2 == 1 )
            {
                newNodes.row( i++ ) = 0.0;
            }

            // Include ± nodes
            IndependentVariableArray uniqueNodes_ = getUniqueNodes( order );
            for ( int j = 0; j < uniqueNodes_.size( ); j++ )
            {
                newNodes.row( i++ ) = -uniqueNodes_[ j ];
                newNodes.row( i++ ) =  uniqueNodes_[ j ];
            }

            nodes_[ order ] = newNodes;
        }

        return nodes_.at( order );
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
            int j = 0;
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

};

//! Object containing nodes/weights for long double Gauss quadrature
static const boost::shared_ptr< GaussQuadratureNodesAndWeights< long double > > longDoubleGaussQuadratureNodesAndWeights =
        boost::make_shared< GaussQuadratureNodesAndWeights< long double > >( );

//! Object containing nodes/weights for double Gauss quadrature
static const boost::shared_ptr< GaussQuadratureNodesAndWeights< double > > doubleGaussQuadratureNodesAndWeights =
        boost::make_shared< GaussQuadratureNodesAndWeights< double > >( );

//! Object containing nodes/weights for float Gauss quadrature
static const boost::shared_ptr< GaussQuadratureNodesAndWeights< float > > floatGaussQuadratureNodesAndWeights =
        boost::make_shared< GaussQuadratureNodesAndWeights< float > >( );

//! Function to create Gauss quadrature node/weight container
/*!
 *  Function to create Gauss quadrature node/weight container, templated by independent variable type
 *  \return Gauss quadrature node/weight container
 */
template< typename IndependentVariableType >
boost::shared_ptr< GaussQuadratureNodesAndWeights< IndependentVariableType > >
getGaussQuadratureNodesAndWeights( )
{
    return boost::make_shared< GaussQuadratureNodesAndWeights< IndependentVariableType > >( );
}

//! Function to create Gauss quadrature node/weight container with long double precision.
/*!
 *  Function to create Gauss quadrature node/weight container with long double precision.
 *  \return Gauss quadrature node/weight container
 */
template< >
boost::shared_ptr< GaussQuadratureNodesAndWeights< long double > >
getGaussQuadratureNodesAndWeights( )
{
    return longDoubleGaussQuadratureNodesAndWeights;
}

//! Function to create Gauss quadrature node/weight container with double precision.
/*!
 *  Function to create Gauss quadrature node/weight container with double precision.
 *  \return Gauss quadrature node/weight container
 */
template< >
boost::shared_ptr< GaussQuadratureNodesAndWeights< double > >
getGaussQuadratureNodesAndWeights( )
{
    return doubleGaussQuadratureNodesAndWeights;
}

//! Function to create Gauss quadrature node/weight container with float precision.
/*!
 *  Function to create Gauss quadrature node/weight container with float precision.
 *  \return Gauss quadrature node/weight container
 */
template< >
boost::shared_ptr< GaussQuadratureNodesAndWeights< float > >
getGaussQuadratureNodesAndWeights( )
{
    return floatGaussQuadratureNodesAndWeights;
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


    typedef Eigen::Array< DependentVariableType, Eigen::Dynamic, 1 > DependentVariableArray;
    typedef Eigen::Array< IndependentVariableType, Eigen::Dynamic, 1 > IndependentVariableArray;

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
        gaussQuadratureNodesAndWeights_ = getGaussQuadratureNodesAndWeights< IndependentVariableType >( );
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


protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     * Function that is called to perform the numerical quadrature. Sets the result in the quadratureResult local
     * variable.
     */
    void performQuadrature( )
    {
        // Determine the values of the auxiliary independent variable (nodes)
        const IndependentVariableArray nodes = gaussQuadratureNodesAndWeights_->getNodes( numberOfNodes_ );

        // Determine the values of the weight factors
        const IndependentVariableArray weights = gaussQuadratureNodesAndWeights_->getWeights( numberOfNodes_ );

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

    boost::shared_ptr< GaussQuadratureNodesAndWeights< IndependentVariableType > > gaussQuadratureNodesAndWeights_;
};

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_GAUSSIAN_QUADRATURE_H
