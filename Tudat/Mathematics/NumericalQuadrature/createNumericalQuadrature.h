/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATENUMERICALQUADRATURE_H
#define TUDAT_CREATENUMERICALQUADRATURE_H

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/lexical_cast.hpp>

#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace numerical_quadrature
{

//! Enum to define available numerical quadratures.
enum AvailableQuadratures
{
    trapezoidal,
    gaussian
};

//! Class to define settings of numerical quadrature.
/*!
 *  Class to define settings of numerical quadrature.
 */
template< typename IndependentVariableType = double >
class QuadratureSettings
{
public:

    //! Constructor
    /*!
     *  Constructor for quadrature settings.
     *  \param quadratureType Type of numerical quadrature.
     */
    QuadratureSettings( const AvailableQuadratures quadratureType ) :
        quadratureType_( quadratureType )
    { }
    
    //! Virtual destructor.
    /*!
     *  Virtual destructor.
     */
    virtual ~QuadratureSettings( ) { }

    //! Type of numerical quadrature
    /*!
     *  Type of numerical quadrature, from enum of available quadratures.
     */
    AvailableQuadratures quadratureType_;

};


//! Class to define settings of gaussian quadrature.
/*!
 *  Class to define settings of gaussian quadrature.
 */
template< typename IndependentVariableType = double >
class GaussianQuadratureSettings: public QuadratureSettings< IndependentVariableType >
{
public:

    //! Default constructor.
    /*!
     *  Constructor for gaussian quadrature settings.
     *  \param initialIndependentVariable Starting independent variable of numerical quadrature.
     *  \param numberOfNodes Number of nodes at which the function to be integrated will be evaluated. Must be an integer value between 2 and 64.
     */
    GaussianQuadratureSettings(
            const IndependentVariableType initialIndependentVariable,
            const unsigned int numberOfNodes ) :
        QuadratureSettings< double >( gaussian ),
        initialIndependentVariable_( initialIndependentVariable ), numberOfNodes_( numberOfNodes )
    {
        if ( ( numberOfNodes_ < 2 ) || ( numberOfNodes_ > 64 ) )
        {
            throw std::runtime_error( "The number of nodes for the Gaussian quadrature is not an integer between 2 and 64." );
        }

    }

    //! Destructor.
    /*!
     *  Destructor.
     */
    ~GaussianQuadratureSettings( ) { }

    //! Starting independent variable of numerical quadrature.
    IndependentVariableType initialIndependentVariable_;

    //! Number of nodes at which the function to be integrated will be evaluated.
    unsigned int numberOfNodes_;

};


//! Class to define settings of trapezoid quadrature.
/*!
 *  Class to define settings of trapezoid quadrature.
 */
template< typename IndependentVariableType = double >
class TrapezoidQuadratureSettings: public QuadratureSettings< IndependentVariableType >
{
public:

    //! Default constructor.
    /*!
     *  Constructor for trapezoid quadrature settings.
     *  \param independentVariables Vector of independent variables at which the values of the derivative function will be given.
     */
    TrapezoidQuadratureSettings(
            const std::vector< IndependentVariableType >& independentVariables ) :
        QuadratureSettings< double >( trapezoidal ),
        independentVariables_( independentVariables ) { }

    //! Destructor.
    /*!
     *  Destructor.
     */
    ~TrapezoidQuadratureSettings( ) { }

    //! Vector of independent variables at which the values of the derivative function will be given.
    const std::vector< IndependentVariableType > independentVariables_;

};


//! Function to create a numerical quadrature.
/*!
 *  Function to create a numerical quadrature from given quadrature settings, and derivative function.
 *  \param derivativeFunction Function returning the derivative from current independent variable.
 *  \param quadratureSettings Settings for numerical quadrature.
 *  \param finalIndependentVariable Final independent variable of numerical quadrature.
 *  \return Numerical quadrature object.
 */
template< typename IndependentVariableType, typename DependentVariableType >
std::shared_ptr< numerical_quadrature::NumericalQuadrature< IndependentVariableType, DependentVariableType > > createQuadrature(
        std::function< DependentVariableType( const IndependentVariableType ) > derivativeFunction,
        std::shared_ptr< QuadratureSettings< IndependentVariableType > > quadratureSettings,
        IndependentVariableType finalIndependentVariable = TUDAT_NAN )
{
    // Declare eventual output.
    std::shared_ptr< NumericalQuadrature < IndependentVariableType, DependentVariableType > > quadrature;

    // Retrieve requested type of quadrature.
    switch( quadratureSettings->quadratureType_ )
    {
    case gaussian:
    {
        // Cast dynamic pointer on gaussian quadrature settings.
        std::shared_ptr< GaussianQuadratureSettings< double > > gaussianQuadratureSettings =
                std::dynamic_pointer_cast< GaussianQuadratureSettings < double > >( quadratureSettings );

        // Create gaussian quadrature.
        quadrature = std::make_shared< GaussianQuadrature < IndependentVariableType, DependentVariableType > >
                ( derivativeFunction, gaussianQuadratureSettings->initialIndependentVariable_, finalIndependentVariable,
                  gaussianQuadratureSettings->numberOfNodes_ ) ;
        break;
    }
    case trapezoidal:
    {
        // Cats dynamic pointer on trapezoid quadrature settings.
        std::shared_ptr< TrapezoidQuadratureSettings< double > > trapezoidQuadratureSettings =
                std::dynamic_pointer_cast< TrapezoidQuadratureSettings< double > >( quadratureSettings );

        // Retrieve independent variables at which the values of the derivative function must be given.
        std::vector< IndependentVariableType > independentVariables = trapezoidQuadratureSettings->independentVariables_;

        // Evaluate derivative function at each of the independent variable values.
        std::vector< DependentVariableType > dependentVariables;
        for ( int currentIndependentValue = 0 ; currentIndependentValue < independentVariables.size() ; currentIndependentValue++ )
        {
            dependentVariables.push_back( derivativeFunction( currentIndependentValue ) );
        }

        // Create Trapezoid quadrature.
        quadrature = std::make_shared< TrapezoidNumericalQuadrature < IndependentVariableType, DependentVariableType > >
                ( independentVariables, dependentVariables ) ;
        break;
    }
    default:
        throw std::runtime_error( "Error, quadrature " +  std::to_string( quadratureSettings->quadratureType_ ) + " not found." );
    }

    // Check that assignment of quadrature went well.
    if ( quadrature == nullptr )
    {
        throw std::runtime_error( "Error while creating quadrature. The resulting quadrature pointer is null." );
    }

    // Return numerical quadrature object.
    return quadrature;
}

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_CREATENUMERICALQUADRATURE_H
