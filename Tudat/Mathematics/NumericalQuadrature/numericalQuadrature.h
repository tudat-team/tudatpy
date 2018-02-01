/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NUMERICAL_QUADRATURE_H
#define TUDAT_NUMERICAL_QUADRATURE_H

#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace numerical_quadrature
{

//! Base class for numerical quadrature.
/*!
 * Base class for the numerical quadrature methods included in Tudat, the dependent and independent variable
 * types are specified as template parameters.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class NumericalQuadrature
{
public:

    //! Destructor.
    virtual ~NumericalQuadrature( ) { }

    //! This function returns the value of the numerical integration performed by quadrature.
    /*!
     * This function returns the value of the numerical integration. It must be implemented in derived classes.
     * \return Integral of the data.
     */
    virtual DependentVariableType getQuadrature( ) = 0;


protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     * Function that is called to perform the numerical quadrature. It must be implemented in derived classes. Any
     * implementation of performQuadrature must be made consistent with that of getQuadrature.
     */
    virtual void performQuadrature( ) = 0;

private:


};

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_NUMERICAL_QUADRATURE_H
