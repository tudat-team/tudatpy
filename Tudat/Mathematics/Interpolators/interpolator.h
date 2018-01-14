/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_INTERPOLATOR_H
#define TUDAT_INTERPOLATOR_H

#include <vector>

namespace tudat
{
namespace interpolators
{

//! Base class for interpolator.
/*!
 * Base class for the interpolators included in Tudat, the dependent and independent variable
 * types are specified as user parameters, as are the number of dimensions. The number of
 * dimensions is chosen as a template parameter, so that a boost multi_array of this
 * dimension can be used as a member variable of derived classes.
 * \tparam IndependentVariableType Type of independent variable(s).
 * \tparam IndependentVariableType Type of dependent variable.
 * \tparam numberOfDimensions Number of independent directions for independent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class Interpolator
{
public:

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~Interpolator( ) { }

    //! Interpolate.
    /*!
     * This function performs the interpolation. It must be implemented in derived classes.
     * \param independentVariableValues Vector of values of independent variables at which
     *          the value of the dependent variable is to be determined.
     * \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType interpolate( const std::vector< IndependentVariableType >&
                                               independentVariableValues ) = 0;

    //! Function to return the number of independent variables of the interpolation.
    /*!
     *  Function to return the number of independent variables of the interpolation, i.e. size
     *  that the vector used as input for Interpolator::interpolate should be.
     *  \return Number of independent variables of the interpolation.
     */
    virtual int getNumberOfDimensions( ) = 0;


};

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_INTERPOLATOR_H
