/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/basics/timeType.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/identityElements.h"

namespace tudat
{

namespace interpolators
{

//! Enumeration for types of boundary interpolation methods.
/*!
 *  Enumeration for types of boundary interpolation methods, i.e., for when the independent variable requested for interpolation
 *  is outside the domain of the independent variables input in the interpolator constructor.
 */
enum BoundaryInterpolationType
{
    throw_exception_at_boundary = 0,
    use_boundary_value = 1,
    use_boundary_value_with_warning = 2,
    extrapolate_at_boundary = 3,
    extrapolate_at_boundary_with_warning = 4,
    use_default_value = 5,
    use_default_value_with_warning = 6,
    use_nan_value = 7,
    use_nan_value_with_warning = 8
};

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
     *  This function performs the interpolation. It must be implemented in derived classes.
     *  \param independentVariableValues Vector of values of independent variables at which
     *      the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable.
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

extern template class Interpolator< double, Eigen::VectorXd >;
extern template class Interpolator< double, Eigen::Vector6d >;
extern template class Interpolator< double, Eigen::MatrixXd >;

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_INTERPOLATOR_H
