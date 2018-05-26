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

#ifndef TUDAT_EXTENDED_KALMAN_FILTER_H
#define TUDAT_EXTENDED_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/filter.h"

namespace tudat
{

namespace filters
{

//! Extended Kalman filter.
/*!
 *
 */
template< typename IndependentVariable = double, typename DependentVariable = double >
class ExtendedKalmanFilter: public FilterCore< IndependentVariable, DependentVariable >
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     */
    ExtendedKalmanFilter( FunctionPointer systemFunction,
                          FunctionPointer measurementFunction,
                          const Eigen::Matrix< DependentVariable, Eigen::Dynamic, Eigen::Dynamic >& systemUncertainty,
                          const Eigen::Matrix< DependentVariable, Eigen::Dynamic, Eigen::Dynamic >& measurementUncertainty,
                          FunctionPointer systemJacobian,
                          FunctionPointer measurementJacobian,
                          const bool isStateToBeIntegrated = false,
                          const IntegratorPointer integrator = NULL ) :
        FilterCore< IndependentVariable, DependentVariable >( isStateToBeIntegrated, integrator ),
        systemFunction_( systemFunction ), measurementFunction_( measurementFunction ),
        systemUncertainty_( systemUncertainty ), measurementUncertainty_( measurementUncertainty ),
        systemJacobian_( systemJacobian ), measurementJacobian_( measurementJacobian )
    { }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~ExtendedKalmanFilter( ){ }

    //!
    DependentVector updateFilter( const DependentVector& currentStateVector )
    {

    }

private:

    //!
    Eigen::Matrix< DependentVariable, Eigen::Dynamic, Eigen::Dynamic > systemUncertainty_;

    //!
    Eigen::Matrix< DependentVariable, Eigen::Dynamic, Eigen::Dynamic > measurementUncertainty_;

    //!
    FunctionPointer systemJacobian_;

    //!
    FunctionPointer measurementJacobian_;

};

} // namespace filters

} // namespace tudat

#endif // TUDAT_EXTENDED_KALMAN_FILTER_H
