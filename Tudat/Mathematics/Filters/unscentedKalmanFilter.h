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

#ifndef TUDAT_UNSCENTED_KALMAN_FILTER_H
#define TUDAT_UNSCENTED_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Unscented Kalman filter.
/*!
 *
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class UnscentedKalmanFilter: public FilterCore< IndependentVariableType, DependentVariableType >
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     */
    UnscentedKalmanFilter( const boost::shared_ptr< FunctionType > systemFunction,
                           const boost::shared_ptr< FunctionType > measurementFunction,
                           const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >& systemUncertainty,
                           const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >& measurementUncertainty,
                           const bool isStateToBeIntegrated = false,
                           const IntegratorPointer integrator = NULL ) :
        FilterCore< IndependentVariableType, DependentVariableType >( isStateToBeIntegrated, integrator ),
        systemFunction_( systemFunction ), measurementFunction_( measurementFunction ),
        systemUncertainty_( systemUncertainty ), measurementUncertainty_( measurementUncertainty )
    { }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~UnscentedKalmanFilter( ){ }

    //!
    DependentVector updateFilter( const DependentVector& currentStateVector )
    {

    }

private:

    //!
    Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic > systemUncertainty_;

    //!
    Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic > measurementUncertainty_;

};

} // namespace filters

} // namespace tudat

#endif // TUDAT_UNSCENTED_KALMAN_FILTER_H
