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

#ifndef TUDAT_FILTER_H
#define TUDAT_FILTER_H

#include <limits>
#include <iostream>

#include <Eigen/Core>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

#include "Tudat/Mathematics/BasicMathematics/function.h"

namespace tudat
{

namespace filters
{

//! Filters class.
/*!
 *
 */
template< typename IndependentVariable = double, typename DependentVariable = double >
class FilterCore
{
public:

    //! Typedef of the function whose root we have to determine.
    typedef boost::shared_ptr< basic_mathematics::Function< IndependentVariable, DependentVariable > > FunctionPointer;

    //! Typedef of the state and measurement vectors.
    typedef Eigen::Matrix< DependentVariable, Eigen::Dynamic, 1 > DependentVector;

    //! Typedef of the integrator.
    typedef boost::shared_ptr< numerical_integrators::NumericalIntegrator< IndependentVariable, DependentVector > > IntegratorPointer;

    //! Constructor.
    /*!
     *  Constructor.
     */
    FilterCore( const bool isStateToBeIntegrated = false,
                const IntegratorPointer integrator = NULL ) :
        isStateToBeIntegrated_( isStateToBeIntegrated ), integrator_( integrator )
    { }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~FilterCore( ){ }

    //!
    virtual DependentVector updateFilter( const DependentVector& currentStateVector ) = 0;

    //! Update state history.
    void updateStateHistory( const Eigen::Matrix< DependentVariable, Eigen::Dynamic, Eigen::Dynamic >& stateHistory )
    {
        stateHistory_ = stateHistory;
    }

    //!
    DependentVector integrateState( const IndependentVariable intervalEnd,
                                    const IndependentVariable initialStepSize,
                                    const IndependentVariable finalTimeTolerance = std::numeric_limits< IndependentVariable >::epsilon( ) )
    {
        return integrator_->integrateTo( intervalEnd, initialStepSize, finalTimeTolerance );
    }

protected:

    //!
    /*!
     *
     */
    FunctionPointer systemFunction_;

    //!
    /*!
     *
     */
    FunctionPointer measurementFunction_;

    //!
    bool isStateToBeIntegrated_;

    //!
    IntegratorPointer integrator_;

    //!
    Eigen::Matrix< IndependentVariable, Eigen::Dynamic, 1 > independentVariable_;

    //!
    Eigen::Matrix< DependentVariable, Eigen::Dynamic, Eigen::Dynamic > stateHistory_;

};

//! Typedef for a filter with double data type.
typedef FilterCore< > Filter;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< Filter > FilterPointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_FILTER_H
