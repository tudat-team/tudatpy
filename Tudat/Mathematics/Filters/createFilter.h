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

#ifndef TUDAT_CREATE_FILTER_H
#define TUDAT_CREATE_FILTER_H

#include <map>
#include <limits>
#include <iostream>

#include <Eigen/Core>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

#include "Tudat/Mathematics/Filters/extendedKalmanFilter.h"
#include "Tudat/Mathematics/Filters/linearKalmanFilter.h"
#include "Tudat/Mathematics/Filters/unscentedKalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Enumeration of available filters.
enum AvailableFilteringTechniques
{
    linear_kalman_filter = 0,
    extended_kalman_filter = 1,
    unscented_kalman_filter = 2
};

//! Filter settings.
/*!
 *  Base class for the creation of a filter object.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class FilterSettings
{
public:

    //! Typedef of the state and measurement vectors.
    typedef Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > DependentVector;

    //! Typedef of the state and measurement matrices.
    typedef Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic > DependentMatrix;

    //! Typedef of the integrator settings.
    typedef numerical_integrators::IntegratorSettings< IndependentVariableType > IntegratorSettings;

    //! Constructor.
    /*!
     *  Constructor.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param isStateToBeIntegrated Boolean defining whether the system function needs to be integrated.
     *  \param integrator Pointer to integrator to be used to propagate state.
     */
    FilterSettings( const AvailableFilteringTechniques filteringTechnique,
                    const DependentMatrix& systemUncertainty,
                    const DependentMatrix& measurementUncertainty,
                    const IndependentVariableType initialTime,
                    const DependentVector& initialStateVector,
                    const DependentMatrix& initialCovarianceMatrix,
                    const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL ) :
        filteringTechnique_( filteringTechnique ), systemUncertainty_( systemUncertainty ),
        measurementUncertainty_( measurementUncertainty ), initialTime_( initialTime ),
        initialStateEstimate_( initialStateVector ), initialCovarianceEstimate_( initialCovarianceMatrix ),
        integratorSettings_( integratorSettings )
    { }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~FilterSettings( ){ }

    //! Enumeration denoting which filtering technique is to be used.
    const AvailableFilteringTechniques filteringTechnique_;

    //! Matrix representing the uncertainty in system modeling.
    const DependentMatrix systemUncertainty_;

    //! Matrix representing the uncertainty in measurement modeling.
    const DependentMatrix measurementUncertainty_;

    //! Scalar representing the initial time.
    const IndependentVariableType initialTime_;

    //! Vector representing the a-posteriori estimated state.
    /*!
     *  Vector representing the a-posteriori estimated state, i.e., the state after the prediction and
     *  update steps of the Kalman filter.
     */
    const DependentVector initialStateEstimate_;

    //! Matrix representing the a-posteriori estimated covariance.
    /*!
     *  Matrix representing the a-posteriori estimated covariance, i.e., the covariance after the prediction and
     *  update steps of the Kalman filter.
     */
    const DependentMatrix initialCovarianceEstimate_;

    //! Pointer to the integrator settings.
    /*!
     *  Pointer to the integrator settings, which are used to create the filter integrator.
     */
    const boost::shared_ptr< IntegratorSettings > integratorSettings_;

};

//! Unscented Kalman filter settings.
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class UnscentedKalmanFilterSettings : public FilterSettings< IndependentVariableType, DependentVariableType >
{
public:

    //! Inherit typedefs from base class.
    typedef typename FilterSettings< IndependentVariableType, DependentVariableType >::DependentVector DependentVector;
    typedef typename FilterSettings< IndependentVariableType, DependentVariableType >::DependentMatrix DependentMatrix;
    typedef typename FilterSettings< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;

    //! Default constructor.
    /*!
     *  Default constructor. This constructor takes the system and measurement functions as models for the simulation.
     *  These functions can be a function of time, state and (for system) control vector.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integratorSettings Pointer to integration settings defining the integrator to be used to propagate the state.
     *  \param constantValueReference Reference to be used for the values of the \f$ \alpha \f$ and \f$ \kappa \f$ parameters. This
     *      variable has to be part of the ConstantParameterReferences enumeration (custom parameters are supported).
     *  \param customConstantParameters Values of the constant parameters \f$ \alpha \f$ and \f$ \kappa \f$, in case the custom_parameters
     *      enumeration is used in the previous field.
     */
    UnscentedKalmanFilterSettings( const DependentMatrix& systemUncertainty,
                                   const DependentMatrix& measurementUncertainty,
                                   const IndependentVariableType initialTime,
                                   const DependentVector& initialStateVector,
                                   const DependentMatrix& initialCovarianceMatrix,
                                   const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL,
                                   const ConstantParameterReferences constantValueReference = reference_Wan_and_Van_der_Merwe,
                                   const std::pair< DependentVariableType, DependentVariableType > customConstantParameters =
            std::make_pair( static_cast< DependentVariableType >( TUDAT_NAN ),
                            static_cast< DependentVariableType >( TUDAT_NAN ) ) ) :
        FilterSettings< IndependentVariableType, DependentVariableType >( unscented_kalman_filter,
                                                                          systemUncertainty, measurementUncertainty,
                                                                          initialTime, initialStateVector,
                                                                          initialCovarianceMatrix, integratorSettings ),
        constantValueReference_( constantValueReference ), customConstantParameters_( customConstantParameters )
    { }

    //! Enumeration denoting the reference to use for the alpha and kappa paramters.
    const ConstantParameterReferences constantValueReference_;

    //! Custom value of the alpha and kappa paramters.
    const std::pair< DependentVariableType, DependentVariableType > customConstantParameters_;

};

//! Function to create a filter object with the use of filter settings.
template< typename IndependentVariableType = double, typename DependentVariableType = double >
boost::shared_ptr< filters::FilterBase< IndependentVariableType, DependentVariableType > >
createFilter( const boost::shared_ptr< FilterSettings< IndependentVariableType, DependentVariableType > > filterSettings,
              const boost::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >& systemFunction,
              const boost::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >& measurementFunction )
{
    // Create an empty filter object
    boost::shared_ptr< filters::FilterBase< IndependentVariableType, DependentVariableType > > createdFilter;

    // Check type of filter
    switch ( filterSettings->filteringTechnique_ )
    {
    case unscented_kalman_filter:
    {
        // Cast filter settings to unscented Kalman filter
        boost::shared_ptr< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >
                unscentedKalmanFilterSettings =
                boost::dynamic_pointer_cast< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >(
                    filterSettings );
        if ( unscentedKalmanFilterSettings == NULL )
        {
            throw std::runtime_error( "Error while creating unscented Kalman filter object. Type of filter "
                                      "settings (unscentedKalmanFilter) not compatible with "
                                      "selected filter (derived class of FilterSettings must be "
                                      "UnscentedKalmanFilterSettings for this type)." );
        }

        // Create filter
        createdFilter = boost::make_shared< UnscentedKalmanFilter< IndependentVariableType, DependentVariableType > >(
                    systemFunction, measurementFunction,
                    unscentedKalmanFilterSettings->systemUncertainty_, unscentedKalmanFilterSettings->measurementUncertainty_,
                    unscentedKalmanFilterSettings->initialTime_, unscentedKalmanFilterSettings->initialStateEstimate_,
                    unscentedKalmanFilterSettings->initialCovarianceEstimate_, unscentedKalmanFilterSettings->integratorSettings_,
                    unscentedKalmanFilterSettings->constantValueReference_, unscentedKalmanFilterSettings->customConstantParameters_ );
        break;
    }
    default:
        throw std::runtime_error( "Error while creating filter obejct. Only the creation of unscented Kalman filters is "
                                  "currently supported." );
    }

    // Check that filter was properly created
    if ( createdFilter == NULL )
    {
        throw std::runtime_error( "Error while creating filter. The resulting filter pointer is null." );
    }

    // Give output
    return createdFilter;
}

} // namespace filters

} // namespace tudat

#endif // TUDAT_CREATE_FILTER_H
