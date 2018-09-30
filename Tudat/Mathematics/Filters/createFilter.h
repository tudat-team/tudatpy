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
     *  \param filteringTechnique Enumeration denoting the type of filtering technique.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param filteringStepSize Scalar representing the value of the constant filtering time step.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integratorSettings Settings for the numerical integrator to be used to propagate state.
     */
    FilterSettings( const AvailableFilteringTechniques filteringTechnique,
                    const DependentMatrix& systemUncertainty,
                    const DependentMatrix& measurementUncertainty,
                    const IndependentVariableType filteringStepSize,
                    const IndependentVariableType initialTime,
                    const DependentVector& initialStateVector,
                    const DependentMatrix& initialCovarianceMatrix,
                    const std::shared_ptr< IntegratorSettings > integratorSettings = nullptr ) :
        filteringTechnique_( filteringTechnique ), systemUncertainty_( systemUncertainty ),
        measurementUncertainty_( measurementUncertainty ), filteringStepSize_( filteringStepSize ), initialTime_( initialTime ),
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

    //! Scalar representing step-size for filtering process.
    const IndependentVariableType filteringStepSize_;

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
    const std::shared_ptr< IntegratorSettings > integratorSettings_;

};

//! Extended Kalman filter settings.
/*!
 *  Extended Kalman filter settings.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class ExtendedKalmanFilterSettings : public FilterSettings< IndependentVariableType, DependentVariableType >
{
public:

    //! Inherit typedefs from base class.
    typedef typename FilterSettings< IndependentVariableType, DependentVariableType >::DependentVector DependentVector;
    typedef typename FilterSettings< IndependentVariableType, DependentVariableType >::DependentMatrix DependentMatrix;
    typedef typename FilterSettings< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;

    //! Default constructor.
    /*!
     *  Default constructor. This constructor takes state and measurement functions and their respective
     *  Jacobian functions as inputs. These functions can be a function of time, state and (for state) control vector.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param filteringStepSize Scalar representing the value of the constant filtering time step.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integratorSettings Pointer to integration settings defining the integrator to be used to propagate the state.
     */
    ExtendedKalmanFilterSettings( const DependentMatrix& systemUncertainty,
                                  const DependentMatrix& measurementUncertainty,
                                  const IndependentVariableType filteringStepSize,
                                  const IndependentVariableType initialTime,
                                  const DependentVector& initialStateVector,
                                  const DependentMatrix& initialCovarianceMatrix,
                                  const std::shared_ptr< IntegratorSettings > integratorSettings = nullptr ) :
        FilterSettings< IndependentVariableType, DependentVariableType >( extended_kalman_filter,
                                                                          systemUncertainty, measurementUncertainty,
                                                                          filteringStepSize, initialTime, initialStateVector,
                                                                          initialCovarianceMatrix, integratorSettings )
    { }

};

//! Unscented Kalman filter settings.
/*!
 *  Unscented Kalman filter settings.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
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
     *  These functions can be a function of time and state vector.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param filteringStepSize Scalar representing the value of the constant filtering time step.
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
                                   const IndependentVariableType filteringStepSize,
                                   const IndependentVariableType initialTime,
                                   const DependentVector& initialStateVector,
                                   const DependentMatrix& initialCovarianceMatrix,
                                   const std::shared_ptr< IntegratorSettings > integratorSettings = nullptr,
                                   const ConstantParameterReferences constantValueReference = reference_Wan_and_Van_der_Merwe,
                                   const std::pair< DependentVariableType, DependentVariableType > customConstantParameters =
            std::make_pair( static_cast< DependentVariableType >( TUDAT_NAN ),
                            static_cast< DependentVariableType >( TUDAT_NAN ) ) ) :
        FilterSettings< IndependentVariableType, DependentVariableType >( unscented_kalman_filter,
                                                                          systemUncertainty, measurementUncertainty,
                                                                          filteringStepSize, initialTime, initialStateVector,
                                                                          initialCovarianceMatrix, integratorSettings ),
        constantValueReference_( constantValueReference ), customConstantParameters_( customConstantParameters )
    { }

    //! Enumeration denoting the reference to use for the alpha and kappa paramters.
    const ConstantParameterReferences constantValueReference_;

    //! Custom value of the alpha and kappa paramters.
    const std::pair< DependentVariableType, DependentVariableType > customConstantParameters_;

};

//! Function to create a filter object with the use of filter settings.
/*!
 *  Function to create a filter object with the use of filter settings.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 *  \param filterSettings Settings for the creation of the filter object.
 *  \param systemFunction Function returning the state as a function of time and state vector. Can be a differential
 *      equation if the integratorSettings is set (i.e., if it is not a nullptr).
 *  \param measurementFunction Function returning the measurement as a function of time and state.
 *  \param stateJacobianFunction Function returning the Jacobian of the system w.r.t. the state. The input values can
 *      be time and state vector.
 *  \param stateNoiseJacobianFunction Function returning the Jacobian of the system function w.r.t. the system noise. The input
 *      values can be time and state vector.
 *  \param measurementJacobianFunction Function returning the Jacobian of the measurement function w.r.t. the state. The input
 *      values can be time and state.
 *  \param measurementNoiseJacobianFunction Function returning the Jacobian of the measurement function w.r.t. the measurement
 *      noise. The input values can be time and state.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
std::shared_ptr< filters::FilterBase< IndependentVariableType, DependentVariableType > >
createFilter( const std::shared_ptr< FilterSettings< IndependentVariableType, DependentVariableType > > filterSettings,
              const std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >& systemFunction,
              const std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >& measurementFunction,
              const std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >&
              stateJacobianFunction = std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >( ),
              const std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >&
              stateNoiseJacobianFunction = std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >( ),
              const std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >&
              measurementJacobianFunction = std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >( ),
              const std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >&
              measurementNoiseJacobianFunction = std::function< Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic >(
                  const IndependentVariableType, const Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >& ) >( ) )
{
    // Create an empty filter object
    std::shared_ptr< filters::FilterBase< IndependentVariableType, DependentVariableType > > createdFilter;

    // Check type of filter
    switch ( filterSettings->filteringTechnique_ )
    {
    case extended_kalman_filter:
    {
        // Cast filter settings to extended Kalman filter
        std::shared_ptr< ExtendedKalmanFilterSettings< IndependentVariableType, DependentVariableType > > extendedKalmanFilterSettings =
                std::dynamic_pointer_cast< ExtendedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >(
                    filterSettings );
        if ( extendedKalmanFilterSettings == nullptr )
        {
            throw std::runtime_error( "Error while creating extended Kalman filter object. Type of filter settings "
                                      "(ExtendedKalmanFilter) not compatible with selected filter (derived class of FilterSettings "
                                      "must be ExtendedKalmanFilterSettings for this type)." );
        }

        // Check that optional inputs are present
        if ( ( stateJacobianFunction == nullptr ) || ( stateNoiseJacobianFunction == nullptr ) ||
             ( measurementJacobianFunction == nullptr ) || ( measurementNoiseJacobianFunction == nullptr ) )
        {
            throw std::runtime_error( "Error while creating extended Kalman filter object. An ExtendedKalmanFilter object "
                                      "requires the input of the four Jacobian functions for state and measurement (including noise)." );
        }

        // Create filter
        createdFilter = std::make_shared< ExtendedKalmanFilter< IndependentVariableType, DependentVariableType > >(
                    systemFunction, measurementFunction, stateJacobianFunction, stateNoiseJacobianFunction,
                    measurementJacobianFunction, measurementNoiseJacobianFunction,
                    extendedKalmanFilterSettings->systemUncertainty_, extendedKalmanFilterSettings->measurementUncertainty_,
                    extendedKalmanFilterSettings->filteringStepSize_, extendedKalmanFilterSettings->initialTime_,
                    extendedKalmanFilterSettings->initialStateEstimate_, extendedKalmanFilterSettings->initialCovarianceEstimate_,
                    extendedKalmanFilterSettings->integratorSettings_ );
        break;
    }
    case unscented_kalman_filter:
    {
        // Cast filter settings to unscented Kalman filter
        std::shared_ptr< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > > unscentedKalmanFilterSettings =
                std::dynamic_pointer_cast< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >(
                    filterSettings );
        if ( unscentedKalmanFilterSettings == nullptr )
        {
            throw std::runtime_error( "Error while creating unscented Kalman filter object. Type of filter settings "
                                      "(UnscentedKalmanFilter) not compatible with selected filter (derived class of FilterSettings "
                                      "must be UnscentedKalmanFilterSettings for this type)." );
        }

        // Create filter
        createdFilter = std::make_shared< UnscentedKalmanFilter< IndependentVariableType, DependentVariableType > >(
                    systemFunction, measurementFunction,
                    unscentedKalmanFilterSettings->systemUncertainty_, unscentedKalmanFilterSettings->measurementUncertainty_,
                    unscentedKalmanFilterSettings->filteringStepSize_, unscentedKalmanFilterSettings->initialTime_,
                    unscentedKalmanFilterSettings->initialStateEstimate_, unscentedKalmanFilterSettings->initialCovarianceEstimate_,
                    unscentedKalmanFilterSettings->integratorSettings_, unscentedKalmanFilterSettings->constantValueReference_,
                    unscentedKalmanFilterSettings->customConstantParameters_ );
        break;
    }
    default:
        throw std::runtime_error( "Error while creating filter obejct. The creation of linear filters is not yet supported." );
    }

    // Check that filter was properly created
    if ( createdFilter == nullptr )
    {
        throw std::runtime_error( "Error while creating filter. The resulting filter pointer is null." );
    }

    // Give output
    return createdFilter;
}

} // namespace filters

} // namespace tudat

#endif // TUDAT_CREATE_FILTER_H
