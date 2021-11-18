/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEINTERPOLATOR_H
#define TUDAT_CREATEINTERPOLATOR_H

#include <iostream>

#include <boost/make_shared.hpp>
#include <memory>

#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/math/interpolators/hermiteCubicSplineInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/math/interpolators/piecewiseConstantInterpolator.h"

#include "tudat/math/interpolators/multiLinearInterpolator.h"

#include "tudat/io/mapTextFileReader.h"

namespace tudat
{

namespace interpolators
{


//! Base class for providing settings for creating an interpolator.
/*!
 *  Base class for providing settings for creating an interpolator using the createInterpolator
 *  function. This base class is not-functional, i.e. a derived class needs to be used.
 */
class InterpolatorSettings
{  
public:

    //! Default constructor.
    /*!
     *  Default constructor. Constructor taking a vector of boundary handling methods. The vector length needs
     *  to be equal to the number of dimensions.
     *  \param interpolatorType Selected type of interpolator.
     *  \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     *  \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     *      time step is a double if false.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     */
    InterpolatorSettings( const InterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                          const bool useLongDoubleTimeStep = false,
                          const std::vector< BoundaryInterpolationType >& boundaryHandling = std::vector< BoundaryInterpolationType >( ) ) :
        interpolatorType_( interpolatorType ), selectedLookupScheme_( selectedLookupScheme ),
        useLongDoubleTimeStep_( useLongDoubleTimeStep ), boundaryHandling_( boundaryHandling )
    {
        // Check that if interpolator type matches with number of dimensions
        std::vector< bool > isMethodOneDimensional = std::vector< bool >( 6, true );
        isMethodOneDimensional.at( static_cast< unsigned int >( multi_linear_interpolator ) ) = false;
        if ( boundaryHandling_.size( ) > 1 && isMethodOneDimensional.at( static_cast< unsigned int >( interpolatorType_ ) ) )
        {
            throw std::runtime_error( "Error while creating interpolator settings. Number of dimensions is greater "
                                      "than 1, but a one-dimensional interpolator has been selected." );
        }

        // Assign default value for boundary handling for one dimension
        if ( boundaryHandling_.empty( ) && isMethodOneDimensional.at( static_cast< unsigned int >( interpolatorType_ ) ) )
        {
            boundaryHandling_ = std::vector< BoundaryInterpolationType >( 1, extrapolate_at_boundary );
        }
    }

    //! Constructor.
    /*!
     *  Constructor. Constructor taking a single boundary handling method.
     *  \param interpolatorType Selected type of interpolator.
     *  \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     *  \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     *      time step is a double if false.
     *  \param boundaryHandling Boundary handling method, in case independent variable is outside the
     *      specified range.
     */
    InterpolatorSettings( const InterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme,
                          const bool useLongDoubleTimeStep,
                          const BoundaryInterpolationType boundaryHandling ) :
        InterpolatorSettings( interpolatorType, selectedLookupScheme, useLongDoubleTimeStep,
                              std::vector< BoundaryInterpolationType >( 1, boundaryHandling ) )
    { }

    //! Virtual destructor.
    virtual ~InterpolatorSettings( ) { }

    //! Function to get the selected type of interpolator.
    /*!
     * Function to get the selected type of interpolator.
     * \return Selected type of interpolator.
     */
    InterpolatorTypes getInterpolatorType( )
    {
        return interpolatorType_;
    }

    //! Function to get the selected type of lookup scheme for independent variables.
    /*!
     * Function to get the selected type of lookup scheme for independent variables.
     * \return Selected type of lookup scheme for independent variables.
     */
    AvailableLookupScheme getSelectedLookupScheme( )
    {
        return selectedLookupScheme_;
    }

    //! Function to reset the use of long double type for time step.
    /*!
     * Function to reset the use of long double type for time step.
     * \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double.
     */
    void resetUseLongDoubleTimeStep( const bool useLongDoubleTimeStep )
    {
        useLongDoubleTimeStep_ = useLongDoubleTimeStep;
    }

    //! Function to get a boolean denoting whether time step is to be a long double.
    /*!
     * Function to get a boolean denoting whether time step is to be a long double.
     * \return Boolean denoting whether time step is to be a long double.
     */
    bool getUseLongDoubleTimeStep( )
    {
        return useLongDoubleTimeStep_;
    }

    //! Function to retrieve boundary handling method.
    /*!
     * Function to retrieve boundary handling method.
     * \return Boundary handling method.
     */
    std::vector< BoundaryInterpolationType > getBoundaryHandling( )
    {
        return boundaryHandling_;
    }

protected:

    //! Selected type of interpolator.
    InterpolatorTypes interpolatorType_;

    //! Selected type of lookup scheme for independent variables.
    AvailableLookupScheme selectedLookupScheme_;

    //! Boolean denoting whether time step is to be a long double.
    bool useLongDoubleTimeStep_;

    //! Boundary handling method.
    std::vector< BoundaryInterpolationType > boundaryHandling_;

};

//! Class for providing settings to creating a Lagrange interpolator.
class LagrangeInterpolatorSettings : public InterpolatorSettings
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param interpolatorOrder Order of the Lagrange interpolator that is to be created.
     * \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     * time step is a double if false.
     * \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     * \param lagrangeBoundaryHandling Lagrange boundary handling method, in case the independent variable is outside the
     * specified range.
     * \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     * specified range.
     */
    LagrangeInterpolatorSettings(
            const int interpolatorOrder,
            const bool useLongDoubleTimeStep = 0,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const LagrangeInterpolatorBoundaryHandling lagrangeBoundaryHandling = lagrange_cubic_spline_boundary_interpolation,
            const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary ):
        InterpolatorSettings( lagrange_interpolator, selectedLookupScheme, useLongDoubleTimeStep, boundaryHandling ),
        interpolatorOrder_( interpolatorOrder ),
        lagrangeBoundaryHandling_( lagrangeBoundaryHandling )
    { }

    //! Destructor
    ~LagrangeInterpolatorSettings( ){ }

    //! Function to get the order of the Lagrange interpolator that is to be created.
    /*!
     * Function to get the order of the Lagrange interpolator that is to be created.
     * \return Order of the Lagrange interpolator that is to be created.
     */
    int getInterpolatorOrder( )
    {
        return interpolatorOrder_;
    }

    //! Function to retrieve Lagrange boundary handling method.
    /*!
     * Function to retrieve Lagrange boundary handling method.
     * \return Lagrange boundary handling method.
     */
    LagrangeInterpolatorBoundaryHandling getLagrangeBoundaryHandling( )
    {
        return lagrangeBoundaryHandling_;
    }

protected:

    //! Order of the Lagrange interpolator that is to be created.
    int interpolatorOrder_;

    //! Lagrange boundary handling method.
    LagrangeInterpolatorBoundaryHandling lagrangeBoundaryHandling_;

};

inline std::shared_ptr< InterpolatorSettings > linearInterpolation(
        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
        const BoundaryInterpolationType boundaryInterpolation = extrapolate_at_boundary_with_warning )
{
    return std::make_shared< InterpolatorSettings >(
                linear_interpolator, selectedLookupScheme, false, boundaryInterpolation );
}

inline std::shared_ptr< InterpolatorSettings > cubicSplineInterpolation(
        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
        const BoundaryInterpolationType boundaryInterpolation = extrapolate_at_boundary_with_warning )
{
    return std::make_shared< InterpolatorSettings >(
                cubic_spline_interpolator, selectedLookupScheme, false, boundaryInterpolation );
}


inline std::shared_ptr< InterpolatorSettings > piecewiseConstantInterpolation(
        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
        const BoundaryInterpolationType boundaryInterpolation = extrapolate_at_boundary_with_warning )
{
    return std::make_shared< InterpolatorSettings >(
                piecewise_constant_interpolator, selectedLookupScheme, false, boundaryInterpolation );
}


inline std::shared_ptr< InterpolatorSettings > hermiteInterpolation(
        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
        const BoundaryInterpolationType boundaryInterpolation = extrapolate_at_boundary_with_warning )
{
    return std::make_shared< InterpolatorSettings >(
                hermite_spline_interpolator, selectedLookupScheme, false, boundaryInterpolation );
}

inline std::shared_ptr< InterpolatorSettings > lagrangeInterpolation(
        const int order,
        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
        const BoundaryInterpolationType boundaryInterpolation = extrapolate_at_boundary_with_warning,
        const LagrangeInterpolatorBoundaryHandling lagrangeBoundaryHandling = lagrange_cubic_spline_boundary_interpolation )
{
    return std::make_shared< LagrangeInterpolatorSettings >(
                order, false, selectedLookupScheme, lagrangeBoundaryHandling, boundaryInterpolation );
}

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief DataMapSettings
 * The class can be used to provided a map directly, or as a base class for loading the map using different procedures.
 */
template< typename IndependentType, typename DependentType >
class DataMapSettings
{
public:

    //! Empty constructor.
    /*!
     * Empty constructor.
     */
    DataMapSettings( ) { }

    //! Constructor with a data map.
    /*!
     * Constructor to be used when the data map is provided directly.
     * \param dataMap The data map containing values for the independent and dependent variables.
     */
    DataMapSettings( const std::map< IndependentType, DependentType >& dataMap ) : dataMap_( dataMap ) { }

    //! Virtual destructor
    virtual ~DataMapSettings( ) { }

    //! Get the data map associated to the current setting object.
    /*!
     * @copybrief getDataMap
     * \return The data map associated to the current setting object.
     */
    virtual std::map< IndependentType, DependentType > getDataMap( ) const
    {
        return dataMap_;
    }

protected:

    //! The data map directly provided by the user in the constructor.
    const std::map< IndependentType, DependentType > dataMap_;

};

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief IndependentDependentDataMapSettings
 * The data map will be created by combining values of independent and dependent variables.
 */
template< typename IndependentType, typename DependentType >
class IndependentDependentDataMapSettings : public DataMapSettings< IndependentType, DependentType >
{
public:

    //! Consturctor.
    /*!
     * Constructor.
     * \param independentVariableValues Vector containing the values of the indepedent variable.
     * \param dependentVariableValues Vector containing the values of the depedent variable.
     */
    IndependentDependentDataMapSettings( const std::vector< IndependentType >& independentVariableValues,
                                         const std::vector< DependentType >& dependentVariableValues ) :
        DataMapSettings< IndependentType, DependentType >( ),
        independentVariableValues_( independentVariableValues ),
        dependentVariableValues_( dependentVariableValues ) { }

    //! Virtual destructor
    virtual ~IndependentDependentDataMapSettings( ) { }

    //! Vector containing the values of the indepedent variable.
    std::vector< IndependentType > independentVariableValues_;

    //! Vector containing the values of the depedent variable.
    std::vector< DependentType > dependentVariableValues_;

    //! Get the data map associated to the current setting object.
    /*!
     * @copybrief getDataMap
     * The map is created from the provided vectors of independent and depedent variables values.
     * \return The data map associated to the current setting object.
     */
    virtual std::map< IndependentType, DependentType > getDataMap( ) const
    {
        if ( independentVariableValues_.size( ) != dependentVariableValues_.size( ) )
        {
            std::cerr << "Could not get data map because the size of the independent and dependent variables values "
                         "is inconsistent." << std::endl;
            throw;
        }
        std::map< IndependentType, DependentType > dataMap;
        for ( unsigned int i = 0; i < independentVariableValues_.size( ); ++i )
        {
            dataMap[ independentVariableValues_.at( i ) ] = dependentVariableValues_.at( i );
        }
        return dataMap;
    }

};

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief FromFileDataMapSettings
 * The data map will be read from a file.
 */
template< typename EigenVectorType >
class FromFileDataMapSettings : public DataMapSettings< typename EigenVectorType::Scalar, EigenVectorType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param relativeFilePath Relative path to the file from which the map is to be loaded.
     */
    FromFileDataMapSettings( const std::string& relativeFilePath ) :
        DataMapSettings< typename EigenVectorType::Scalar, EigenVectorType >( ),
        relativeFilePath_( relativeFilePath ) { }

    //! Virtual destructor
    virtual ~FromFileDataMapSettings( ) { }

    //! Relative path to the file from which the map is to be loaded.
    std::string relativeFilePath_;

    //! Get the data map associated to the current setting object.
    /*!
     * @copybrief getDataMap
     * The map is loaded from the specified relativeFilePath_.
     * \return The data map associated to the current setting object.
     */
    virtual std::map< typename EigenVectorType::Scalar, EigenVectorType > getDataMap( ) const
    {
        return input_output::readEigenVectorMapFromFile< EigenVectorType >( relativeFilePath_ );
    }

};

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief HermiteDataSettings
 * The data map will be created directly, in addition to the first derivatives.
 */
template< typename IndependentType, typename DependentType >
class HermiteDataSettings : public DataMapSettings< IndependentType, DependentType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param dataToInterpolate The data map containing values for the independent and dependent variables.
     * \param firstDerivativeOfDependentVariables Vector containing the first derivatives of the depedent variables.
     */
    HermiteDataSettings( const std::map< IndependentType, DependentType >& dataToInterpolate,
                         const std::vector< DependentType >& firstDerivativeOfDependentVariables ) :
        DataMapSettings< IndependentType, DependentType >( dataToInterpolate ),
        firstDerivativeOfDependentVariables_( firstDerivativeOfDependentVariables ) { }

    //! Virtual destructor
    virtual ~HermiteDataSettings( ) { }

    //! Vector containing the first derivatives of the depedent variables.
    std::vector< DependentType > firstDerivativeOfDependentVariables_;

};

//! Class containing (the settings to create) the data needed for the interpolation and the settings to create the
//! interpolator.
/*!
 * @copybrief DataInterpolationSettings
 */
template< typename IndependentType, typename DependentType >
class DataInterpolationSettings
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param dataSettings Object containing (the settings to create) the data needed for the interpolation.
     * \param interpolatorSettings Object containing the settings to create the interpolator to be used.
     */
    DataInterpolationSettings(
            const std::shared_ptr< DataMapSettings< IndependentType, DependentType > >& dataSettings,
            const std::shared_ptr< InterpolatorSettings >& interpolatorSettings ) :
        dataSettings_( dataSettings ), interpolatorSettings_( interpolatorSettings ) { }

    //! Virtual destructor
    virtual ~DataInterpolationSettings( ){ }

    //! Object containing (the settings to create) the data needed for the interpolation.
    std::shared_ptr< DataMapSettings< IndependentType, DependentType > > dataSettings_;

    //! Object containing the settings to create the interpolator to be used.
    std::shared_ptr< InterpolatorSettings > interpolatorSettings_;
};

//! Function to create a one-dimensional interpolator
/*!
 *  Function to create a one-dimensional interpolator from the data that is to be interpolated,
 *  as well as the settings that are to be used to create the interpolator.
 *  \param dataToInterpolate Map providing data that is to be interpolated (key = independent
 *      variables, value = dependent variables).
 *  \param interpolatorSettings Settings that are to be used to create interpolator.
 *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
        of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
 *  \param firstDerivativeOfDependentVariables First derivative of dependent variables w.r.t. independent variable at
 *      independent variables values in values of dataToInterpolate. By default, this vector is empty, it only needs to
 *      be supplied if the selected interpolator requires this data (e.g. Hermite spline).
 *  \return Interpolator created from dataToInterpolate using interpolatorSettings.
 */
template< typename IndependentVariableType, typename DependentVariableType >
std::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
createOneDimensionalInterpolator(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const std::shared_ptr< InterpolatorSettings > interpolatorSettings,
        const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
        std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                        IdentityElement::getAdditionIdentity< DependentVariableType >( ) ),
        const std::vector< DependentVariableType > firstDerivativeOfDependentVariables =
        std::vector< DependentVariableType >( ) )
{
    std::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
            createdInterpolator;

    // Check that boundary handling is one-dimensional
    if ( interpolatorSettings->getBoundaryHandling( ).size( ) != 1 )
    {
        throw std::runtime_error( "Error while creating interpolator of type: " +
                                  std::to_string( interpolatorSettings->getInterpolatorType( ) ) +
                                  ". The interpolator is one-dimensional, but more than one boundary "
                                  "handling methods have been defined." );
    }

    // Check type of interpolator.
    switch( interpolatorSettings->getInterpolatorType( ) )
    {
    case linear_interpolator:
        createdInterpolator = std::make_shared< LinearInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryHandling( ).at( 0 ), defaultExtrapolationValue );
        break;
    case cubic_spline_interpolator:
    {
        if( !interpolatorSettings->getUseLongDoubleTimeStep( ) )
        {
            createdInterpolator = std::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType > >(
                        dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                        interpolatorSettings->getBoundaryHandling( ).at( 0 ) );
        }
        else
        {
            createdInterpolator = std::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType, long double > >(
                        dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                        interpolatorSettings->getBoundaryHandling( ).at( 0 ), defaultExtrapolationValue );
        }
        break;
    }
    case lagrange_interpolator:
    {
        // Check consistency of input
        std::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                std::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );
        if( lagrangeInterpolatorSettings != nullptr )
        {
            // Create Lagrange interpolator with requested time step type
            if( !lagrangeInterpolatorSettings->getUseLongDoubleTimeStep( ) )
            {
                createdInterpolator = std::make_shared< LagrangeInterpolator
                        < IndependentVariableType, DependentVariableType, double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ),
                            lagrangeInterpolatorSettings->getLagrangeBoundaryHandling( ),
                            interpolatorSettings->getBoundaryHandling( ).at( 0 ), defaultExtrapolationValue );
            }
            else
            {
                createdInterpolator = std::make_shared< LagrangeInterpolator
                        < IndependentVariableType, DependentVariableType, long double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ),
                            lagrangeInterpolatorSettings->getLagrangeBoundaryHandling( ),
                            interpolatorSettings->getBoundaryHandling( ).at( 0 ), defaultExtrapolationValue );
            }
        }
        else
        {
            throw std::runtime_error( "Error, did not recognize lagrange interpolator settings" );

        }
        break;
    }
    case hermite_spline_interpolator:
    {
        if( firstDerivativeOfDependentVariables.size( ) != dataToInterpolate.size( ) )
        {
            throw std::runtime_error(
                        "Error when creating hermite spline interpolator, derivative size is inconsistent" );
        }
        createdInterpolator = std::make_shared< HermiteCubicSplineInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, firstDerivativeOfDependentVariables,
                    interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryHandling( ).at( 0 ), defaultExtrapolationValue );
        break;
    }
    case piecewise_constant_interpolator:
        createdInterpolator = std::make_shared< PiecewiseConstantInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryHandling( ).at( 0 ), defaultExtrapolationValue );
        break;
    default:
        throw std::runtime_error( "Error when making interpolator, function cannot be used to create interplator of type " +
                                  std::to_string( interpolatorSettings->getInterpolatorType( ) ) );
    }
    return createdInterpolator;
}

//! Function to create an interpolator from DataInterpolationSettings
/*!
 *  Function to create an interpolator from DataInterpolationSettings
 *  \param dataInterpolationSettings Object containing the data that is to be interpolated and settings that are to be
 *      used to create the interpolator.
 *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
        of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
 *  \return Interpolator created from dataToInterpolate using interpolatorSettings.
 */
template< typename IndependentType, typename DependentType >
std::shared_ptr< OneDimensionalInterpolator< IndependentType, DependentType > > createOneDimensionalInterpolator(
        const std::shared_ptr< DataInterpolationSettings< IndependentType, DependentType > > dataInterpolationSettings,
        const std::pair< DependentType, DependentType >& defaultExtrapolationValue =
        std::make_pair( IdentityElement::getAdditionIdentity< DependentType >( ),
                        IdentityElement::getAdditionIdentity< DependentType >( ) ) )
{
    std::vector< DependentType > firstDerivativeOfDependentVariables;
    std::shared_ptr< HermiteDataSettings< IndependentType, DependentType > > hermiteDataSettings =
            std::dynamic_pointer_cast< HermiteDataSettings< IndependentType, DependentType > >(
                dataInterpolationSettings->dataSettings_ );
    if ( hermiteDataSettings )
    {
        firstDerivativeOfDependentVariables = hermiteDataSettings->firstDerivativeOfDependentVariables_;
    }
    return createOneDimensionalInterpolator( dataInterpolationSettings->dataSettings_->getDataMap( ),
                                             dataInterpolationSettings->interpolatorSettings_,
                                             defaultExtrapolationValue,
                                             firstDerivativeOfDependentVariables );
}

//! Function to create a multi-dimensional interpolator
/*!
 *  Function to create a multi-dimensional interpolator from the data that is to be interpolated,
 *  as well as the settings that are to be used to create the interpolator.
 *  \param independentValues Vector of vectors containing data points of independent variables,
 *      each must be sorted in ascending order.
 *  \param dependentData Multi-dimensional array of dependent data at each point of
 *      hyper-rectangular grid formed by independent variable points.
 *  \param interpolatorSettings Settings that are to be used to create interpolator.
 *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
        of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
 *  \return Interpolator created from independentValues and dependentData using interpolatorSettings.
 */
template< typename IndependentVariableType, typename DependentVariableType, unsigned int NumberOfDimensions >
std::shared_ptr< MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions > >
createMultiDimensionalInterpolator(
        const std::vector< std::vector< IndependentVariableType > >& independentValues,
        const boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) >& dependentData,
        const std::shared_ptr< InterpolatorSettings > interpolatorSettings,
        const std::vector< std::pair< DependentVariableType, DependentVariableType > >& defaultExtrapolationValue =
        std::vector< std::pair< DependentVariableType, DependentVariableType > >
        ( NumberOfDimensions, std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                                              IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ) )
{
    std::shared_ptr< MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions > >
            createdInterpolator;

    // Assign default values
    std::vector< BoundaryInterpolationType > boundaryHandlingVector = interpolatorSettings->getBoundaryHandling( );
    if ( boundaryHandlingVector.empty( ) )
    {
        boundaryHandlingVector = std::vector< BoundaryInterpolationType >( NumberOfDimensions, extrapolate_at_boundary );
    }
    // Check size of boundary handling methods
    else if ( boundaryHandlingVector.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error while creating multi-dimensional interpolator. The number of boundary handling methods does not "
                                  "match the number of dimensions." );
    }

    // Check type of interpolator.
    switch ( interpolatorSettings->getInterpolatorType( ) )
    {
    case multi_linear_interpolator:
    {
        createdInterpolator = std::make_shared< MultiLinearInterpolator
                < IndependentVariableType, DependentVariableType, NumberOfDimensions > >(
                    independentValues, dependentData, interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryHandling( ), defaultExtrapolationValue );
        break;
    }
    default:
        throw std::runtime_error( "Error when making interpolator, function cannot be used to create interplator of type " +
                                  std::to_string( interpolatorSettings->getInterpolatorType( ) ) );
    }
    return createdInterpolator;
}

template< typename TimeType, typename StateScalarType, int InputRows, int InputColumns, int OutputRows, int OutputColumns >
std::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, OutputRows, OutputColumns > > >
convertBetweenStaticDynamicEigenTypeInterpolators(
        const std::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, InputRows, InputColumns > > > inputInterpolator )
{
    if( ( InputRows > 0 && OutputRows > 0 && InputRows != OutputRows ) ||
            ( InputColumns > 0 && OutputColumns > 0 && InputColumns != OutputColumns ) )
    {
        throw std::runtime_error( "Error when converting interpolator Eigen type; sizes are inconsistent, input columns, "
                                  "cannot convert between different fixed sizes" );
    }
    typedef Eigen::Matrix< StateScalarType, InputRows, InputColumns > InputState;
    typedef Eigen::Matrix< StateScalarType, OutputRows, OutputColumns > OutputState;

    std::shared_ptr< OneDimensionalInterpolator< TimeType, OutputState > > outputInterpolator;

    std::vector< TimeType > times = inputInterpolator->getIndependentValues( );
    std::vector< InputState > inputStates = inputInterpolator->getDependentValues( );
    std::vector< OutputState > outputStates;

    if( OutputColumns >= 0 && inputStates.at( 0 ).cols( ) != OutputColumns  )
    {
        throw std::runtime_error( "Error when converting interpolator Eigen type; sizes are inconsistent, input columns: " +
                                  std::to_string( inputStates.at( 0 ).cols( ) ) + "; output columns:" +
                                  std::to_string( OutputColumns ) );
    }
    else
    {
        for( unsigned int i = 0; i < inputStates.size( ); i++ )
        {
            outputStates.push_back( inputStates.at( i ) );
        }
        std::pair< InputState, InputState > inputDefaultExtrapolationValues =
                inputInterpolator->getDefaultExtrapolationValue( );

        std::pair< OutputState, OutputState > outputDefaultExtrapolationValues;
        if( inputDefaultExtrapolationValues.first.rows( ) > 0 &&
                inputDefaultExtrapolationValues.first.cols( ) > 0 )
        {
            outputDefaultExtrapolationValues =
                    std::make_pair( inputDefaultExtrapolationValues.first, inputDefaultExtrapolationValues.second );
        }
        else
        {
            outputDefaultExtrapolationValues = std::make_pair(
                        IdentityElement::getAdditionIdentity< OutputState >( ),
                        IdentityElement::getAdditionIdentity< OutputState >( ) );
        }

        switch ( inputInterpolator->getInterpolatorType( ) )
        {
        case linear_interpolator:
            outputInterpolator = std::make_shared< LinearInterpolator< TimeType, OutputState > >(
                        times, outputStates,
                        inputInterpolator->getSelectedLookupScheme( ),
                        inputInterpolator->getBoundaryHandling( ),
                        outputDefaultExtrapolationValues );
            break;
        case cubic_spline_interpolator:

            outputInterpolator = std::make_shared< CubicSplineInterpolator< TimeType, OutputState > >(
                        times, outputStates,
                        inputInterpolator->getSelectedLookupScheme( ),
                        inputInterpolator->getBoundaryHandling( ),
                        outputDefaultExtrapolationValues );

            break;
        case piecewise_constant_interpolator:
            outputInterpolator = std::make_shared< PiecewiseConstantInterpolator< TimeType, OutputState > >(
                        times, outputStates, inputInterpolator->getSelectedLookupScheme( ),
                        inputInterpolator->getBoundaryHandling( ),
                        outputDefaultExtrapolationValues );
            break;
        case lagrange_interpolator:
        {
            std::shared_ptr< LagrangeInterpolator< TimeType, Eigen::Matrix< StateScalarType, InputRows, InputColumns > > > lagrangeInputInterpolator =
                    std::dynamic_pointer_cast< LagrangeInterpolator< TimeType, Eigen::Matrix< StateScalarType, InputRows, InputColumns > > >( inputInterpolator );
            if( lagrangeInputInterpolator == nullptr )
            {
                throw std::runtime_error( "Error when changing vector size type of Lagrange interpolator; input type is inconsistent" );
            }
            else
            {
                outputInterpolator = std::make_shared< LagrangeInterpolator< TimeType, OutputState > >(
                            times, outputStates,
                            lagrangeInputInterpolator->getNumberOfStages( ),
                            inputInterpolator->getSelectedLookupScheme( ),
                            lagrangeInputInterpolator->getLagrangeBoundaryHandling( ),
                            inputInterpolator->getBoundaryHandling( ),
                            outputDefaultExtrapolationValues );
            }
            break;
        }
        case hermite_spline_interpolator:
        {
            std::shared_ptr< HermiteCubicSplineInterpolator< TimeType, Eigen::Matrix< StateScalarType, InputRows, InputColumns > > > hermiteInputInterpolator =
                    std::dynamic_pointer_cast< HermiteCubicSplineInterpolator< TimeType, Eigen::Matrix< StateScalarType, InputRows, InputColumns > > >( inputInterpolator );
            if( hermiteInputInterpolator == nullptr )
            {
                throw std::runtime_error( "Error when changing vector size type of Hermite interpolator; input type is inconsistent" );
            }
            else
            {
                std::vector< InputState > inputStateDerivatives = hermiteInputInterpolator->getDerivativeValues( );
                std::vector< OutputState > outputStateDerivatives;
                for( unsigned int i = 0; i < inputStateDerivatives.size( ); i++ )
                {
                    outputStateDerivatives.push_back( inputStateDerivatives.at( i ) );
                }

                outputInterpolator = std::make_shared< HermiteCubicSplineInterpolator< TimeType, OutputState > >(
                            times, outputStates, outputStateDerivatives, inputInterpolator->getSelectedLookupScheme( ),
                            inputInterpolator->getBoundaryHandling( ),
                            outputDefaultExtrapolationValues );
            }
            break;
        }
        default:
            throw std::runtime_error( "Error when changing vector size type of interpolator; input  interpolator type is not recognized" );
            break;
        }
    }
    return outputInterpolator;
}
} // namespace interpolators

} // namespace tudat

#endif // TUDAT_CREATEINTERPOLATOR_H
