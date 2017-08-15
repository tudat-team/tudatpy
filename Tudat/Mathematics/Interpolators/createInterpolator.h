/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/hermiteCubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"

#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{

namespace interpolators
{

//! Enum of available interpolator types.
enum OneDimensionalInterpolatorTypes
{
    linear_interpolator = 1,
    cubic_spline_interpolator = 2,
    lagrange_interpolator = 3,
    hermite_spline_interpolator = 4,
    piecewise_constant_interpolator = 5
};

//! Base class for providing settings for creating an interpolator.
/*!
 *  Base class for providing settings for creating an interpolator using the createInterpolator
 *  function. This base class is functional, i.e. may be used directly for interpolator types
 *  requiring no additional information. For interpolators that do require more information,
 *  a derived class is provided.
 */
class InterpolatorSettings
{
public:

    //! Constructor
    /*!
     *  Constructor
     * \param interpolatorType Selected type of interpolator.
     * \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     * \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     * time step is a double if false.
     */
    InterpolatorSettings( const OneDimensionalInterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                          const bool useLongDoubleTimeStep = 0 ):
        interpolatorType_( interpolatorType ), selectedLookupScheme_( selectedLookupScheme ),
        useLongDoubleTimeStep_( useLongDoubleTimeStep ){ }

    //! Virtual destructor
    virtual ~InterpolatorSettings( ){ }

    //! Function to get the selected type of interpolator.
    /*!
     * Function to get the selected type of interpolator.
     * \return Selected type of interpolator.
     */
    OneDimensionalInterpolatorTypes getInterpolatorType( )
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


    //! Function to get a boolean denoting whether time step is to be a long double
    /*!
     * Function to get a boolean denoting whether time step is to be a long double
     * \return Boolean denoting whether time step is to be a long double
     */
    bool getUseLongDoubleTimeStep( )
    {
        return useLongDoubleTimeStep_;
    }

protected:

    //! Selected type of interpolator.
    OneDimensionalInterpolatorTypes interpolatorType_;

    //! Selected type of lookup scheme for independent variables.
    AvailableLookupScheme selectedLookupScheme_;

    //!  Boolean denoting whether time step is to be a long double.
    bool useLongDoubleTimeStep_;

};

//! Class for providing settings to creating a Lagrange interpolator.
class LagrangeInterpolatorSettings: public InterpolatorSettings
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param interpolatorOrder Order of the Lagrange interpolator that is to be created.
     * \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     * time step is a double if false.
     * \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     * \param boundaryHandling Variable denoting the method by which the boundary interpolation is handled.
     */
    LagrangeInterpolatorSettings(
            const int interpolatorOrder,
            const bool useLongDoubleTimeStep = 0,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const LagrangeInterpolatorBoundaryHandling boundaryHandling =
            lagrange_cubic_spline_boundary_interpolation ):
        InterpolatorSettings( lagrange_interpolator, selectedLookupScheme, useLongDoubleTimeStep ),
        interpolatorOrder_( interpolatorOrder ),
        boundaryHandling_( boundaryHandling )
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

    LagrangeInterpolatorBoundaryHandling getBoundaryHandling( )
    {
        return boundaryHandling_;
    }


protected:

    //! Order of the Lagrange interpolator that is to be created.
    int interpolatorOrder_;

    LagrangeInterpolatorBoundaryHandling boundaryHandling_;

};


//! -DOC
template< typename IndependentType, typename DependentType >
class DataMapSettings
{
public:
    //! -DOC
    DataMapSettings( ) { }

    //! -DOC
    DataMapSettings( const std::map< IndependentType, DependentType >& dataMap ) : dataMap_( dataMap ) { }

    //! Virtual destructor
    virtual ~DataMapSettings( ) { }

    //! -DOC
    virtual std::map< IndependentType, DependentType > getDataMap( ) const
    {
        return dataMap_;
    }

protected:
    //! -DOC
    std::map< IndependentType, DependentType > dataMap_;
};

//! -DOC
template< typename IndependentType, typename DependentType >
class IndependentDependentDataMapSettings : public DataMapSettings< IndependentType, DependentType >
{
public:
    //! -DOC
    IndependentDependentDataMapSettings( const std::vector< IndependentType >& independentVariableValues,
                                         const std::vector< DependentType >& dependentVariableValues ) :
        DataMapSettings< IndependentType, DependentType >( ),
        independentVariableValues_( independentVariableValues ),
        dependentVariableValues_( dependentVariableValues ) { }

    //! Virtual destructor
    virtual ~IndependentDependentDataMapSettings( ) { }

    //! -DOC
    std::vector< IndependentType > independentVariableValues_;

    //! -DOC
    std::vector< DependentType > dependentVariableValues_;

    //! -DOC
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

//! -DOC
template< typename EigenVectorType >
class FromFileDataMapSettings : public DataMapSettings< typename EigenVectorType::Scalar, EigenVectorType >
{
public:
    //! -DOC
    FromFileDataMapSettings( const std::string& relativeFilePath ) :
        DataMapSettings< typename EigenVectorType::Scalar, EigenVectorType >( ),
        relativeFilePath_( relativeFilePath ) { }

    //! Virtual destructor
    virtual ~FromFileDataMapSettings( ) { }

    //! -DOC
    std::string relativeFilePath_;

    //! -DOC
    virtual std::map< typename EigenVectorType::Scalar, EigenVectorType > getDataMap( ) const
    {
        return input_output::readEigenVectorMapFromFile< EigenVectorType >( relativeFilePath_ );
    }
};

//! -DOC
template< typename IndependentType, typename DependentType >
class HermiteDataSettings : public DataMapSettings< IndependentType, DependentType >
{
public:
    //! -DOC
    HermiteDataSettings( const std::map< IndependentType, DependentType >& dataToInterpolate,
                         const std::vector< DependentType >& firstDerivativeOfDependentVariables ) :
        DataMapSettings< IndependentType, DependentType >( dataToInterpolate ),
        firstDerivativeOfDependentVariables_( firstDerivativeOfDependentVariables ) { }

    //! Virtual destructor
    virtual ~HermiteDataSettings( ) { }

    //! -DOC
    std::vector< DependentType > firstDerivativeOfDependentVariables_;
};


//! -DOC
template< typename IndependentType, typename DependentType >
class DataInterpolationSettings
{
public:
    //! -DOC
    DataInterpolationSettings(
            const boost::shared_ptr< DataMapSettings< IndependentType, DependentType > >& dataSettings,
            const boost::shared_ptr< InterpolatorSettings >& interpolatorSettings ) :
        dataSettings_( dataSettings ), interpolatorSettings_( interpolatorSettings ) { }

    //! Virtual destructor
    virtual ~DataInterpolationSettings( ){ }

    //! -DOC
    boost::shared_ptr< DataMapSettings< IndependentType, DependentType > > dataSettings_;

    //! -DOC
    boost::shared_ptr< InterpolatorSettings > interpolatorSettings_;
};


//! Function to create an interpolator
/*!
 *  Function to create an interpolator from the data that is to be interpolated, as well as the
 *  settings that are to be used to create the interpolator.
 *  \param dataToInterpolate Map providing data that is to be interpolated (key = independent
 *  variables, value = dependent variables)
 *  \param interpolatorSettings Settings that are to be used to create interpolator
 *  \param firstDerivativeOfDependentVariables First derivative of dependent variables w.r.t. independent variable at
 *  independent variables values in values of dataToInterpolate. By default, this vector is empty, it only needs to
 *  be supplied if the selected interpolator requires this data (e.g. Hermite spline).
 *  \return Interpolator created from dataToInterpolate using interpolatorSettings.
 */
template< typename IndependentVariableType, typename DependentVariableType >
boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
createOneDimensionalInterpolator(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const boost::shared_ptr< InterpolatorSettings > interpolatorSettings,
        const std::vector< DependentVariableType > firstDerivativeOfDependentVariables =
        std::vector< DependentVariableType >( ) )
{
    boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
            createdInterpolator;

    // Check type of interpolator.
    switch( interpolatorSettings->getInterpolatorType( ) )
    {
    case linear_interpolator:
        createdInterpolator = boost::make_shared< LinearInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        break;
    case cubic_spline_interpolator:
    {
        if( !interpolatorSettings->getUseLongDoubleTimeStep( ) )
        {
            createdInterpolator = boost::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType > >(
                        dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        }
        else
        {
            createdInterpolator = boost::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType, long double > >(
                        dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        }
        break;
    }
    case lagrange_interpolator:
    {
        // Check consistency of input
        boost::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                boost::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );
        if( lagrangeInterpolatorSettings != NULL )
        {
            // Create Lagrange interpolator with requested time step type
            if( !lagrangeInterpolatorSettings->getUseLongDoubleTimeStep( ) )
            {
                createdInterpolator = boost::make_shared< LagrangeInterpolator
                        < IndependentVariableType, DependentVariableType, double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ),
                            lagrangeInterpolatorSettings->getBoundaryHandling( ) );
            }
            else
            {
                createdInterpolator = boost::make_shared< LagrangeInterpolator
                        < IndependentVariableType, DependentVariableType, long double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ),
                            lagrangeInterpolatorSettings->getBoundaryHandling( ) );
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
        createdInterpolator = boost::make_shared< HermiteCubicSplineInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, firstDerivativeOfDependentVariables,
                    interpolatorSettings->getSelectedLookupScheme( ) );
        break;
    }
    case piecewise_constant_interpolator:
        createdInterpolator = boost::make_shared< PiecewiseConstantInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        break;
    default:
        throw std::runtime_error(
                    "Error when making interpolator, function cannot be used to create interplator of type " +
                    boost::lexical_cast< std::string >(
                        interpolatorSettings->getInterpolatorType( ) ) );
    }
    return createdInterpolator;
}

//! -DOC
template< typename IndependentType, typename DependentType >
boost::shared_ptr< OneDimensionalInterpolator< IndependentType, DependentType > > createOneDimensionalInterpolator(
        const boost::shared_ptr< DataInterpolationSettings< IndependentType, DependentType > >
        dataInterpolationSettings )
{
    std::vector< DependentType > firstDerivativeOfDependentVariables;
    boost::shared_ptr< HermiteDataSettings< IndependentType, DependentType > > hermiteDataSettings =
            boost::dynamic_pointer_cast< HermiteDataSettings< IndependentType, DependentType > >(
                dataInterpolationSettings->dataSettings_ );
    if ( hermiteDataSettings )
    {
        firstDerivativeOfDependentVariables = hermiteDataSettings->firstDerivativeOfDependentVariables_;
    }
    return createOneDimensionalInterpolator( dataInterpolationSettings->dataSettings_->getDataMap( ),
                                             dataInterpolationSettings->interpolatorSettings_,
                                             firstDerivativeOfDependentVariables );
}

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_CREATEINTERPOLATOR_H
