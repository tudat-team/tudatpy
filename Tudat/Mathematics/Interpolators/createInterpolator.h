/*    Copyright (c) 2010-2016, Delft University of Technology
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
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace interpolators
{

//! Enum of available interpolator types.
enum OneDimensionalInterpolatorTypes
{
    linear_interpolator = 1,
    cubic_spline_interpolator = 2,
    lagrange_interpolator = 3
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
     */
    InterpolatorSettings( const OneDimensionalInterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm ):
        interpolatorType_( interpolatorType ), selectedLookupScheme_( selectedLookupScheme ){ }

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

protected:

    //! Selected type of interpolator.
    OneDimensionalInterpolatorTypes interpolatorType_;

    //! Selected type of lookup scheme for independent variables.
    AvailableLookupScheme selectedLookupScheme_;

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
            const LagrangeInterpolatorBoundaryHandling boundaryHandling = lagrange_cubic_spline_boundary_interpolation ):
        InterpolatorSettings( lagrange_interpolator, selectedLookupScheme ),
        interpolatorOrder_( interpolatorOrder ), useLongDoubleTimeStep_( useLongDoubleTimeStep ),
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

    //! Function to get a boolean denoting whether time step is to be a long double
    /*!
     * Function to get a boolean denoting whether time step is to be a long double
     * \return Boolean denoting whether time step is to be a long double
     */
    bool getUseLongDoubleTimeStep( )
    {
        return useLongDoubleTimeStep_;
    }

    LagrangeInterpolatorBoundaryHandling getBoundaryHandling( )
    {
        return boundaryHandling_;
    }


protected:

    //! Order of the Lagrange interpolator that is to be created.
    int interpolatorOrder_;

    //!  Boolean denoting whether time step is to be a long double.
    bool useLongDoubleTimeStep_;

    LagrangeInterpolatorBoundaryHandling boundaryHandling_;

};

//! Function to create an interpolator
/*!
 *  Function to create an interpolator from the data that is to be interpolated, as well as the
 *  settings that are to be used to create the interpolator.
 *  \param dataToInterpolate Map providing data that is to be interpolated (key = independent
 *  variables, value = dependent variables)
 *  \param interpolatorSettings Settings that are to be used to create interpolator
 *  \return Interpolator created from dataToInterpolate using interpolatorSettings.
 */
template< typename IndependentVariableType, typename DependentVariableType >
boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
createOneDimensionalInterpolator(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const boost::shared_ptr< InterpolatorSettings > interpolatorSettings )
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
        createdInterpolator = boost::make_shared< CubicSplineInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        break;
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
    default:
        throw std::runtime_error(
                    "Error when making interpolator, function cannot be used to create interplator of type " +
                    boost::lexical_cast< std::string >(
                        interpolatorSettings->getInterpolatorType( ) ) );
    }
    return createdInterpolator;
}

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_CREATEINTERPOLATOR_H
