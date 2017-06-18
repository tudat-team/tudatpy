#ifndef SLERPINTERPOLATOR_H
#define SLERPINTERPOLATOR_H

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

namespace interpolators
{

template< typename IndependentVariableType, typename ScalarType = IndependentVariableType >
class SlerpInterpolator : public OneDimensionalInterpolator< IndependentVariableType, Eigen::Quaterniond >
{
public:

    //! Using statements to prevent having to put 'this' everywhere in the code.
    using OneDimensionalInterpolator< IndependentVariableType, Eigen::Quaterniond >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, Eigen::Quaterniond >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, Eigen::Quaterniond >::lookUpScheme_;


    //! Constructor from vectors of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from two vectors containing the independent
     *  variables and dependent variables. A look-up scheme can be provided to
     *  override the given default.
     *  \param independentValues. Vector of values of independent variables that are used.
     *  \param dependentValues. Vector of values of dependent variables that are used.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *  to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     */
    SlerpInterpolator( const std::vector< IndependentVariableType >& independentVariables,
                       const std::vector< Eigen::Quaterniond >& dependentVariables,
                       const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
    {
        // Set data vectors.
        independentValues_ = independentVariables;
        dependentValues_ = dependentVariables;

        // Verify that the initialization variables are not empty.
        if ( independentValues_.size( ) == 0 || dependentValues_.size( ) == 0 )
        {
            boost::throw_exception( boost::enable_error_info( std::runtime_error(
                                                                  "The vectors used in the slerp interpolator initialization are empty." ) ) );
        }

        // Check consistency of input data.
        if( static_cast< int >( independentValues_.size( ) ) != dependentValues_.size( ) )
        {
            boost::throw_exception( boost::enable_error_info( std::runtime_error(
                                                                  "Warning: independent and slerp variables not of same size in lagrange interpolator constrcutor" ) ) );
        }

        // Create lookup scheme from independent variable values
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Constructor from map of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from a map containing independent variables
     *  as key and dependent variables as value. A look-up scheme can be provided to override the
     *  given default.
     *  \param dataMap. Map containing independent variables as key and dependent variables as value
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *  to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     */
    SlerpInterpolator( const std::map< IndependentVariableType, Eigen::Quaterniond >& dataMap,
                       const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
    {

        // Verify that the initialization variables are not empty.
        if ( dataMap.size( ) == 0 )
        {
            boost::throw_exception( boost::enable_error_info( std::runtime_error(
                                                                  "The vectors used in the slerp interpolator initialization are empty." ) ) );
        }

        // Resize data vectors of independent/dependent values.
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );

        // Fill data vectors with data from map.
        int counter = 0;
        for( typename std::map< IndependentVariableType, Eigen::Quaterniond >::const_iterator mapIterator = dataMap.begin( );
             mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_[ counter ] = mapIterator->first;
            dependentValues_[ counter ] = mapIterator->second;
            counter++;
        }

        // Create lookup scheme from independent variable data points.
        this->makeLookupScheme( selectedLookupScheme );
    }

    // Using statement to prevent compiler warning.
    using Interpolator< IndependentVariableType, Eigen::Quaterniond >::interpolate;

    //! Function interpolates dependent variable value at given independent variable value.
    /*!
     *  Function interpolates dependent variable value at given independent variable value.
     *  \param independentVariableValue Value of independent variable at which interpolation is to take place.
     *  \return Interpolated value of dependent variable.
     */
    Eigen::Quaterniond interpolate( const IndependentVariableType targetIndependentVariableValue )
    {
        // Check if requested time is within time range specified by tabulated data.
        if( targetIndependentVariableValue < independentValues_[ 0 ] ||
                targetIndependentVariableValue > independentValues_[ independentValues_.size( ) - 1 ] )
        {
            std::cerr<<"Error when performing quaternion interpolation, requested time outsie available range "<<
                       targetIndependentVariableValue<<" "<<independentValues_[ 0 ]<<" "<<independentValues_[ independentValues_.size( ) - 1 ]<<std::endl;
        }

        // Get interval of tabulated data which requested time is  in
        int nearestLowerIndex = lookUpScheme_->findNearestLowerNeighbour( targetIndependentVariableValue );
        IndependentVariableType lowerIndependentValue = independentValues_[ nearestLowerIndex ];
        IndependentVariableType upperIndependentValue = independentValues_[ nearestLowerIndex + 1 ];

        // Get fraction in curretn interval where for requested time
        previousDifference_ = static_cast< ScalarType >( targetIndependentVariableValue - lowerIndependentValue );
        ScalarType fractionInCurrentStep =
                previousDifference_ / static_cast< ScalarType >( upperIndependentValue - lowerIndependentValue );


        // Perform slerp interpolation to get rotation value at requested time and return.
        return dependentValues_[ nearestLowerIndex ].slerp( fractionInCurrentStep, dependentValues_[ nearestLowerIndex + 1 ] );
    }

    ScalarType getPreviousDifference( )
    {
        return previousDifference_;
    }

protected:

private:
    ScalarType previousDifference_;
};

}

}

#endif // SLERPINTERPOLATOR_H
