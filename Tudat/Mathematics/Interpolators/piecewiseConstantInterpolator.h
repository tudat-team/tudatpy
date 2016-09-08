#ifndef PIECEWISECONSTANTINTERPOLATOR_H
#define PIECEWISECONSTANTINTERPOLATOR_H

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

namespace interpolators
{

template< typename IndependentVariableType, typename DependentVariableType >
class PiecewiseConstantInterpolator : public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::lookUpScheme_;

    PiecewiseConstantInterpolator( const std::vector< IndependentVariableType > independentVariables,
                                   const std::vector< DependentVariableType > dependentVariables,
                                   const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
    {
        independentValues_ = independentVariables;
        dependentValues_ = dependentVariables;

        if( dependentValues_.size( ) != independentValues_.size( ) )
        {
            std::cerr<<"Warning: independent and dependent variables not of same size in piecewise constant interpolator constrcutor"<<std::endl;
        }

        this->makeLookupScheme( selectedLookupScheme );
    }

    PiecewiseConstantInterpolator( const std::map< IndependentVariableType, DependentVariableType > dataMap,
                                   const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
    {
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );

        int counter = 0;
        for( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator mapIterator = dataMap.begin( );
             mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_[ counter ] = mapIterator->first;
            dependentValues_[ counter ] = mapIterator->second;
            counter++;
        }
        this->makeLookupScheme( selectedLookupScheme );
    }

    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;


    DependentVariableType interpolate( const IndependentVariableType targetIndependentVariableValue )
    {
        // Determine the lower entry in the table corresponding to the target independent variable value.
        int lowerEntry;
        if( targetIndependentVariableValue <= independentValues_.at( 0 ) )
        {
            lowerEntry = 0;
        }
        else if( targetIndependentVariableValue >= independentValues_.at( independentValues_.size( ) -1 ) )
        {
            lowerEntry = independentValues_.size( ) - 1;
        }
        else
        {
            lowerEntry = lookUpScheme_->findNearestLowerNeighbour( targetIndependentVariableValue );
        }

        return dependentValues_.at( lowerEntry );
    }

protected:

private:

};

}

}

#endif // PIECEWISECONSTANTINTERPOLATOR_H
