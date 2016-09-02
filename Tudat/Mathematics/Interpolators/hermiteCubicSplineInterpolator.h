#ifndef TUDAT_HERMITE_CUBIC_SPLINE_INTERPOLATOR_H
#define TUDAT_HERMITE_CUBIC_SPLINE_INTERPOLATOR_H


#include <iostream> // cout sometimes needs this

#include <Eigen/Core>
#include <vector>
#include <tudat/Mathematics/Interpolators/interpolator.h>
#include <tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>

#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"


namespace tudat
{

namespace interpolators
{

//! Hermite Cubic Spline Interpolator
template< typename VariableType >
class HermiteCubicSplineInterpolator :
        public OneDimensionalInterpolator< VariableType, VariableType >
{
public:

    using OneDimensionalInterpolator< VariableType, VariableType >::
    dependentValues_;
    using OneDimensionalInterpolator< VariableType, VariableType >::
    independentValues_;
    using OneDimensionalInterpolator< VariableType, VariableType >::
    lookUpScheme_;

    //! Constructor
    HermiteCubicSplineInterpolator(
            std::vector< VariableType > independentValues,
            std::vector< VariableType > dependentValues,
            std::vector< VariableType > derivativeValues,
            AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
    {
        // save data
        independentValues_ = independentValues   ;
        dependentValues_ = dependentValues   ;
        derivativeValues_  = derivativeValues    ;

        // compute coefficients
        computeCoefficients();

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Destructor
    ~HermiteCubicSplineInterpolator( ){ }

    //! Get coefficients
    std::vector< std::vector< VariableType > > GetCoefficients( )
    {
        return coefficients_;
    }

    //! Interpolate
    VariableType interpolate( const VariableType targetIndependentVariableValue )
    {
        // Determine the lower entry in the table corresponding to the target independent variable
        // value.
        int lowerEntry_ = lookUpScheme_->findNearestLowerNeighbour(
                    targetIndependentVariableValue );

        // p(x) = a((x-x0)/(x1-x0))^3 + b((x-x0)/(x1-x0))^2 + c((x-x0)/(x1-x0)) + d
        double factor = ( targetIndependentVariableValue - independentValues_[lowerEntry_] )
                /( independentValues_[lowerEntry_+1] - independentValues_[lowerEntry_] );

        VariableType targetValue =
                coefficients_[0][ lowerEntry_ ] * factor * factor * factor
                + coefficients_[1][ lowerEntry_ ] * factor * factor
                + coefficients_[2][ lowerEntry_ ] * factor
                + coefficients_[3][ lowerEntry_ ] ;

        return targetValue;
    }

protected:

    //! Compute coefficients of the splines
    void computeCoefficients()
    {
        // Initialize vector
        std::vector< VariableType > zeroVect( independentValues_.size() - 1 );

        for( int i = 0 ; i < 4 ; i++ )
        {
            coefficients_.push_back( zeroVect );
        }

        for( unsigned int i = 0 ; i < ( independentValues_.size() - 1 ) ; i++ )
        {
            // p(x) = a((x-x0)/(x1-x0))^3 + b((x-x0)/(x1-x0))^2 + c((x-x0)/(x1-x0)) + d
            // a
            coefficients_[0][i] = 2.0*dependentValues_[i] - 2.0*dependentValues_[i+1] + derivativeValues_[i]*(independentValues_[i+1]-independentValues_[i]) + derivativeValues_[i+1]*(independentValues_[i+1]-independentValues_[i])        ;

            // b
            coefficients_[1][i] = -3.0*dependentValues_[i] + 3.0*dependentValues_[i+1] - 2.0*derivativeValues_[i]*(independentValues_[i+1]-independentValues_[i]) - derivativeValues_[i+1]*(independentValues_[i+1]-independentValues_[i])   ;

            // c
            coefficients_[2][i] = derivativeValues_[i]*(independentValues_[i+1]-independentValues_[i])    ;

            // d
            coefficients_[3][i] = dependentValues_[i]   ;
        }
    }


private:

    //! Derivatives of dependent variable to independent variable
    std::vector< VariableType > derivativeValues_ ;

    //! Coefficients of splines
    std::vector< std::vector< VariableType > > coefficients_ ;
};

typedef HermiteCubicSplineInterpolator< double > HermiteCubicSplineInterpolatorDouble;

} // Close Namespace Interpolators


} // Close Namespace tudat

#endif // TUDAT_HERMITE_CUBIC_SPLINE_INTERPOLATOR_H
