#ifndef CREATEINTERPOLATOR_H
#define CREATEINTERPOLATOR_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace interpolators
{

enum OneDimensionalInterpolatorTypes
{
    linear_interpolator = 1,
    cubic_spline_interpolator = 2,
    lagrange_interpolator = 3
};

class InterpolatorSettings
{
public:

    InterpolatorSettings( const OneDimensionalInterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme ):
        interpolatorType_( interpolatorType ), selectedLookupScheme_( selectedLookupScheme ){ }

    virtual ~InterpolatorSettings( ){ }

    OneDimensionalInterpolatorTypes getInterpolatorType( )
    {
        return interpolatorType_;
    }

    AvailableLookupScheme getSelectedLookupScheme( )
    {
        return selectedLookupScheme_;
    }

protected:

    OneDimensionalInterpolatorTypes interpolatorType_;

    AvailableLookupScheme selectedLookupScheme_;

};

class LagrangeInterpolatorSettings: public InterpolatorSettings
{
public:
    LagrangeInterpolatorSettings( const int interpolatorOrder,
                                  const bool useLongDoubleTimeStep,
                                  const AvailableLookupScheme selectedLookupScheme ):
        InterpolatorSettings( lagrange_interpolator, selectedLookupScheme ),
        interpolatorOrder_( interpolatorOrder ), useLongDoubleTimeStep_( useLongDoubleTimeStep )
        { }

    ~LagrangeInterpolatorSettings( ){ }

    int getInterpolatorOrder( )
    {
        return interpolatorOrder_;
    }

    bool getUseLongDoubleTimeStep( )
    {
        return useLongDoubleTimeStep_;
    }


protected:

    int interpolatorOrder_;

    bool useLongDoubleTimeStep_;



};


template< typename IndependentVariableType, typename DependentVariableType >
boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > > createOneDimensionalInterpolator(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const boost::shared_ptr< InterpolatorSettings > interpolatorSettings )
{
    boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > > createdInterpolator;
    switch( interpolatorSettings->getInterpolatorType( ) )
    {
    case linear_interpolator:
        createdInterpolator = boost::make_shared< LinearInterpolator< IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        break;
    case cubic_spline_interpolator:
        createdInterpolator = boost::make_shared< CubicSplineInterpolator< IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ) );
        break;
    case lagrange_interpolator:
    {
        boost::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                boost::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );

        if( lagrangeInterpolatorSettings != NULL )
        {
            if( !lagrangeInterpolatorSettings->getUseLongDoubleTimeStep( ) )
            {
                createdInterpolator = boost::make_shared< LagrangeInterpolator< IndependentVariableType, DependentVariableType, double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ) );
            }
            else
            {
                createdInterpolator = boost::make_shared< LagrangeInterpolator< IndependentVariableType, DependentVariableType, long double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ) );
            }
        }
        else
        {
            std::cerr<<"Error, did not recognize lagrange interpolator settings"<<std::endl;

        }
        break;
    }
    default:
        std::cerr<<"Error when making interpolator, function cannot be used to create interplator of type "<<interpolatorSettings->getInterpolatorType( )<<std::endl;
    }
    return createdInterpolator;
}

}

}

#endif // CREATEINTERPOLATOR_H
