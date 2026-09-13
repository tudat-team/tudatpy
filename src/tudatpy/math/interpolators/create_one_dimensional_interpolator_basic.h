#ifndef TUDATPY_CREATE_ONE_DIMENSIONAL_INTERPOLATOR_BASIC_H
#define TUDATPY_CREATE_ONE_DIMENSIONAL_INTERPOLATOR_BASIC_H

#include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{
namespace interpolators
{

template< typename IndependentVariableType, typename DependentVariableType >
std::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > > createOneDimensionalInterpolatorBasic(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const std::shared_ptr< InterpolatorSettings > interpolatorSettings,
        const std::vector< DependentVariableType > firstDerivativesOfDataToIntepolate = std::vector< DependentVariableType >( ) )
{
    return createOneDimensionalInterpolator< IndependentVariableType, DependentVariableType >(
            dataToInterpolate,
            interpolatorSettings,
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ),
            firstDerivativesOfDataToIntepolate );
}

}  // namespace interpolators
}  // namespace tudat

#endif  // TUDATPY_CREATE_ONE_DIMENSIONAL_INTERPOLATOR_BASIC_H
