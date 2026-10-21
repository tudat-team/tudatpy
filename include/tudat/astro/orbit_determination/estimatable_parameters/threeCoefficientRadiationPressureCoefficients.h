/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_THREECOEFFICIENTRADIATIONPRESSURECOEFFICIENTS_H
#define TUDAT_THREECOEFFICIENTRADIATIONPRESSURECOEFFICIENTS_H

#include "tudat/astro/electromagnetism/threeCoefficientRadiationPressureAcceleration.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{
namespace estimatable_parameters
{

//! Three constant area coefficients of one or more three-coefficient radiation-pressure models.
class ThreeCoefficientRadiationPressureCoefficients : public EstimatableParameter< Eigen::VectorXd >
{
public:
    ThreeCoefficientRadiationPressureCoefficients(
            const std::vector< std::shared_ptr< electromagnetism::ThreeCoefficientRadiationPressureAcceleration > >& accelerationModels,
            const std::string& acceleratedBody,
            const std::string& sourceBody ):
        EstimatableParameter< Eigen::VectorXd >( three_coefficient_radiation_pressure_coefficients, acceleratedBody, sourceBody ),
        accelerationModels_( accelerationModels )
    {
        if( accelerationModels_.empty( ) )
        {
            throw std::runtime_error(
                    "Error when creating three-coefficient radiation-pressure parameter: no acceleration models were provided." );
        }
    }

    Eigen::VectorXd getParameterValue( ) override
    {
        const Eigen::Vector3d referenceCoefficients = accelerationModels_.front( )->getCoefficients( );
        for( const auto& accelerationModel : accelerationModels_ )
        {
            if( accelerationModel->getCoefficients( ) != referenceCoefficients )
            {
                throw std::runtime_error(
                        "Error when retrieving three-coefficient radiation-pressure parameter: model values are inconsistent." );
            }
        }
        return referenceCoefficients;
    }

    void setParameterValue( const Eigen::VectorXd parameterValue ) override
    {
        if( parameterValue.size( ) != 3 )
        {
            throw std::runtime_error(
                    "Error when resetting three-coefficient radiation-pressure parameter: parameter must have size three." );
        }
        for( const auto& accelerationModel : accelerationModels_ )
        {
            accelerationModel->resetCoefficients( parameterValue );
        }
    }

    int getParameterSize( ) override
    {
        return 3;
    }

    std::string getParameterDescription( ) override
    {
        return "three-coefficient radiation-pressure areas A1, A2, and A3";
    }

private:
    std::vector< std::shared_ptr< electromagnetism::ThreeCoefficientRadiationPressureAcceleration > > accelerationModels_;
};

}  // namespace estimatable_parameters
}  // namespace tudat

#endif  // TUDAT_THREECOEFFICIENTRADIATIONPRESSURECOEFFICIENTS_H
