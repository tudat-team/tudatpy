#ifndef TUDAT_PANELLED_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_PANELLED_AERODYNAMIC_COEFFICIENT_INTERFACE_H

#include <vector>
#include <map>

#include <memory>
#include <Eigen/Core>

#include "tudat/astro/system_models/selfShadowing.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/aerodynamics/gasSurfaceInteractionModel.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"

namespace tudat
{
namespace aerodynamics
{

class PanelledAerodynamicCoefficientInterface : public AerodynamicCoefficientInterface
{
public:
    PanelledAerodynamicCoefficientInterface( const std::shared_ptr< GasSurfaceInteractionModel > gasSurfaceInteractionModel,
                                             const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames,
                                             const double referenceArea,
                                             const AerodynamicCoefficientFrames coefficientFrame ):
        AerodynamicCoefficientInterface( TUDAT_NAN, referenceArea, Eigen::Vector3d::Zero( ), independentVariableNames, coefficientFrame ),
        gasSurfaceInteractionModel_( gasSurfaceInteractionModel )
    {
        numberOfIndependentVariables_ = independentVariableNames_.size( );
        gasSurfaceInteractionModelType_ = gasSurfaceInteractionModel_->getGasSurfaceInteractionModelType( );
    }

    void updateCurrentCoefficients( const std::vector< double >& independentVariables, const double currentTime = TUDAT_NAN ) override
    {
        gasSurfaceInteractionModel_->setIncomingDirection(
                Eigen::Vector3d( -std::cos( independentVariables[ 0 ] ) * std::cos( independentVariables[ 1 ] ),
                                 -std::sin( independentVariables[ 1 ] ),
                                 -std::sin( independentVariables[ 0 ] ) * std::cos( independentVariables[ 1 ] ) ) );
        if( gasSurfaceInteractionModelType_ == sentman || gasSurfaceInteractionModelType_ == cook )
        {
            gasSurfaceInteractionModel_->setFreeStreamTemperature( independentVariables[ 2 ] );
        }
        if( gasSurfaceInteractionModelType_ == sentman )
        {
            gasSurfaceInteractionModel_->setAirSpeed( independentVariables[ 3 ] );
        }
        currentForceCoefficients_ = gasSurfaceInteractionModel_->computeAerodynamicCoefficients( );

        referenceArea_ = gasSurfaceInteractionModel_->getReferenceArea( );
    }

    GasSurfaceInteractionModelType getGasSurfaceInteractionModelType( ) const
    {
        return gasSurfaceInteractionModelType_;
    }

    std::shared_ptr< GasSurfaceInteractionModel > getGasSurfaceInteractionModel( )
    {
        return gasSurfaceInteractionModel_;
    }

private:
    GasSurfaceInteractionModelType gasSurfaceInteractionModelType_;

    std::shared_ptr< GasSurfaceInteractionModel > gasSurfaceInteractionModel_;
};

//! Function to retrieve the panelled coefficient interface, also when it is wrapped in scaled coefficients
//! (nullptr if the coefficients are not panelled).
inline std::shared_ptr< PanelledAerodynamicCoefficientInterface > getPanelledCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface )
{
    auto scaledCoefficientInterface = std::dynamic_pointer_cast< ScaledAerodynamicCoefficientInterface >( coefficientInterface );
    return std::dynamic_pointer_cast< PanelledAerodynamicCoefficientInterface >(
            scaledCoefficientInterface != nullptr ? scaledCoefficientInterface->getBaseCoefficientInterface( ) : coefficientInterface );
}

}  // namespace aerodynamics

}  // namespace tudat

#endif  // TUDAT_PANELLED_AERODYNAMIC_COEFFICIENT_INTERFACE_H
