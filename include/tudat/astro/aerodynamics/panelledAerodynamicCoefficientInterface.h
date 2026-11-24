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
        if( gasSurfaceInteractionModelType_ == sentman )
        {
            gasSurfaceInteractionModel_->setFreeStreamTemperature( independentVariables[ 2 ] );
            gasSurfaceInteractionModel_->setAirSpeed( independentVariables[ 3 ] );
            gasSurfaceInteractionModel_->setSpecifiGasConstant( independentVariables[ 4 ] );
        }
        else if( gasSurfaceInteractionModelType_ == cook )
        {
            gasSurfaceInteractionModel_->setAirSpeed( independentVariables[ 2 ] );
            gasSurfaceInteractionModel_->setSpecifiGasConstant( independentVariables[ 3 ] );
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

}  // namespace aerodynamics

}  // namespace tudat

#endif  // TUDAT_PANELLED_AERODYNAMIC_COEFFICIENT_INTERFACE_H
