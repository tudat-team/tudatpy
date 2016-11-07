#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter(
        const double singleBodyLightTimeCorrection, const double bodyGravitationalParameter )
{
    return singleBodyLightTimeCorrection / bodyGravitationalParameter;
}

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtPpnParameterGamma(
        const double totalLightTimeCorrection, const double ppnParameterGamma )
{
    return totalLightTimeCorrection / ( ppnParameterGamma + 1.0 );
}


FirstOrderRelativisticLightTimeCorrectionPartial::FirstOrderRelativisticLightTimeCorrectionPartial(
        const boost::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > correctionCalculator ):
    LightTimeCorrectionPartial( observation_models::first_order_relativistic ),
    correctionCalculator_( correctionCalculator )
{
    perturbingBodies_ = correctionCalculator_->getPerturbingBodyNames( );
    perturbingBodyGravitationalParameterFunctions_ = correctionCalculator_->getPerturbingBodyGravitationalParameterFunctions( );

}

}

}
