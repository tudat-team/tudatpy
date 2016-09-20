#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter(
        const double singleBodyLightTimeCorrection, const double bodyGravitationalParameter )
{
    return singleBodyLightTimeCorrection / bodyGravitationalParameter;
}

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
