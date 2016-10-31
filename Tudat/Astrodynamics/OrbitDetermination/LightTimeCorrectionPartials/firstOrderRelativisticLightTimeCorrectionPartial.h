#ifndef FIRSTORDERRELATIVISTICLIGHTTIMECORRECTIONPARTIAL_H
#define FIRSTORDERRELATIVISTICLIGHTTIMECORRECTIONPARTIAL_H

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter(
        const double singleBodyLightTimeCorrection, const double bodyGravitationalParameter );

double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtPpnParameterGamma(
        const double totalLightTimeCorrection, const double ppnParameterGamma );

class FirstOrderRelativisticLightTimeCorrectionPartial: public LightTimeCorrectionPartial
{
public:

    FirstOrderRelativisticLightTimeCorrectionPartial(
            const boost::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > correctionCalculator );

    ~FirstOrderRelativisticLightTimeCorrectionPartial( ){ }

    SingleOneWayRangePartialReturnType wrtBodyGravitationalParameter(
            const std::vector< basic_mathematics::Vector6d >& states, const std::vector< double >& times, const int bodyIndex )
    {
        double partialValue = getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter(
                    correctionCalculator_->getCurrentLightTimeCorrectionComponent( bodyIndex ),
                    perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( ) );
        return std::make_pair(
                    ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partialValue ).finished( ), ( times[ 0 ] + times[ 1 ] ) / 2.0 );
    }

    SingleOneWayRangePartialReturnType wrtPpnParameterGamma(
            const std::vector< basic_mathematics::Vector6d >& states, const std::vector< double >& times )
    {
        double partialValue = getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtPpnParameterGamma(
                    correctionCalculator_->getCurrentTotalLightTimeCorrection( ), correctionCalculator_->getPpnParameterGammaFunction_( )( ) );
        return std::make_pair(
                    ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partialValue ).finished( ), ( times[ 0 ] + times[ 1 ] ) / 2.0 );
    }

    std::vector< std::string > getPerturbingBodies( )
    {
        return perturbingBodies_;
    }


protected:

    boost::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > correctionCalculator_;

    std::vector< std::string > perturbingBodies_;

    std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;

};

}

}

#endif // FIRSTORDERRELATIVISTICLIGHTTIMECORRECTIONPARTIAL_H
