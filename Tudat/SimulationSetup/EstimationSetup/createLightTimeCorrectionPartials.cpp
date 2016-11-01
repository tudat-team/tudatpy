#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"


namespace tudat
{

namespace observation_partials
{

std::vector< boost::shared_ptr< LightTimeCorrectionPartial > > createLightTimeCorrectionPartials(
        const std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrectionList )
{
    std::vector< boost::shared_ptr< LightTimeCorrectionPartial > > partialList;

    for( unsigned int i = 0; i < lightTimeCorrectionList.size( ); i++ )
    {
        switch( lightTimeCorrectionList.at( i )->getLightTimeCorrectionType( ) )
        {
        case observation_models::first_order_relativistic:
        {
            boost::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > currentCorrection =
                    boost::dynamic_pointer_cast< observation_models::FirstOrderLightTimeCorrectionCalculator >( lightTimeCorrectionList.at( i ) );
            if( currentCorrection == NULL )
            {
                std::cerr<<"Error when making first order light time correction partial, type id "<<observation_models::first_order_relativistic<<
                           " not consistent with class type."<<std::endl;
            }
            else
            {
                partialList.push_back( boost::make_shared< FirstOrderRelativisticLightTimeCorrectionPartial >( currentCorrection ) );
            }

            break;
        }
        default:
            break;
        }
    }

    return partialList;
}

}

}

