#ifndef CREATELIGHTTIMECORRECTIONPARTIALS_H
#define CREATELIGHTTIMECORRECTIONPARTIALS_H

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

std::vector< boost::shared_ptr< LightTimeCorrectionPartial > > createLightTimeCorrectionPartials(
        const std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrectionList );

}

}

#endif // CREATELIGHTTIMECORRECTIONPARTIALS_H
