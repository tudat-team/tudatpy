#ifndef LIGHTTIMECORRECTIONPARTIAL_H
#define LIGHTTIMECORRECTIONPARTIAL_H

#include <boost/smart_ptr/enable_shared_from_this.hpp>

#include <vector>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

class LightTimeCorrectionPartial
{
public:
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    LightTimeCorrectionPartial( const observation_models::LightTimeCorrectionType correctionType ):
        correctionType_( correctionType ){ }

    virtual ~LightTimeCorrectionPartial( ){ }

    observation_models::LightTimeCorrectionType getCorrectionType( )
    {
        return correctionType_;
    }

protected:

    observation_models::LightTimeCorrectionType correctionType_;

};

std::pair< boost::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
        const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool > getLightTimeParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
        const boost::shared_ptr< LightTimeCorrectionPartial > lightTimeCorrectionPartial );

}

}

#endif // LIGHTTIMECORRECTIONPARTIAL_H
