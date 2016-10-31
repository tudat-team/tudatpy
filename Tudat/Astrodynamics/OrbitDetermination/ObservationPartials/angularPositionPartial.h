#ifndef ANGULARPOSITIONPARTIAL_H
#define ANGULARPOSITIONPARTIAL_H

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver );

Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

Eigen::Matrix< double, 2, 3 > calculatePartialOfAngularPositionWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

class AngularPositionScaling: public PositionPartialScaling
{
public:
    void update( const std::vector< basic_mathematics::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd );

    Eigen::Matrix< double, 2, 3 > getTransmitterScalingFactor_( )
    {
        return -scalingFactor_;
    }
    Eigen::Matrix< double, 2, 3 > getReceiverScalingFactor_( )
    {
        return scalingFactor_;
    }

    Eigen::Matrix< double, 2, 3 > getScalingFactor( const observation_models::LinkEndType linkEndType );

private:
    Eigen::Matrix< double, 2, 3 > scalingFactor_;
};

class AngularPositionPartial: public ObservationPartial< 2 >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > AngularPositionPartialReturnType;

    AngularPositionPartial( const boost::shared_ptr< AngularPositionScaling > angularPositionScaler,
                            const std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >& positionPartialList,
                            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        ObservationPartial< 2 >( parameterIdentifier ),
        angularPositionScaler_( angularPositionScaler ), positionPartialList_( positionPartialList ){ }

    ~AngularPositionPartial( ){ }

    AngularPositionPartialReturnType calculatePartial(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime );

protected:
    boost::shared_ptr< AngularPositionScaling > angularPositionScaler_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartialList_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >::iterator positionPartialIterator_;
};

}

}
#endif // ANGULARPOSITIONPARTIAL_H
