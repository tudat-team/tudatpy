#include <iostream>

#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/acceleration_partials/panelledRadiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for updating partial w.r.t. the bodies' positions.
void PanelledRadiationPressurePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        Eigen::Vector3d unitVectorToSource =
            - ( radiationPressureAcceleration_->getTargetRotationFromGlobalToLocalFrame( ) *
                radiationPressureAcceleration_->getTargetPositionWrtSource( ).normalized( ) );
        double distanceToSource = unitVectorToSource.norm( );
        Eigen::Vector3d currentAcceleration = radiationPressureAcceleration_->getAcceleration( );
        double currentRadiationPressure = panelledTargetModel_->getRadiationPressure( );
        double currentMass = radiationPressureAcceleration_->getCurrentTargetMass( );

        currentPartialWrtPosition_.setZero( );

        if( currentRadiationPressure > 0.0  )
        {
            Eigen::Matrix3d currentSourceUnitVectorPartial =  -1.0 / distanceToSource * (
                        Eigen::Matrix3d::Identity( ) - unitVectorToSource * unitVectorToSource.transpose( ) );
            std::cout<<"Unit vector partial "<<currentSourceUnitVectorPartial<<std::endl;
            Eigen::Matrix< double, 1, 3 > currentRadiationPressurePositionPartial =
                    2.0 * currentRadiationPressure * unitVectorToSource.transpose( ) / ( distanceToSource );

            Eigen::Matrix< double, 1, 3 > currentCosineAnglePartial = Eigen::Matrix< double, 1, 3 >::Zero( );
            Eigen::Vector3d currentPanelReactionVector = Eigen::Vector3d::Zero( );
            Eigen::Vector3d currentPanelNormal = Eigen::Vector3d::Zero( );
            Eigen::Matrix3d currentPanelPartialContribution = Eigen::Matrix3d::Zero( );


            double currentPanelArea = 0.0, currentPanelEmissivity = 0.0, cosineOfPanelInclination = 0.0;

            for( int i = 0; i < panelledTargetModel_->getTotalNumberOfPanels( ); i++ )
            {
                currentPanelNormal = panelledTargetModel_->getSurfaceNormals( ).at( i );
                cosineOfPanelInclination = currentPanelNormal.dot( unitVectorToSource.normalized( ) );

                currentPanelPartialContribution.setZero( );
                if( cosineOfPanelInclination > 0.0 )
                {
                    currentCosineAnglePartial = currentPanelNormal.transpose( ) * currentSourceUnitVectorPartial;

                    currentPanelArea = panelledTargetModel_->getBodyFixedPanels( ).at( i )->getPanelArea( );
                    currentPanelReactionVector = panelledTargetModel_->getPanelForces( ).at( i ) / ( currentRadiationPressure * currentPanelArea );

                    currentPanelPartialContribution =
                        currentPanelReactionVector * currentCosineAnglePartial / cosineOfPanelInclination;
                    std::cout<<"Panel contribution "<<currentPanelPartialContribution<<std::endl;
                    currentPanelPartialContribution +=  cosineOfPanelInclination * panelledTargetModel_->getFullPanels( ).at( i )->getReflectionLaw( )->
                        evaluateReactionVectorDerivativeWrtTargetPosition(
                            currentPanelNormal, -unitVectorToSource, cosineOfPanelInclination, currentPanelReactionVector,
                            currentSourceUnitVectorPartial, currentCosineAnglePartial );
                    std::cout<<"Panel contribution "<<currentPanelPartialContribution<<std::endl;

                    currentPanelPartialContribution *= currentRadiationPressure * currentPanelArea;
                    std::cout<<"Panel contribution "<<currentPanelPartialContribution<<std::endl;


                }
                currentPartialWrtPosition_ += currentPanelPartialContribution;
            }

            Eigen::Matrix3d rotationToInertialFrame =
                Eigen::Matrix3d( radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrame( ) );
            currentPartialWrtPosition_ = rotationToInertialFrame * currentPartialWrtPosition_ * rotationToInertialFrame.transpose( );
            currentPartialWrtPosition_ /= currentMass;
            std::cout<<"Total partial: "<<currentPartialWrtPosition_<<std::endl;

            currentPartialWrtPosition_ += currentAcceleration / currentRadiationPressure * currentRadiationPressurePositionPartial;
            std::cout<<"Total partial: "<<currentPartialWrtPosition_<<std::endl;

        }
        currentTime_ = currentTime;

    }
}


}

}

