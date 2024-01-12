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
        Eigen::Vector3d inertialVectorFromSource = radiationPressureAcceleration_->getTargetPositionWrtSource( );
        Eigen::Vector3d bodyFixedUnitVectorToSource =
            - ( radiationPressureAcceleration_->getTargetRotationFromGlobalToLocalFrame( ) *
                inertialVectorFromSource.normalized( ) );
        double distanceToSource = inertialVectorFromSource.norm( );
        Eigen::Vector3d currentAcceleration = radiationPressureAcceleration_->getAcceleration( );
        double currentRadiationPressure = radiationPressureAcceleration_->getCurrentRadiationPressure();
        double currentMass = radiationPressureAcceleration_->getCurrentTargetMass( );

        currentPartialWrtPosition_.setZero( );

        if( currentRadiationPressure > 0.0  )
        {
            currentSourceUnitVectorPartial_ =  -1.0 / distanceToSource * (
                        Eigen::Matrix3d::Identity( ) - bodyFixedUnitVectorToSource * bodyFixedUnitVectorToSource.transpose( ) );
            currentRadiationPressurePositionPartial_ =
                    2.0 * currentRadiationPressure * bodyFixedUnitVectorToSource.transpose( ) / ( distanceToSource );

            currentCosineAnglePartial_ = Eigen::Matrix< double, 1, 3 >::Zero( );
            Eigen::Vector3d currentPanelReactionVector = Eigen::Vector3d::Zero( );
            Eigen::Vector3d currentPanelNormal = Eigen::Vector3d::Zero( );
            Eigen::Matrix3d currentPanelPartialContribution = Eigen::Matrix3d::Zero( );


            double currentPanelArea = 0.0, currentPanelEmissivity = 0.0, cosineOfPanelInclination = 0.0;

            for( int i = 0; i < panelledTargetModel_->getTotalNumberOfPanels( ); i++ )
            {
                currentPanelNormal = panelledTargetModel_->getSurfaceNormals( ).at( i );
                cosineOfPanelInclination = panelledTargetModel_->getSurfacePanelCosines( ).at( i );

                currentPanelPartialContribution.setZero( );
                if( cosineOfPanelInclination > 0.0 )
                {
                    currentCosineAnglePartial_ = currentPanelNormal.transpose( ) * currentSourceUnitVectorPartial_;

                    currentPanelArea = panelledTargetModel_->getBodyFixedPanels( ).at( i )->getPanelArea( );
                    currentPanelReactionVector = panelledTargetModel_->getPanelForces( ).at( i ) / ( currentRadiationPressure * currentPanelArea );

                    currentPanelPartialContribution += panelledTargetModel_->getFullPanels( ).at( i )->getReflectionLaw( )->
                        evaluateReactionVectorDerivativeWrtTargetPosition(
                            currentPanelNormal, -bodyFixedUnitVectorToSource, cosineOfPanelInclination, currentPanelReactionVector,
                            currentSourceUnitVectorPartial_, currentCosineAnglePartial_ );

                    currentPanelPartialContribution *= currentRadiationPressure * currentPanelArea;


                }
                currentPartialWrtPosition_ += currentPanelPartialContribution;
            }

            Eigen::Matrix3d rotationToInertialFrame =
                Eigen::Matrix3d( radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrame( ) );
            currentPartialWrtPosition_ = rotationToInertialFrame * currentPartialWrtPosition_ * rotationToInertialFrame.transpose( );
            currentPartialWrtPosition_ /= currentMass;

            currentPartialWrtPosition_ += currentAcceleration / currentRadiationPressure *
                ( currentRadiationPressurePositionPartial_ * rotationToInertialFrame.transpose( ) );

        }
        currentTime_ = currentTime;

    }
}


}

}

