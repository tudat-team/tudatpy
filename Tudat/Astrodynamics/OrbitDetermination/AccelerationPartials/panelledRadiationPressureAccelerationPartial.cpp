#include <iostream>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/panelledRadiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for updating partial w.r.t. the bodies' positions.
void PanelledRadiationPressurePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        Eigen::Vector3d unitVectorToSource = radiationPressureAcceleration_->getCurrentVectorToSource( );
        double distanceToSource = radiationPressureAcceleration_->getCurrentDistanceToSource( );
        Eigen::Vector3d currentAcceleration = radiationPressureAcceleration_->getAcceleration( );
        double currentRadiationPressure = radiationPressureAcceleration_->getCurrentRadiationPressure( );
        double currentMass = radiationPressureAcceleration_->getCurrentMass( );

        currentPartialWrtPosition_.setZero( );

        if( currentRadiationPressure > 0.0 && currentAcceleration.norm( ) > 0.0 )
        {
            Eigen::Matrix3d currentSourceUnitVectorPartial =  -1.0 / distanceToSource * (
                        Eigen::Matrix3d::Identity( ) - unitVectorToSource * unitVectorToSource.transpose( ) );
            Eigen::Matrix< double, 1, 3 > currentRadiationPressurePositionPartial =
                    2.0 * currentRadiationPressure * unitVectorToSource.transpose( ) / ( distanceToSource );

            Eigen::Matrix< double, 1, 3 > currentCosineAnglePartial = Eigen::Matrix< double, 1, 3 >::Zero( );
            Eigen::Vector3d currentPanelAcceleration = Eigen::Vector3d::Zero( );
            Eigen::Vector3d currentPanelNormal = Eigen::Vector3d::Zero( );


            double currentPanelArea = 0.0, currentPanelEmissivity = 0.0, cosineOfPanelInclination = 0.0;

            for( int i = 0; i < radiationPressureInterface_->getNumberOfPanels( ); i++ )
            {
                currentPanelNormal = radiationPressureInterface_->getCurrentSurfaceNormal( i );
                cosineOfPanelInclination = currentPanelNormal.dot( unitVectorToSource );

                if( cosineOfPanelInclination > 0.0 )
                {
                    currentCosineAnglePartial = currentPanelNormal.transpose( ) * currentSourceUnitVectorPartial;

                    currentPanelAcceleration = radiationPressureAcceleration_->getCurrentPanelAcceleration( i );
                    currentPanelArea = radiationPressureInterface_->getArea( i );
                    currentPanelEmissivity = radiationPressureInterface_->getEmissivity( i );

                    currentPartialWrtPosition_ +=
                            currentPanelAcceleration * currentCosineAnglePartial / cosineOfPanelInclination;
                    currentPartialWrtPosition_ -=
                            currentRadiationPressure / currentMass * currentPanelArea * cosineOfPanelInclination * (
                                ( 1.0 - currentPanelEmissivity ) * currentSourceUnitVectorPartial +
                                2.0 * currentPanelEmissivity * currentPanelNormal * currentCosineAnglePartial );

                }
            }

            currentPartialWrtPosition_ += currentAcceleration / currentRadiationPressure * currentRadiationPressurePositionPartial;
        }
        currentTime_ = currentTime;

    }
}


}

}

