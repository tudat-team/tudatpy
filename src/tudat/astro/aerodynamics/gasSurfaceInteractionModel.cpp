/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
#include <functional>
#include <memory>
#include <Eigen/Core>
#include "tudat/astro/system_models/selfShadowing.h"
#include "tudat/astro/aerodynamics/gasSurfaceInteractionModel.h"

namespace tudat
{
namespace aerodynamics
{

void GasSurfaceInteractionModel::updateMembers( )
{
    if ( maximumNumberOfPixels_ == 0 )
    {   
        // SSH off
        illuminatedPanelFractions_ = unityIlluminationFraction_;
    }
    else
    {
        // SSH on
        aerodynamicSelfShadowing_->reset( );
        aerodynamicSelfShadowing_->updateIlluminatedPanelFractions( incomingDirection_ );
        illuminatedPanelFractions_ = aerodynamicSelfShadowing_->getIlluminatedPanelFractions( );
    }
}

Eigen::Vector3d ConstantInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );

    double cosineDelta, sineDelta;
    double actualCrossSectionArea = 0;
    Eigen::Vector3d panelNormal;

    for ( int i=0; i<totalNumberOfPanels_; i++ )
    {   
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if ( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if ( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        actualCrossSectionArea += allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
    }
    referenceArea_ = actualCrossSectionArea;
    return constantAerodynamicCoefficients_;
    
}

Eigen::Vector3d NewtonGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );

    double cosineDelta, sineDelta;
    double Ct = 0;
    double Cp, panelArea;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;

    for ( int i=0; i<totalNumberOfPanels_; i++ )
    {   
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if ( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if ( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        // Cp
        Cp = 2 * cosineDelta * cosineDelta;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface 
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame += -Cp * panelNormal * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if ( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
    
}

Eigen::Vector3d StorchGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );
    double cosineDelta, sineDelta;
    double Cp, Ct, panelArea;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;
    for ( int i=0; i<totalNumberOfPanels_; i++ )
    {   
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if ( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if ( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        sineDelta = std::sqrt(std::max(0.0, 1 - cosineDelta * cosineDelta));
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        //Cp 
        Cp = 2 * cosineDelta * ( 
            allPanels_[ i ]->getNormalAccomodationCoefficient( ) * allPanels_[ i ]->getNormalVelocityAtWallRatio( ) + 
            ( 2 - allPanels_[ i ]->getNormalAccomodationCoefficient( ) ) * cosineDelta );
        // Ct
        Ct = 2 * allPanels_[ i ]->getTangentialAccomodationCoefficient( ) * sineDelta * cosineDelta;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface 
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame += ( -Cp * panelNormal - Ct * (
             incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if ( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

Eigen::Vector3d SentmanGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );
    speedRatio_ = airSpeed_ / std::sqrt( 2 * specificGasConstant_ * freeStreamTemperature_ );
    double cosineDelta, sineDelta;
    double Cp, Ct, panelArea;
    double erf, exp;
    double sqrtPi = std::sqrt( mathematical_constants::PI );
    incidentTemperature_ = 2.0/3.0 * speedRatio_ * speedRatio_ * freeStreamTemperature_;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;

    for ( int i=0; i<totalNumberOfPanels_; i++ )
    {   
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if ( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if ( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        erf = std::erf( speedRatio_ * cosineDelta );
        exp = std::exp( -speedRatio_ * speedRatio_ * cosineDelta * cosineDelta );
        sineDelta = std::sqrt(std::max(0.0, 1 - cosineDelta * cosineDelta));
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        //Cp 
        Cp = ( cosineDelta * cosineDelta ) * ( 1 + erf ) + 
            cosineDelta / ( speedRatio_ * sqrtPi ) * exp +
            0.5 * std::sqrt( 2.0/3.0 * ( 1 + (allPanels_[ i ]->getEnergyAccomodationCoefficient( ) * 
            allPanels_[ i ]->getPanelTemperature( ) ) / ( incidentTemperature_ - 1) ) ) * ( sqrtPi * cosineDelta * ( 1 + erf ) +
            1.0 / speedRatio_ * exp );
        //Ct
        Ct = sineDelta * cosineDelta * ( 1 + erf ) + sineDelta / ( speedRatio_ * sqrtPi ) * exp;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface 
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame += ( -Cp * panelNormal - Ct * (
             incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if ( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

Eigen::Vector3d CookGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );
    double cosineDelta, sineDelta;
    double Cd, Cl, Cp, Ct, panelArea;
    double sqrt;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;

    for ( int i=0; i<totalNumberOfPanels_; i++ )
    {   
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if ( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if ( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        sineDelta = std::sqrt(std::max(0.0, 1 - cosineDelta * cosineDelta));
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        sqrt = std::sqrt(  1 + ( allPanels_[ i ]->getEnergyAccomodationCoefficient( ) * 
            allPanels_[ i ]->getPanelTemperature( ) ) / ( freeStreamTemperature_ - 1) );
        //Cd
        Cd = 2 * cosineDelta * ( 1 + 2.0/3.0 * cosineDelta * sqrt );
        // Cl
        Cl = 4.0/3.0 * sineDelta * cosineDelta * sqrt;
        // convert cd, cd to cp, ct
        Cp = cosineDelta * Cd + sineDelta * Cl;
        Ct = sineDelta * Cd - cosineDelta * Cl;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface 
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame += ( -Cp * panelNormal - Ct * (
             incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if ( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

}

}
