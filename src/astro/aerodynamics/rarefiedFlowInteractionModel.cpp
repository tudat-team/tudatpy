/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h"

#include <Eigen/Core>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace aerodynamics
{

// RarefiedFlowInteractionModel::RarefiedFlowInteractionModel() = default;

/*! Computes aerodynamic force coefficients for a single panel
*  \param cosineOfNormalDragAngle Cosine of the angle between the freestream velocity and the panel normal
*  \param cosineOfNormalLiftAngle Cosine of the angle between the freestream velocity and the panel normal
*  \param panelSurfaceArea Area of the panel
*  \param panelTemperature Temperature of the panel
*  \param liftUnitVector Lift vector
*  \param dragUnitVecotr Drag vector
*  \param freestreamVelocity Freestream velocity
*  \param atmosphericTemperature Atmospheric temperature
*  \param numberDensities Number densities of species
*  \param totalNumberDensity Total number density
*  \param referenceArea Reference area
*  \return Force coefficient vector
*/
Eigen::Vector3d RarefiedFlowInteractionModel::computePanelForceCoefficientVector( 
    double cosineOfNormalDragAngle, //cosineOfNormalDragAngle in Doornbos
    double cosineOfNormalLiftAngle, //cosineOfNormalLiftAngle in Doornbos
    double panelSurfaceArea,
    double panelTemperature,
    Eigen::Vector3d liftUnitVector,
    Eigen::Vector3d dragUnitVecotr,
    double freestreamVelocity,
    double atmosphericTemperature,
    std::vector<double> numberDensities,
    double totalNumberDensity,
    double referenceArea)
{
    
    // To be moved to correct location (mathematicalConstants.h maybe?)
    std::vector<double> nrlmsise00SpeciesAtomicMasses{4.002602, 15.999, 28.0134, 31.9988, 39.948, 1.008, 14.0067, 15.999}; // [g/mol] atomic mass of species in NRLMSISE-00 model

    // Initialize force coefficient vector
    Eigen::Vector3d panelForceCoefficientVector = Eigen::Vector3d::Zero();

    if (cosineOfNormalDragAngle > 0){ // check if panel is pointing into the flow
        for (int j_species = 0; j_species < nrlmsise00SpeciesAtomicMasses.size(); j_species++) {
            double Cdij = get_Cd_ij(
                freestreamVelocity, atmosphericTemperature, nrlmsise00SpeciesAtomicMasses[j_species], 
                cosineOfNormalDragAngle, panelSurfaceArea, panelTemperature, referenceArea) * numberDensities[j_species] / totalNumberDensity;
            
            double Clij = get_Cl_ij(
                freestreamVelocity, atmosphericTemperature, nrlmsise00SpeciesAtomicMasses[j_species], 
                cosineOfNormalDragAngle, cosineOfNormalLiftAngle, panelSurfaceArea, panelTemperature, referenceArea) * numberDensities[j_species] / totalNumberDensity;

            panelForceCoefficientVector += Cdij * dragUnitVecotr + Clij * liftUnitVector;
        }
    }
    
    return panelForceCoefficientVector;
}

Eigen::Vector3d RarefiedFlowInteractionModel::computePanelMomentCoefficientVector( 
    Eigen::Vector3d panelForceCoefficientVector,
    Eigen::Vector3d panelPositionVector,
    double referenceLength)
{
    return panelPositionVector.cross(panelForceCoefficientVector) / referenceLength;
}

// private:

double RarefiedFlowInteractionModel::get_Cd_ij(double freestreamVelocity, double atmosphericTemperature, double speciesAtomicMass, double cosineOfNormalDragAngle, double panelSurfaceArea, double panelTemperature, double referenceArea)
    {
    // Function to calculate the drag coefficient for a single species j on a single panel i
    // freestreamVelocity: freestream velocity
    // atmosphericTemperature: freestream temperature
    // speciesAtomicMass: molecular mass of species j
    // cosineOfNormalDragAngle: cosine of the angle between the freestream velocity and the panel normal
    // panelSurfaceArea: area of panel i
    // referenceArea: reference area
    // rhoO: atomic oxygen density

    double Sj = get_Sj(freestreamVelocity, atmosphericTemperature, speciesAtomicMass);
    double Pij = get_Pij(Sj, cosineOfNormalDragAngle);
    double Gj = get_Gj(Sj);
    double Qj = 1.0 + Gj;
    double Zij = get_Zij(Sj, cosineOfNormalDragAngle);
    double alpha = get_alpha(atmosphericTemperature);
    double reboundVelocityFreestreamVelocityRatio = get_reboundVelocityFreestreamVelocityRatio(alpha, panelTemperature, freestreamVelocity);

    double Cd_ij = Pij / pow(mathematical_constants::PI, 0.5);
    Cd_ij += cosineOfNormalDragAngle * Qj * Zij;
    Cd_ij += cosineOfNormalDragAngle / 2.0 * reboundVelocityFreestreamVelocityRatio * (cosineOfNormalDragAngle * pow(mathematical_constants::PI, 0.5) * Zij + Pij);
    Cd_ij *= panelSurfaceArea / referenceArea;

    return Cd_ij;
}

double RarefiedFlowInteractionModel::get_Cl_ij(double freestreamVelocity, double atmosphericTemperature, double speciesAtomicMass, double cosineOfNormalDragAngle, double cosineOfNormalLiftAngle, double panelSurfaceArea, double panelTemperature, double referenceArea){
    // Function to calculate the drag coefficient for a single species j on a single panel i
    // freestreamVelocity: freestream velocity
    // atmosphericTemperature: freestream temperature
    // speciesAtomicMass: molecular mass of species j
    // cosineOfNormalDragAngle: cosine of the angle between the freestream velocity and the panel normal
    // cosineOfNormalLiftAngle: sine of the angle between the freestream velocity and the panel normal
    // panelSurfaceArea: area of panel i
    // referenceArea: reference area
    // rhoO: atomic oxygen density

    double Sj = get_Sj(freestreamVelocity, atmosphericTemperature, speciesAtomicMass);
    double Pij = get_Pij(Sj, cosineOfNormalDragAngle);
    double Gj = get_Gj(Sj);
    double Zij = get_Zij(Sj, cosineOfNormalDragAngle);
    double alpha = get_alpha(atmosphericTemperature);
    double reboundVelocityFreestreamVelocityRatio = get_reboundVelocityFreestreamVelocityRatio(alpha, panelTemperature, freestreamVelocity);

    double Cl_ij = cosineOfNormalLiftAngle * Gj * Zij;
    Cl_ij += cosineOfNormalLiftAngle / 2.0 * reboundVelocityFreestreamVelocityRatio * (cosineOfNormalDragAngle * pow(mathematical_constants::PI, 0.5) * Zij + Pij);
    Cl_ij *= panelSurfaceArea / referenceArea;

    return Cl_ij;
}

double RarefiedFlowInteractionModel::get_Sj(double freestreamVelocity, double atmosphericTemperature, double speciesAtomicMass){

    // To be moved to correct location (mathematicalConstants.h maybe?)
    double k = 1.38064852e-23; // [m^2 kg s^-2 K^-1] Boltzmann constant

    return freestreamVelocity / pow(2.0 * k * atmosphericTemperature/speciesAtomicMass, 0.5);
}

double RarefiedFlowInteractionModel::get_Gj(double Sj){
    return 1.0 / (2.0*pow(Sj, 2.0));;
}

double RarefiedFlowInteractionModel::get_Pij(double Sj, double cosineOfNormalDragAngle){
    return (1.0/Sj) * exp(-pow(cosineOfNormalDragAngle, 2.0) * pow(Sj, 2.0));
}

double RarefiedFlowInteractionModel::get_Zij(double Sj, double cosineOfNormalDragAngle){
    return 1.0 + erf(cosineOfNormalDragAngle * Sj);
}

double RarefiedFlowInteractionModel::get_alpha(double atmosphericTemperature){
    // Function to calculate the energy accomodation coefficient, Miyata et. al. 2018
    // rhoO: atomic oxygen density
    // atmosphericTemperature: freestream temperature
    // return (7.5e-17 * rhoO * atmosphericTemperature) / (1 + 7.5e-17 * rhoO * atmosphericTemperature);
    return 1.0;
}

double RarefiedFlowInteractionModel::get_reboundVelocityFreestreamVelocityRatio(double alpha, double panelTemperature, double freestreamVelocity){
    // Function to calculate the ratio of rebound velocity to freestream velocity
    // alpha: energy accomodation coefficient
    // panelTemperature: temperature of panel i
    // freestreamVelocity: freestream velocity

    // To be moved to correct location (mathematicalConstants.h maybe?)
    double R = 8.314462618; // [J K^-1 mol^-1] ideal gas constant

    return pow(0.5 * ( 1.0 + alpha * ( (4.0*R*panelTemperature)/pow(freestreamVelocity, 2.0) - 1.0 ) ), 0.5);
}


} // tudat
} // electromagnetism
