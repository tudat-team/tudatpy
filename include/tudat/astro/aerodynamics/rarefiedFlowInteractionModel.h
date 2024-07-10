/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_RAREFIEDFLOWMODEL_H
#define TUDAT_RAREFIEDFLOWMODEL_H

#include <Eigen/Core>

namespace tudat
{
namespace aerodynamics
{

class RarefiedFlowInteractionModel{
public:
    RarefiedFlowInteractionModel() = default;

    virtual ~RarefiedFlowInteractionModel() = default;

    Eigen::Vector3d computePanelForceCoefficientVector( 
        double CosineOfNormalDragAngle, //gammai in Doornbos
        double CosineOfNormalLiftAngle, //li in Doornbos
        double panelSurfaceArea,
        double panelTemperature,
        Eigen::Vector3d liftUnitVector,
        Eigen::Vector3d dragUnitVecotr,
        double Vinf,
        double T_atm,
        std::vector<double> number_densities,
        double total_number_density,
        double Aref);

    Eigen::Vector3d computePanelMomentCoefficientVector( 
        Eigen::Vector3d panelForceCoefficientVector,
        Eigen::Vector3d panelPositionVector,
        double lref
        );

    double get_Cd_ij(double Vinf, double Tinf, double mj, double gammai, double Ai, double panelTemperature, double Aref);

    double get_Cl_ij(double Vinf, double Tinf, double mj, double gammai, double li, double Ai, double panelTemperature, double Aref);

    double get_Sj(double Vinf, double Tinf, double mj);

    double get_Gj(double Sj);

    double get_Pij(double Sj, double gammai);

    double get_Zij(double Sj, double gammai);

    double get_alpha(double Tinf);

    double get_reboundVelocityFreestreamVelocityRatio(double alpha, double Ti, double Vinf);

};
} // tudat
} // aerodynamic

#endif //TUDAT_RAREFIEDFLOWMODEL_H
