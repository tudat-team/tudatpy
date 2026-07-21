/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITYDEFORMATIONSETTINGS_H
#define TUDAT_GRAVITYDEFORMATIONSETTINGS_H

#include <functional>
#include <memory>
#include <map>
#include <string>
#include "tudat/astro/basic_astro/gravityDeformationModelTypes.h"

namespace tudat
{

namespace simulation_setup
{

// Class for providing settings for gravity deformation models.
/*
 *  Class for providing settings for gravity deformation models.
 */
class GravityDeformationSettings
{
public:
    // Constructor, sets type of deformation.
    /*
     *  Constructor, sets type of deformation.
     *  \param deformationType Type of acceleration from GravityDeformationType enum.
     */
    GravityDeformationSettings( const basic_astrodynamics::GravityDeformationType deformationType ): deformationType_( deformationType ) {}

    // Destructor.
    virtual ~GravityDeformationSettings( ) {}

    // Type of acceleration from AvailableAcceleration enum.
    basic_astrodynamics::GravityDeformationType deformationType_;
};

class MaxwellDeformationSettings : public GravityDeformationSettings
{
public:
    // Constructor, sets type of acceleration.
    /*
     *  Constructor, sets type of acceleration.
     *  \param accelerationType Type of acceleration from AvailableAcceleration enum.
     */
    MaxwellDeformationSettings( const double maxwellRelaxationTime,
                                const double globalRelaxationTime,
                                const double loveNumber,
                                const int maximumDegree,
                                const int maximumOrder,
                                const std::vector< std::string > perturbingBody,
                                const Eigen::VectorXd staticCoefficients = Eigen::VectorXd::Zero( 5 ),
                                const bool includeOrder1 = true,
                                const bool includeCentrifugalPotential = false ):
        GravityDeformationSettings( basic_astrodynamics::maxwell_deformation ), maxwellRelaxationTime_( maxwellRelaxationTime ),
        globalRelaxationTime_( globalRelaxationTime ), loveNumber_( loveNumber ), maximumDegree_( maximumDegree ),
        maximumOrder_( maximumOrder ), perturbingBody_( perturbingBody ), staticCoefficients_( staticCoefficients ),
        includeOrder1_( includeOrder1 ), includeCentrifugalPotential_( includeCentrifugalPotential )
    {}

    // Destructor.
    virtual ~MaxwellDeformationSettings( ) {}

    const double maxwellRelaxationTime_;
    const double globalRelaxationTime_;
    const double loveNumber_;
    const int maximumDegree_;
    const int maximumOrder_;
    const std::vector< std::string > perturbingBody_;
    Eigen::VectorXd staticCoefficients_;
    const bool includeOrder1_;
    const bool includeCentrifugalPotential_;
};

typedef std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > SelectedGravityDeformationModelMap;

inline std::shared_ptr< GravityDeformationSettings > maxwellDeformationSettings(
        const double maxwellRelaxationTime,
        const double globalRelaxationTime,
        const double loveNumber,
        const int maximumDegree,
        const int maximumOrder,
        const std::string perturbingBody,
        const Eigen::VectorXd staticCoefficients = Eigen::VectorXd::Zero( 5 ),
        const bool includeOrder1 = true,
        const bool includeCentrifugalPotential = false )
{
    std::vector< std::string > perturbingBodies = { perturbingBody };
    return std::make_shared< MaxwellDeformationSettings >( maxwellRelaxationTime,
                                                           globalRelaxationTime,
                                                           loveNumber,
                                                           maximumDegree,
                                                           maximumOrder,
                                                           perturbingBodies,
                                                           staticCoefficients,
                                                           includeOrder1,
                                                           includeCentrifugalPotential );
}

inline std::shared_ptr< GravityDeformationSettings > maxwellDeformationSettings(
        const double maxwellRelaxationTime,
        const double globalRelaxationTime,
        const double loveNumber,
        const int maximumDegree,
        const int maximumOrder,
        const std::vector< std::string > perturbingBodies,
        const Eigen::VectorXd staticCoefficients = Eigen::VectorXd::Zero( 5 ),
        const bool includeOrder1 = true,
        const bool includeCentrifugalPotential = false )
{
    return std::make_shared< MaxwellDeformationSettings >( maxwellRelaxationTime,
                                                           globalRelaxationTime,
                                                           loveNumber,
                                                           maximumDegree,
                                                           maximumOrder,
                                                           perturbingBodies,
                                                           staticCoefficients,
                                                           includeOrder1,
                                                           includeCentrifugalPotential );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_GRAVITYDEFORMATIONSETTINGS_H
