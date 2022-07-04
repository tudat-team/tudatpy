/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#ifndef TUDAT_BASICTIDALBODYDEFORMATION_H
#define TUDAT_BASICTIDALBODYDEFORMATION_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>
#include <memory>

#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"

namespace tudat
{

namespace basic_astrodynamics
{



//! Calculate displacement due to degree 2 tide, band-independet love and shida numbers.
Eigen::Vector3d calculateDegreeTwoBasicTidalDisplacement( const double gravitationalParameterRatio,
                                                          const Eigen::Vector3d& stationPositionUnitVector,
                                                          const Eigen::Vector3d& relativeBodyState,
                                                          const double bodyEquatorialRadius,
                                                          const double degreeTwoLoveNumber,
                                                          const double degreeTwoShidaNumber );


//! Calculate displacement due to degree 3 tide, band-independet love and shida numbers.
Eigen::Vector3d calculateDegreeThreeBasicTidalDisplacement( const double gravitationalParameterRatio,
                                                            const Eigen::Vector3d& stationPositionUnitVector,
                                                            const Eigen::Vector3d& relativeBodyState,
                                                            const double bodyEquatorialRadius,
                                                            const double degreeThreeLoveNumber,
                                                            const double degreeThreeShidaNumber );


class BasicTidalBodyDeformation: public BodyDeformationModel
{
public:
    BasicTidalBodyDeformation( const std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction,
                               const std::vector< std::function< Eigen::Vector6d( const double ) > >& deformingBodyStateFunctions,
                               const std::function< Eigen::Quaterniond( const double ) > deformedBodyRotationFunction,
                               const std::function< double( ) > gravitionalParameterOfDeformedBody,
                               const std::vector< std::function< double( ) > >& gravitionalParametersOfDeformingBodies,
                               const double deformedBodyEquatorialRadius,
                               const std::map< int, std::pair< double, double > >& displacementLoveNumbers );

    virtual Eigen::Vector3d calculateDisplacement(
            const double time,
            const Eigen::Vector3d& bodyFixedPosition );

    virtual Eigen::Vector3d calculateDisplacement(
            const double time,
            const std::shared_ptr< ground_stations::GroundStationState > stationState )
    {
        return calculateDisplacement( time, stationState->getNominalCartesianPosition( ) );
    }

protected:

    void updateBodyProperties( const double time );


    Eigen::Vector3d calculateBasicTicalDisplacement(
            const double time,
            const Eigen::Vector3d& bodyFixedPosition,
            const std::map< int, std::pair< double, double > >& loveNumbers );

    std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction_;

    std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyStateFunctions_;

    std::function< Eigen::Quaterniond( const double ) > deformedBodyRotationFunction_;

    std::function< double( ) > gravitionalParameterOfDeformedBody_;

    std::vector< std::function< double( ) > > gravitionalParametersOfDeformingBodies_;

    double deformedBodyEquatorialRadius_;

    std::map< int, std::pair< double, double > > displacementLoveNumbers_;



    int numberOfBodies_;

    double currentTime_;

    std::vector< Eigen::Vector3d > relativeBodyStates_;

    std::vector< double > gravitationalParameterRatios_;
};

} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_BASICTIDALBODYDEFORMATION_H
