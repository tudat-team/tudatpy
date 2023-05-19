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

#include <tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h>

namespace tudat
{

namespace aerodynamics
{

Eigen::Vector3d AerodynamicMomentContributionInterface::getMomentCoefficientsCorrection(
    const Eigen::Vector3d& momentReferencePoint,
    const Eigen::Vector3d& forceCoefficients,
    const double referenceLength )
{
    return signMultiplier_ * ( bodyFixedToMomentFrameRotation_( ) * ( momentReferencePoint - centerOfMassPosition_( ) ) ).cross(
        forceToMomentFrameRotation_( ) * forceCoefficients ) / referenceLength;
}

void AerodynamicCoefficientInterface::updateFullCurrentCoefficients(
    const std::vector< double >& independentVariables,
    const std::map< std::string, std::vector< double > >& controlSurfaceIndependentVariables,
    const double currentTime,
    const bool addMomentContributionIfPresent )
{
    updateCurrentCoefficients( independentVariables, currentTime );

    if( controlSurfaceIndependentVariables.size( ) != 0 )
    {
        currentControlSurfaceFreeForceCoefficients_ = currentForceCoefficients_;
        currentControlSurfaceFreeMomentCoefficients_ = currentMomentCoefficients_;

        for ( std::map<std::string, std::vector<double> >::const_iterator controlSurfaceIterator =
            controlSurfaceIndependentVariables.begin( );
              controlSurfaceIterator != controlSurfaceIndependentVariables.end( );
              controlSurfaceIterator++ )
        {
            updateCurrentControlSurfaceCoefficientsCoefficients(
                controlSurfaceIterator->first, controlSurfaceIterator->second );
        }
    }

    if( momentContributionInterface_ != nullptr && addMomentContributionIfPresent )
    {
        currentMomentCoefficients_ += momentContributionInterface_->getMomentCoefficientsCorrection(
            momentReferencePoint_, currentForceCoefficients_, referenceLength_ );
    }
}

} // namespace aerodynamics

} // namespace tudat
