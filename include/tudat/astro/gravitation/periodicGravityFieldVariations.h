/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PERIODICGRAVITYFIELDVARIATIONS_H
#define TUDAT_PERIODICGRAVITYFIELDVARIATIONS_H

#include <Eigen/Geometry>

#include "tudat/astro/gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace gravitation
{

class PeriodicGravityFieldVariations: public GravityFieldVariations
{
public:
    PeriodicGravityFieldVariations(
            const std::vector< Eigen::MatrixXd >& cosineShAmplitudesCosineTime,
            const std::vector< Eigen::MatrixXd >& cosineShAmplitudesSineTime,
            const std::vector< Eigen::MatrixXd >& sineShAmplitudesCosineTime,
            const std::vector< Eigen::MatrixXd >& sineShAmplitudesSineTime,
            const std::vector< double >& frequencies,
            const double referenceEpoch,
            const int minimumDegree = 2,
            const int minimumOrder = 0 );

    virtual ~PeriodicGravityFieldVariations( ){ }

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections(
            const double time );


    std::vector< Eigen::MatrixXd > getCosineShAmplitudesCosineTime( )
    {
        return cosineShAmplitudesCosineTime_;
    }

    std::vector< Eigen::MatrixXd > getCosineShAmplitudesSineTime( )
    {
        return cosineShAmplitudesSineTime_;
    }

    std::vector< Eigen::MatrixXd > getSineShAmplitudesCosineTime( )
    {
        return sineShAmplitudesCosineTime_;
    }

    std::vector< Eigen::MatrixXd > getSineShAmplitudesSineTime( )
    {
        return sineShAmplitudesSineTime_;
    }

    void resetCosineShAmplitudesCosineTime( const std::vector< Eigen::MatrixXd >& cosineShAmplitudesCosineTime )
    {
        cosineShAmplitudesCosineTime_ = cosineShAmplitudesCosineTime;
    }

    void resetCosineShAmplitudesSineTime( const std::vector< Eigen::MatrixXd >& cosineShAmplitudesSineTime )
    {
        cosineShAmplitudesSineTime_ = cosineShAmplitudesSineTime;
    }

    void resetSineShAmplitudesCosineTime( const std::vector< Eigen::MatrixXd >& sineShAmplitudesCosineTime )
    {
        sineShAmplitudesCosineTime_ = sineShAmplitudesCosineTime;
    }

    void resetSineShAmplitudesSineTime( const std::vector< Eigen::MatrixXd >& sineShAmplitudesSineTime )
    {
        sineShAmplitudesSineTime_ = sineShAmplitudesSineTime;
    }

    std::vector< double > getFrequencies( )
    {
        return frequencies_;
    }

    double getReferenceEpoch( )
    {
        return referenceEpoch_;
    }

protected:

    std::vector< Eigen::MatrixXd > cosineShAmplitudesCosineTime_;

    std::vector< Eigen::MatrixXd > cosineShAmplitudesSineTime_;

    std::vector< Eigen::MatrixXd > sineShAmplitudesCosineTime_;

    std::vector< Eigen::MatrixXd > sineShAmplitudesSineTime_;

    const std::vector< double > frequencies_;

    const double referenceEpoch_;

};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_PERIODICGRAVITYFIELDVARIATIONS_H
