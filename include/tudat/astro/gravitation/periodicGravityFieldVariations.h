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
            const std::vector< Eigen::MatrixXd >& cosineAmplitudes,
            const std::vector< Eigen::MatrixXd >& sineAmplitudes,
            const std::vector< double >& frequencies,
            const std::vector< double >& phases,
            const double referenceEpoch,
            const int minimumDegree = 2,
            const int minimumOrder = 0 );

    virtual ~PeriodicGravityFieldVariations( ){ }

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections(
            const double time );


protected:

    const std::vector< Eigen::MatrixXd > cosineAmplitudes_;

    const std::vector< Eigen::MatrixXd > sineAmplitudes_;

    const std::vector< double > frequencies_;

    const std::vector< double > phases_;

    const double referenceEpoch_;

};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_PERIODICGRAVITYFIELDVARIATIONS_H
