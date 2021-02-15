/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CUSTOMGRAVITYFIELDVARIATIONS_H
#define TUDAT_CUSTOMGRAVITYFIELDVARIATIONS_H

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{

namespace gravitation
{

//! Class for variations in spherical harmonic coefficients, directly provided as a function
class CustomGravityFieldVariations: public GravityFieldVariations
{

public:


    CustomGravityFieldVariations(
            std::function< double, std::pair< Eigen::MatrixXd >( const double ) > customCorrectionFunction,
            const int minimumDegree, const int minimumOrder ):
    GravityFieldVariations( minimumDegree, minimumOrder, -1, -1 ),
    customCorrectionFunction_( customCorrectionFunction ){ }



    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections(
            const double time )
    {
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > currentCorrections =
                customCorrectionFunction_( time );
        if( numberOfDegrees_ == -1 || numberOfOrders_ == -1 )
        {
            maximumDegree_ = minimumDegree_ + currentCorrections.first.rows( ) - 1;
            maximumOrder_ = minimumOrder_ + currentCorrections.first.cols( ) - 1;
            numberOfDegrees_ = maximumDegree_ - minimumDegree_ + 1;
            numberOfOrders_ = maximumOrder_ - minimumOrder_ + 1;
            lastCosineCorrection_.setZero( maximumDegree_ + 1, maximumOrder_ + 1 );
            lastSineCorrection_.setZero( maximumDegree_ + 1, maximumOrder_ + 1 );
        }
        return currentCorrections;
    }


private:

    std::function< double, std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( const double ) > customCorrectionFunction_;

};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_CUSTOMGRAVITYFIELDVARIATIONS_H
