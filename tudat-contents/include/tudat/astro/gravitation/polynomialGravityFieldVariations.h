/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PPOLYNOMIALGRAVITYFIELDVARIATIONS_H
#define TUDAT_PPOLYNOMIALGRAVITYFIELDVARIATIONS_H

#include <Eigen/Geometry>

#include "tudat/astro/gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace gravitation
{

class PolynomialGravityFieldVariations : public GravityFieldVariations
{
public:
    PolynomialGravityFieldVariations( const std::map< int, Eigen::MatrixXd >& cosineAmplitudes,
                                      const std::map< int, Eigen::MatrixXd >& sineAmplitudes,
                                      const double referenceEpoch,
                                      const int minimumDegree = 2,
                                      const int minimumOrder = 0 );

    virtual ~PolynomialGravityFieldVariations( ) { }

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( const double time );

    std::map< int, Eigen::MatrixXd > getCosineAmplitudes( )
    {
        return cosineAmplitudes_;
    }

    std::map< int, Eigen::MatrixXd > getSineAmplitudes( )
    {
        return sineAmplitudes_;
    }

    void resetCosineAmplitudes( const std::map< int, Eigen::MatrixXd > cosineAmplitudes )
    {
        cosineAmplitudes_ = cosineAmplitudes;
    }

    void resetSineAmplitudes( const std::map< int, Eigen::MatrixXd > sineAmplitudes )
    {
        sineAmplitudes_ = sineAmplitudes;
    }

    std::map< int, Eigen::MatrixXd >& getCosineAmplitudesReference( )
    {
        return cosineAmplitudes_;
    }

    std::map< int, Eigen::MatrixXd >& getSineAmplitudesReference( )
    {
        return sineAmplitudes_;
    }

    double getReferenceEpoch( )
    {
        return referenceEpoch_;
    }

protected:
    std::map< int, Eigen::MatrixXd > cosineAmplitudes_;

    std::map< int, Eigen::MatrixXd > sineAmplitudes_;

    const double referenceEpoch_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_PPOLYNOMIALGRAVITYFIELDVARIATIONS_H
