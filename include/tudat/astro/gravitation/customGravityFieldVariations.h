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

#include <functional>
#include <string>

#include "tudat/astro/gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace gravitation
{

//! Class for variations in spherical harmonic coefficients, directly provided as a function of time
class CustomGravityFieldVariations : public GravityFieldVariations
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param customCorrectionFunction Function returning the (cosine, sine) coefficient corrections as a function of time, as
     * matrix blocks starting at minimumDegree and minimumOrder
     * \param minimumDegree Degree where the correction blocks start
     * \param minimumOrder Order where the correction blocks start
     * \param maximumDegree Degree where the correction blocks end
     * \param maximumOrder Order where the correction blocks end
     * \param customCorrectionTimeDerivativeFunction Function returning the time derivatives of the (cosine, sine) coefficient
     * corrections as a function of time, with the same blocks as customCorrectionFunction (if empty, zero time derivatives are used,
     * with a warning)
     */
    CustomGravityFieldVariations(
            const std::function< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( const double ) > customCorrectionFunction,
            const int minimumDegree,
            const int minimumOrder,
            const int maximumDegree,
            const int maximumOrder,
            const std::function< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( const double ) > customCorrectionTimeDerivativeFunction =
                    nullptr ):
        GravityFieldVariations( minimumDegree, minimumOrder, maximumDegree, maximumOrder ),
        customCorrectionFunction_( customCorrectionFunction ),
        customCorrectionTimeDerivativeFunction_( customCorrectionTimeDerivativeFunction )
    {
        if( customCorrectionFunction_ == nullptr )
        {
            throw std::runtime_error( "Error when creating custom gravity field variations: no correction function provided." );
        }
    }

    //! Function for calculating spherical harmonic coefficient corrections.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( const double time ) override
    {
        return checkCorrectionBlockSizes( customCorrectionFunction_( time ), "corrections" );
    }

    //! Function for calculating time derivatives of spherical harmonic coefficient corrections (zero with warning if not defined).
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrectionsTimeDerivative( const double time ) override
    {
        if( customCorrectionTimeDerivativeFunction_ == nullptr )
        {
            return GravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative( time );
        }
        return checkCorrectionBlockSizes( customCorrectionTimeDerivativeFunction_( time ), "correction time derivatives" );
    }

private:
    //! Function to check whether the size of the user-defined correction blocks matches the degrees and orders of the variation
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > checkCorrectionBlockSizes(
            const std::pair< Eigen::MatrixXd, Eigen::MatrixXd >& corrections,
            const std::string& correctionType )
    {
        if( corrections.first.rows( ) != numberOfDegrees_ || corrections.first.cols( ) != numberOfOrders_ ||
            corrections.second.rows( ) != numberOfDegrees_ || corrections.second.cols( ) != numberOfOrders_ )
        {
            throw std::runtime_error( "Error in custom gravity field variations: size of " + correctionType +
                                      " is inconsistent with minimum and maximum degree and order." );
        }
        return corrections;
    }

    //! Function returning the (cosine, sine) coefficient corrections as a function of time
    std::function< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( const double ) > customCorrectionFunction_;

    //! Function returning the time derivatives of the (cosine, sine) coefficient corrections as a function of time
    std::function< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( const double ) > customCorrectionTimeDerivativeFunction_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_CUSTOMGRAVITYFIELDVARIATIONS_H
