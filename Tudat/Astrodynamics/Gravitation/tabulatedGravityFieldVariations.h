/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TABULATEDGRAVITYFIELDVARIATIONS_H
#define TUDAT_TABULATEDGRAVITYFIELDVARIATIONS_H

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"

namespace tudat
{

namespace gravitation
{

//! Class for variations in spherical harmonic coefficients, directly interpolated from tables.
/*!
 *  Class for variations in spherical harmonic coefficients, directly interpolated from tables.
 *  Maps of variatons of cosine and sine coefficients as a function of time are provided, as well
 *  as the type of interpolator that is to be used, from which interpolated variations are then
 *  returned at any time. The class may also be used for defining the complete gravity field
 *  (i.e. not only the variations) as a function of time when the constant coefficient contribution
 *  (i.e. not from a GravityFieldVariations) is zero for all degrees and orders.
 */
class TabulatedGravityFieldVariations: public GravityFieldVariations
{

public:

    //! Class constructor
    /*!
     *  Class constructor, takes the maps of variations, as well as the minimum degree and order
     *  of the corrections, the corrections are applied as a rectangular block in degree and order,
     *  with the the size of the blocks determined by the contents of the maps. The matrix size of
     *  all map entries must be equal for both sine and cosine values.
     *  \param cosineCoefficientCorrections Map of cosine coefficient variations, with associated
     *  times as map key.
     *  \param sineCoefficientCorrections Map of sine coefficient variations, with associated
     *  times as map key.
     *  \param minimumDegree Minimum degree of correction block.
     *  \param minimumOrder Minimum order of correction block.
     *  \param interpolatorType Type of interpolator to use for calculating coefficients at any
     *  time between tabulated times.
     */
    TabulatedGravityFieldVariations(
            const std::map< double, Eigen::MatrixXd >& cosineCoefficientCorrections,
            const std::map< double, Eigen::MatrixXd >& sineCoefficientCorrections,
            const int minimumDegree, const int minimumOrder,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorType =
            std::make_shared< interpolators::InterpolatorSettings >(
                interpolators::linear_interpolator, interpolators::huntingAlgorithm ) );

    //! Function to (re)set the tabulated spherical harmonic coefficients.
    /*!
     *  Function to (re)set the tabulated spherical harmonic coefficients.
     *  Creates interpolator for returning the coefficients as a function of time.
     *  \param cosineCoefficientCorrections Map of cosine coefficient variations,
     *  with associated times as map key.
     *  \param sineCoefficientCorrections Map of sine coefficient variations,
     *  with associated times as map key.
     */
    void resetCoefficientInterpolator(
            const std::map< double, Eigen::MatrixXd >& cosineCoefficientCorrections,
            const std::map< double, Eigen::MatrixXd >& sineCoefficientCorrections );

    //! Function for calculating corrections by interpolating tabulated corrections.
    /*!
     *  Function for calculating corrections by interpolating tabulated corrections.
     *  \param time Time at which variations are to be calculated.
     *  \return Pair of  matrices containing variations in (cosine, sine) coefficients at
     *  block positions in total matrices defined by
     *  minimumDegree_, minimumOrder_, numberOfDegrees_, numberOfOrders_;
     */
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections(
            const double time );

    //! Function to return map of cosine coefficient variations, with associated times as map key.
    /*!
     *  Function to return map of cosine coefficient variations, with associated times as map key.
     * \return Map of cosine coefficient variations, with associated times as map key.
     */
    std::map< double, Eigen::MatrixXd > getCosineCoefficientCorrections( )
    {
        return cosineCoefficientCorrections_;
    }

    //! Function to return map of sine coefficient variations, with associated times as map key.
    /*!
     *  Function to return map of sine coefficient variations, with associated times as map key.
     * \return Map of sine coefficient variations, with associated times as map key.
     */
    std::map< double, Eigen::MatrixXd > getSineCoefficientCorrections( )
    {
        return sineCoefficientCorrections_;
    }

    //! Function to return type of interpolator to use for calculating coefficients.
    /*!
     *  Function to return type of interpolator to use for calculating coefficients at any time
     *  between tabulated times.
     * \return Type of interpolator to use for calculating coefficients
     */
    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( )
    {
        return interpolatorType_;
    }

private:

    //! Type of interpolator to use for calculating coefficients at any time.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorType_;

    //! Map of cosine coefficient variations, with associated times as map key.
    std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections_;

    //! Map of sine coefficient variations, with associated times as map key.
    std::map< double, Eigen::MatrixXd > sineCoefficientCorrections_;

    //! Interpolator for determining concatenated [cosine|sine] coefficients from tabulated values
    //! at discrete times.
    /*!
     *  Interpolator for determining concatenated [cosine|sine] coefficients from tabulated values
     *  at discrete times. The coefficients are interpolated concurrently and then split to reduce
     *  the number of calls to the interpolator.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    variationInterpolator_;

};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_TABULATEDGRAVITYFIELDVARIATIONS_H
