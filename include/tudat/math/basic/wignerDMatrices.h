/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Eberly, D. Spherical Harmonics. Help documentation of Geometric Tools, 2008. Available at
 *        URL http://www.geometrictools.com/Documentation/SphericalHarmonics.pdf. Last access:
 *        09-09-2012.
 *      Heiskanen, W.A., Moritz, H. Physical geodesy. Freeman, 1967.
 *      Holmes, S.A., Featherstone, W.E. A unified approach to the Clenshaw summation and the
 *        recursive computation of very high degree and order normalised associated Legendre
 *        functions. Journal of Geodesy, 76(5):279-299, 2002.
 *      Vallado, D. and McClain, W. Fundamentals of astrodynammics and applications. Microcosm
 *        Press, 2001.
 *      Weisstein, E.W. Associated Legendre Polynomial, 2012.
 *        URL http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html. Last access:
 *        12-09-2012.
 *
 *    Notes
 *      For information on how the caching mechanism works, please contact S. Billemont
 *      (S.Billemont@studelft.tudelft.nl).
 *
 */

#ifndef TUDAT_WIGNER_D_MATRIXRES_H
#define TUDAT_WIGNER_D_MATRIXRES_H

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{

class WignerDMatricesCache
{
public:
    WignerDMatricesCache( const int maximumDegree );

    ~WignerDMatricesCache( ){ }

    void updateMatrices( const std::complex< double > cayleyKleinA, const std::complex< double > cayleyKleinB );

    std::complex< double > getWignerDCoefficient( const int degree, const int originalOrder, const int newOrder )
    {
        return wignerDMatrices_[ degree ]( originalOrder + degree, newOrder + degree );
    }

    Eigen::MatrixXcd& getWignerDMatrix( const int degree )
    {
        return wignerDMatrices_[ degree ];
    }

    std::vector< Eigen::MatrixXcd >& getWignerDMatrices( )
    {
        return wignerDMatrices_;
    }

private:

    void computeCoefficients( );

    std::vector< Eigen::MatrixXd > coefficientsIndexMinusOne_;

    std::vector< Eigen::MatrixXd > coefficientsIndexZero_;

    std::vector< Eigen::MatrixXd > coefficientsIndexOne_;

    int maximumDegree_;

    std::vector< Eigen::MatrixXcd > wignerDMatrices_;

    std::complex< double > currentCayleyKleinA_;

    std::complex< double > currentCayleyKleinB_;

    std::complex< double > currentCayleyKleinAConjugate_;

    std::complex< double > currentCayleyKleinBConjugate_;
};

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_WIGNER_D_MATRIXRES_H
