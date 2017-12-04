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
 *      Varschalovich, et al. Quantum Theory of Angular Momentum. World Scientific, February 1988
 *      Boué, (2017) The two rigid body interaction using angular momentum theory formulae, CMDA 128(2-3):261-273
 *
 */

#ifndef TUDAT_WIGNER_D_MATRIXRES_H
#define TUDAT_WIGNER_D_MATRIXRES_H

#include <Eigen/Core>

namespace tudat
{

namespace basic_mathematics
{

//! Class to compute Wigner D-Matrices, used to transform spherical harmonic coefficients
/*!
 *  Class to compute Wigner D-Matrices, used to transform spherical harmonic coefficients under a rotation of the reference
 *  frame. The formulation is in terms of Cayley-Klein parameters, and follows the same recursive formulation as Boué (2016)
 */
class WignerDMatricesCache
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param maximumDegree Maximum degree for which Wigner-D matrix is to be computed
     */
    WignerDMatricesCache( const int maximumDegree );

    //! Destructor
    ~WignerDMatricesCache( ){ }

    //! Function to update contents of this object to new orientation
    /*!
     * Function to update contents of this object to new orientation, defined by Cayley-Klein parameters a and b. The
     * parameters a and b are related to teh quaternion q by a = q0 - i q3, b = q1 - i q2.
     * Calling this function updates the wignerDMatrices_ member variable
     * \param cayleyKleinA Cayley-Klein parameters a
     * \param cayleyKleinB Cayley-Klein parameters b
     */
    void updateMatrices( const std::complex< double > cayleyKleinA, const std::complex< double > cayleyKleinB );

    //! Function to retrieve single component of Wigner D-matrix
    /*!
     * Function to retrieve single component of Wigner D-matrix
     * \param degree Degree for which coefficient is to be retrieved (denoted l by Boué)
     * \param originalOrder Original order for which coefficient is to be retrieved (denoted m by Boué)
     * \param newOrder New order for which coefficient is to be retrieved (denoted m' by Boué)
     * \return Single component of Wigner D-matrix
     */
    std::complex< double > getWignerDCoefficient( const int degree, const int originalOrder, const int newOrder )
    {
        return wignerDMatrices_[ degree ]( originalOrder + degree, newOrder + degree );
    }

    //! Function to retrieve single Wigner D-matrix
    /*!
     * Function to retrieve single Wigner D-matrix
     * \param degree  Degree for which matrix is to be retrieved (denoted l by Boué)
     * \return Wigner D-matrix at requested degree
     */
    Eigen::MatrixXcd& getWignerDMatrix( const int degree )
    {
        return wignerDMatrices_[ degree ];
    }

    //! Function to retrieve list of Wigner D-matrices
    /*!
     *  Function to retrieve list of Wigner D-matrices
     *  \return List of Wigner D-matrices
     */
    std::vector< Eigen::MatrixXcd >& getWignerDMatrices( )
    {
        return wignerDMatrices_;
    }

private:

    //! Function to precompute the coefficients used on the recursive formulation for Wigner D-matrices
    void computeCoefficients( );

    //! Coefficients used in the recursive formulation for Wigner D-matrices
    std::vector< Eigen::MatrixXd > coefficientsIndexMinusOne_;

    //! Coefficients used in the recursive formulation for Wigner D-matrices
    std::vector< Eigen::MatrixXd > coefficientsIndexZero_;

    //! Coefficients used in the recursive formulation for Wigner D-matrices
    std::vector< Eigen::MatrixXd > coefficientsIndexOne_;

    //! Maximum degree for which Wigner-D matrix is to be computed
    int maximumDegree_;

    //! List of Wigner D-matrices
    std::vector< Eigen::MatrixXcd > wignerDMatrices_;

    //! Cayley-Klein parameter a for current orientation.
    std::complex< double > currentCayleyKleinA_;

    //! Cayley-Klein parameter b for current orientation.
    std::complex< double > currentCayleyKleinB_;

    //! Conjugate of Cayley-Klein parameter a for current orientation.
    std::complex< double > currentCayleyKleinAConjugate_;

    //! Conjugate of Cayley-Klein parameter b for current orientation.
    std::complex< double > currentCayleyKleinBConjugate_;
};

} // namespace basic_mathematics

} // namespace tudat

#endif // TUDAT_WIGNER_D_MATRIXRES_H
