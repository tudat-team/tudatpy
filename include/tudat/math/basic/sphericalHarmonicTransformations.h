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

#ifndef TUDAT_SPHERICALHARMONICTRANSFORMATIONS_H
#define TUDAT_SPHERICALHARMONICTRANSFORMATIONS_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "Tudat/Mathematics/BasicMathematics/wignerDMatrices.h"

namespace tudat
{

namespace basic_mathematics
{

//! Class to compute the transformation of spherical harmonic coefficients under a change of reference frame orientation
class SphericalHarmonicTransformationCache
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param maxiumDegree Maximum degree of spherical harmonic expansion that will be handled by object
     * \param maximumOrder Maximum order of spherical harmonic expansion that will be handled by object
     */
    SphericalHarmonicTransformationCache(
            const int maxiumDegree,
            const int maximumOrder ):
        maximumDegree_( maxiumDegree ), maximumOrder_( maximumOrder ), updatePartials_( false )
    {
        wignerDMatricesCache_ = boost::make_shared< WignerDMatricesCache >( maxiumDegree );
    }

    //! Function to update Wigner D-matrices for current orientation, parameterized by Cayley-Klein parameters
    /*!
     * Function to update Wigner D-matrices for current orientation, parameterized by Cayley-Klein parameters
     * \param cayleyKleinA Cayley-Klein parameter a denoting current orientation
     * \param cayleyKleinB Cayley-Klein parameter b denoting current orientation
     */
    void updateFromCayleyKleinParameters(
            const std::complex< double > cayleyKleinA,
            const std::complex< double > cayleyKleinB );

    //! Function to update Wigner D-matrices for current orientation, parameterized by quaternion
    /*!
     *  Function to update Wigner D-matrices for current orientation, parameterized by quaternion
     * \param rotationQuaternion Quaternion denoting current orientation
     */
    void updateFromQuaternion(
            const Eigen::Quaterniond& rotationQuaternion );

    //! Function to update Wigner D-matrices for current orientation, parameterized by 3-1-3 Euler angles
    /*!
     * Function to update Wigner D-matrices for current orientation, parameterized by 3-1-3 Euler angles
     * \param firstZRotation First z-axis rotation
     * \param XRotation X-axis rotation
     * \param secondZRotation Second z-axis rotation
     */
    void updateFrom313EulerAngles(
            const double firstZRotation, const double XRotation, const double secondZRotation )
    {
        updateFromQuaternion(
                    Eigen::Quaterniond( Eigen::AngleAxisd( -secondZRotation, Eigen::Vector3d::UnitZ( ) ) *
                                        Eigen::AngleAxisd( -XRotation, Eigen::Vector3d::UnitX( ) ) *
                                        Eigen::AngleAxisd( -firstZRotation, Eigen::Vector3d::UnitZ( ) ) ) );
    }

    //! Function to transform spherical harmonic coefficients using current wignerDMatricesCache_
    /*!
     * Function to transform spherical harmonic coefficients using current wignerDMatricesCache_
     * \param originalCosineCoefficients Original cosine coefficients
     * \param originalSineCoefficients Original sine coefficients
     * \param currentCosineCoefficients Transformed cosine coefficients (returned by reference)
     * \param currentSineCoefficients Transformed sine coefficients (returned by reference)
     * \param areCoefficientsNormalized Boolean denoting whether coefficients ar fully normalized or unnormalized
     */
    void transformCoefficientsAtDegree(
            const Eigen::MatrixXd& originalCosineCoefficients,
            const Eigen::MatrixXd& originalSineCoefficients,
            Eigen::MatrixXd& currentCosineCoefficients,
            Eigen::MatrixXd& currentSineCoefficients,
            const bool areCoefficientsNormalized = 1 );

    //! Function to set boolean denoting whether partial derivatives of coefficients are to be computed
    /*!
     * Function to set boolean denoting whether partial derivatives of coefficients are to be computed
     */
    void setUpdatePartials(
            const bool updatePartials = 1 )
    {
        updatePartials_ = updatePartials;
    }


private:

    //! Object used to compute Wigner D-matrices
    boost::shared_ptr< WignerDMatricesCache > wignerDMatricesCache_;

    //! Maximum degree of spherical harmonic expansion that will be handled by object
    int maximumDegree_;

    //! Maximum degree of spherical harmonic expansion that will be handled by object
    int maximumOrder_;

    //! Boolean denoting whether partial derivatives of coefficients are to be computed
    bool updatePartials_;

};

}

}

#endif // TUDAT_SPHERICALHARMONICTRANSFORMATIONS_H
