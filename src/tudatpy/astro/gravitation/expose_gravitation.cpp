/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_gravitation.h"

#include <tuple>

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/gravitation/sphericalHarmonicsGravityField.h>
#include <tudat/math/basic/legendrePolynomials.h>
#include <tudat/math/basic/sphericalHarmonicTransformations.h>
#include <tudat/math/basic/wignerDMatrices.h>

namespace py = pybind11;
namespace tg = tudat::gravitation;
namespace tbm = tudat::basic_mathematics;

namespace tudat
{
namespace gravitation
{
std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > getDegreeTwoSphericalHarmonicCoefficientsPy( const Eigen::Matrix3d inertiaTensor,
                                                                                                    const double bodyGravitationalParameter,
                                                                                                    const double referenceRadius,
                                                                                                    const bool useNormalizedCoefficients )
{
    return tg::getDegreeTwoSphericalHarmonicCoefficients(
            inertiaTensor, bodyGravitationalParameter, referenceRadius, 2, useNormalizedCoefficients );
}

std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > transformSphericalHarmonicCoefficientsWithWignerDPy(
        const Eigen::MatrixXd& originalCosineCoefficients,
        const Eigen::MatrixXd& originalSineCoefficients,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        const bool areCoefficientsNormalized )
{
    Eigen::MatrixXd transformedCosineCoefficients =
            Eigen::MatrixXd::Zero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
    Eigen::MatrixXd transformedSineCoefficients =
            Eigen::MatrixXd::Zero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );
    basic_mathematics::transformSphericalHarmonicCoefficientsWithWignerD( originalCosineCoefficients,
                                                                          originalSineCoefficients,
                                                                          wignerMatrices,
                                                                          transformedCosineCoefficients,
                                                                          transformedSineCoefficients,
                                                                          areCoefficientsNormalized );
    return std::make_tuple( transformedCosineCoefficients, transformedSineCoefficients );
}

std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > transformSphericalHarmonicCoefficientsWithCachePy(
        basic_mathematics::SphericalHarmonicTransformationCache& transformationCache,
        const Eigen::MatrixXd& originalCosineCoefficients,
        const Eigen::MatrixXd& originalSineCoefficients,
        const bool areCoefficientsNormalized )
{
    Eigen::MatrixXd transformedCosineCoefficients =
            Eigen::MatrixXd::Zero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
    Eigen::MatrixXd transformedSineCoefficients =
            Eigen::MatrixXd::Zero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );
    transformationCache.transformCoefficientsAtDegree( originalCosineCoefficients,
                                                       originalSineCoefficients,
                                                       transformedCosineCoefficients,
                                                       transformedSineCoefficients,
                                                       areCoefficientsNormalized );
    return std::make_tuple( transformedCosineCoefficients, transformedSineCoefficients );
}

void updateSphericalHarmonicTransformationCacheFromQuaternionVectorPy(
        basic_mathematics::SphericalHarmonicTransformationCache& transformationCache,
        const Eigen::Vector4d& rotationQuaternion )
{
    transformationCache.updateFromQuaternion(
            Eigen::Quaterniond( rotationQuaternion( 0 ), rotationQuaternion( 1 ), rotationQuaternion( 2 ), rotationQuaternion( 3 ) ) );
}
}  // namespace gravitation
}  // namespace tudat

namespace tudatpy
{
namespace astro
{
namespace gravitation
{

void expose_gravitation( py::module& m )
{
    py::class_< tbm::WignerDMatricesCache, std::shared_ptr< tbm::WignerDMatricesCache > >( m,
                                                                                           "WignerDMatricesCache",
                                                                                           R"doc(

         Cache for Wigner D-matrices used to rotate spherical harmonic coefficients.

         This class computes and stores Wigner D-matrices up to a requested maximum degree for a rotation
         parameterized by Cayley-Klein parameters. These matrices are used in the spherical-harmonic coefficient
         transformation algorithms of :cite:t:`boue2017`, using the angular-momentum conventions of
         :cite:t:`varshalovich1988`.

         The implemented algorithm is the recursive Cayley-Klein formulation used by :cite:t:`boue2017`.
         With Cayley-Klein parameters :math:`a` and :math:`b`, the recursion starts from

         .. math::

            D^0_{0,0}=1

         and the degree-one seed

         .. math::

            \begin{aligned}
            D^1_{1,1}&=a^2,&
            D^1_{1,0}&=-\sqrt{2}a b^*,&
            D^1_{1,-1}&=(b^*)^2,\\
            D^1_{0,1}&=\sqrt{2}ab,&
            D^1_{0,0}&=|a|^2-|b|^2,&
            D^1_{0,-1}&=-\sqrt{2}a^*b^*,\\
            D^1_{-1,1}&=b^2,&
            D^1_{-1,0}&=\sqrt{2}a^*b,&
            D^1_{-1,-1}&=(a^*)^2 .
            \end{aligned}

         For :math:`l>1`, non-negative row orders are computed from

         .. math::

            D^l_{m,k} =
            q^-_{lmk}D^1_{1,1}D^{l-1}_{m-1,k-1}
            +q^0_{lmk}D^1_{1,0}D^{l-1}_{m-1,k}
            +q^+_{lmk}D^1_{1,-1}D^{l-1}_{m-1,k+1}

         with

         .. math::

            q^-_{lmk}=\sqrt{\frac{(l+k)(l+k-1)}{(l+m)(l+m-1)}},\quad
            q^0_{lmk}=\sqrt{\frac{2(l+k)(l-k)}{(l+m)(l+m-1)}},\quad
            q^+_{lmk}=\sqrt{\frac{(l-k)(l-k-1)}{(l+m)(l+m-1)}} .

         The remaining rows are filled from the conjugate symmetry

         .. math::

            D^l_{m,k}=(-1)^{m-k}\left(D^l_{-m,-k}\right)^* .

         Cartesian angular-momentum-operator values are computed from the cached matrices as

         .. math::

            \hat{\mathcal{J}}_zD^l_{m,k}=-i mD^l_{m,k}

         and, with
         :math:`A^+_{lm}=i\sqrt{(l(l+1)-m(m-1))/2}\,D^l_{m-1,k}` and
         :math:`A^-_{lm}=-i\sqrt{(l(l+1)-m(m+1))/2}\,D^l_{m+1,k}`,

         .. math::

            \hat{\mathcal{J}}_xD^l_{m,k}=\frac{A^-_{lm}-A^+_{lm}}{\sqrt{2}},
            \qquad
            \hat{\mathcal{J}}_yD^l_{m,k}=i\frac{A^-_{lm}+A^+_{lm}}{\sqrt{2}} .





      )doc" )
            .def( py::init< const int >( ), py::arg( "maximum_degree" ) )
            .def( "update_matrices",
                  &tbm::WignerDMatricesCache::updateMatrices,
                  py::arg( "cayley_klein_a" ),
                  py::arg( "cayley_klein_b" ),
                  R"doc(

 Updates the cached Wigner D-matrices from Cayley-Klein parameters.

 Parameters
 ----------
 cayley_klein_a : complex
     Cayley-Klein parameter :math:`a`.
 cayley_klein_b : complex
     Cayley-Klein parameter :math:`b`.





     )doc" )
            .def( "get_coefficient",
                  &tbm::WignerDMatricesCache::getWignerDCoefficient,
                  py::arg( "degree" ),
                  py::arg( "original_order" ),
                  py::arg( "new_order" ),
                  R"doc(

 Returns a single Wigner D-matrix coefficient.

 Parameters
 ----------
 degree : int
     Spherical harmonic degree.
 original_order : int
     Original spherical harmonic order.
 new_order : int
     Transformed spherical harmonic order.
 Returns
 -------
 complex
     Requested Wigner D coefficient.





     )doc" )
            .def(
                    "get_matrix",
                    []( tbm::WignerDMatricesCache& cache, const int degree ) {
                        return Eigen::MatrixXcd( cache.getWignerDMatrix( degree ) );
                    },
                    py::arg( "degree" ),
                    R"doc(

 Returns the Wigner D-matrix for a single degree.

 Parameters
 ----------
 degree : int
     Spherical harmonic degree.
 Returns
 -------
 numpy.ndarray
     Complex Wigner D-matrix for the requested degree.





     )doc" )
            .def(
                    "get_matrices",
                    []( tbm::WignerDMatricesCache& cache ) { return cache.getWignerDMatrices( ); },
                    R"doc(

 Returns all cached Wigner D-matrices.

 Returns
 -------
 list[numpy.ndarray]
     Complex Wigner D-matrices from degree zero up to the maximum degree.





     )doc" )
            .def( "get_angular_momentum_operator_on_coefficient",
                  &tbm::WignerDMatricesCache::getAngularMomentumOperatorOnWignerCoefficient,
                  py::arg( "degree" ),
                  py::arg( "original_order" ),
                  py::arg( "new_order" ),
                  R"doc(

 Returns the Cartesian angular-momentum operator applied to a Wigner D coefficient.

 Parameters
 ----------
 degree : int
     Spherical harmonic degree.
 original_order : int
     Original spherical harmonic order.
 new_order : int
     Transformed spherical harmonic order.
 Returns
 -------
 numpy.ndarray
     Complex vector with the three Cartesian angular-momentum-operator components.





     )doc" );

    py::class_< tbm::SphericalHarmonicTransformationCache, std::shared_ptr< tbm::SphericalHarmonicTransformationCache > >(
            m,
            "SphericalHarmonicTransformationCache",
            R"doc(

         Cache for rotating spherical harmonic coefficients between reference frames.

         This class stores Wigner D-matrices for a specified maximum degree/order and applies them to real-valued
         spherical harmonic coefficient matrices. The implementation follows the coefficient transformation approach
         of :cite:t:`boue2017`.

         A call to one of the update functions first computes the Cayley-Klein parameters of the requested frame
         rotation and updates the internal :class:`~tudatpy.astro.gravitation.WignerDMatricesCache`. A subsequent
         call to :func:`~tudatpy.astro.gravitation.SphericalHarmonicTransformationCache.transform_coefficients`
         applies the cached Wigner D-matrices degree by degree to the real cosine/sine coefficient matrices.

         For each degree :math:`l` and target order :math:`m\geq 0`, define

         .. math::

            \eta_m =
            \begin{cases}
            1, & m=0\\
            1/\sqrt{2}, & m>0
            \end{cases}

         for normalized coefficients, and
         :math:`\eta_m=\sqrt{(l-m)!/(l+m)!}` for unnormalized coefficients. For :math:`k>0`, define
         :math:`\eta_{mk}=\sqrt{2}\eta_m` for normalized coefficients and
         :math:`\eta_{mk}=\sqrt{(l+k)!/(l-k)!}\eta_m` for unnormalized coefficients. The intermediate transformed
         coefficients are computed as

         .. math::

            \tilde{C}_{lm} =
            \eta_m\Re(D^l_{m,0})C_{l0}
            +\frac{1}{2}\sum_{k=1}^{l}\eta_{mk}
            \left[
            \left((-1)^k\Re(D^l_{m,k})+\Re(D^l_{m,-k})\right)C_{lk}
            +\left((-1)^k\Im(D^l_{m,k})-\Im(D^l_{m,-k})\right)S_{lk}
            \right]

         .. math::

            \tilde{S}_{lm} =
            \eta_m\Im(D^l_{m,0})C_{l0}
            +\frac{1}{2}\sum_{k=1}^{l}\eta_{mk}
            \left[
            \left((-1)^k\Im(D^l_{m,k})+\Im(D^l_{m,-k})\right)C_{lk}
            +\left(-(-1)^k\Re(D^l_{m,k})+\Re(D^l_{m,-k})\right)S_{lk}
            \right].

         The returned real coefficients are then

         .. math::

            C'_{lm}=(-1)^m(2-\delta_{0m})\tilde{C}_{lm},\qquad
            S'_{lm}=(-1)^{m+1}(2-\delta_{0m})\tilde{S}_{lm}.





      )doc" )
            .def( py::init< const int, const int >( ), py::arg( "maximum_degree" ), py::arg( "maximum_order" ) )
            .def( "update_from_cayley_klein_parameters",
                  &tbm::SphericalHarmonicTransformationCache::updateFromCayleyKleinParameters,
                  py::arg( "cayley_klein_a" ),
                  py::arg( "cayley_klein_b" ),
                  R"doc(

 Updates the transformation cache from Cayley-Klein parameters.

 Parameters
 ----------
 cayley_klein_a : complex
     Cayley-Klein parameter :math:`a`.
 cayley_klein_b : complex
     Cayley-Klein parameter :math:`b`.





     )doc" )
            .def( "update_from_quaternion",
                  &tg::updateSphericalHarmonicTransformationCacheFromQuaternionVectorPy,
                  py::arg( "rotation_quaternion" ),
                  R"doc(

 Updates the transformation cache from a quaternion.

 Parameters
 ----------
 rotation_quaternion : numpy.ndarray
     Quaternion as a scalar-first vector ``[w, x, y, z]``.





     )doc" )
            .def( "update_from_313_euler_angles",
                  &tbm::SphericalHarmonicTransformationCache::updateFrom313EulerAngles,
                  py::arg( "first_z_rotation" ),
                  py::arg( "x_rotation" ),
                  py::arg( "second_z_rotation" ),
                  R"doc(

 Updates the transformation cache from 3-1-3 Euler angles.

 Parameters
 ----------
 first_z_rotation : float
     First rotation angle about the z-axis, in radians.
 x_rotation : float
     Rotation angle about the x-axis, in radians.
 second_z_rotation : float
     Second rotation angle about the z-axis, in radians.





     )doc" )
            .def( "transform_coefficients",
                  &tg::transformSphericalHarmonicCoefficientsWithCachePy,
                  py::arg( "cosine_coefficients" ),
                  py::arg( "sine_coefficients" ),
                  py::arg( "are_coefficients_normalized" ) = true,
                  R"doc(

 Rotates spherical harmonic coefficients using the current cache orientation.

 Parameters
 ----------
 cosine_coefficients : numpy.ndarray
     Matrix of cosine spherical harmonic coefficients.
 sine_coefficients : numpy.ndarray
     Matrix of sine spherical harmonic coefficients.
 are_coefficients_normalized : bool, default=True
     Boolean denoting whether the input coefficients are fully normalized.
 Returns
 -------
 tuple[numpy.ndarray, numpy.ndarray]
     Tuple containing the rotated cosine and sine coefficient matrices.





     )doc" )
            .def( "set_update_partials",
                  &tbm::SphericalHarmonicTransformationCache::setUpdatePartials,
                  py::arg( "update_partials" ) = true,
                  R"doc(

 Sets whether coefficient partials are updated by the cache.

 Parameters
 ----------
 update_partials : bool, default=True
     Boolean denoting whether coefficient partials should be updated.





     )doc" )
            .def( "get_wigner_matrices_cache",
                  &tbm::SphericalHarmonicTransformationCache::getWignerDMatricesCache,
                  R"doc(

 Returns the Wigner D-matrix cache used internally.

 Returns
 -------
 WignerDMatricesCache
     Wigner D-matrix cache used by this spherical harmonic transformation cache.





     )doc" );

    m.def( "legendre_normalization_factor",
           &tbm::calculateLegendreGeodesyNormalizationFactor,
           py::arg( "degree" ),
           py::arg( "order" ),
           R"doc(

 Function to calculate the normalization factor for spherical harmonics at a given degree and order

 Function to calculate the normalization factor for spherical harmonics at a given degree and order.
 Specifically, this function returns the value :math:`\mathcal{N}_{lm}`, computed from:

 .. math::
     \mathcal{N}_{lm}=\sqrt{\frac{(2-\delta_{0m})(2l+1)(l-m)!)}{(l+m)!}}

 with :math:`\delta_{0m}` is the Kronecker Delta function. The following can be used such that the conversion between unnormalized and fully
 normalized spherical harmonic coefficients and Legendre polynomials can be computed from:

 .. math::
     [C,S]_{lm}&=\mathcal{N}_{lm}[\bar{C},\bar{S}]_{lm}\\
     \bar{P}_{lm}&=\mathcal{N}_{lm}{P}_{lm}

 with :math:`[C,S]_{lm}` denoting the unnormalized cosine or sine spherical harmonic coefficients at degree :math:`l` and order :math:`m`,
 and :math:`P_{lm}` and :math:`\bar{P}_{lm}` representing the unnormalized and normalized associated Legendre polynomials at degree :math:`l` and order :math:`m`.


 Parameters
 ----------
 degree : int
     Spherical harmonic degree :math:`l`
 order : int
     Spherical harmonic order :math:`m`
 Returns
 -------
 float
     Normalization factor :math:`\mathcal{N}_{lm}`






     )doc" );

    m.def( "normalize_spherical_harmonic_coefficients",
           py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd& >( &tbm::convertUnnormalizedToGeodesyNormalizedCoefficients ),
           py::arg( "unnormalized_cosine_coefficients" ),
           py::arg( "unnormalized_sine_coefficients" ),
           R"doc(

 Function to normalize spherical harmonic coefficients

 Function to normalize spherical harmonic coefficients, using the equations provided in the :func:`~tudatpy.astro.gravitation.legendre_normalization_factor` function.

 Parameters
 ----------
 unnormalized_cosine_coefficients : numpy.ndarray
     Matrix for which entry :math:`(i,j)` contains the unnormalized cosine coefficient :math:`C_{lm}`
 unnormalized_sine_coefficients : numpy.ndarray
     Matrix for which entry :math:`(i,j)` contains the unnormalized sine coefficient :math:`S_{lm}`
 Returns
 -------
 tuple[numpy.ndarray, numpy.ndarray]
     Tuple of two matrices, containing the normalized coefficients :math:`\bar{C}_{lm}` (first) and :math:`\bar{S}_{lm}` (second)






     )doc" );

    m.def( "unnormalize_spherical_harmonic_coefficients",
           py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd& >( &tbm::convertGeodesyNormalizedToUnnormalizedCoefficients ),
           py::arg( "normalized_cosine_coefficients" ),
           py::arg( "normalized_sine_coefficients" ),
           R"doc(

 Function to unnormalize spherical harmonic coefficients

 Function to unnormalize spherical harmonic coefficients, using the equations provided in the :func:`~tudatpy.gravitation.astro.legendre_normalization_factor` function.

 Parameters
 ----------
 normalized_cosine_coefficients : numpy.ndarray
     Matrix for which entry :math:`(i,j)` contains the normalized cosine coefficient :math:`\bar{C}_{lm}`
 normalized_sine_coefficients : numpy.ndarray
     Matrix for which entry :math:`(i,j)` contains the normalized sine coefficient :math:`\bar{S}_{lm}`
 Returns
 -------
 tuple[numpy.ndarray, numpy.ndarray]
     Tuple of two matrices, containing the unnormalized coefficients :math:`{C}_{lm}` (first) and :math:`{S}_{lm}` (second)






     )doc" );

    m.def( "spherical_harmonic_coefficients_from_inertia",
           tg::getDegreeTwoSphericalHarmonicCoefficientsPy,
           py::arg( "inertia_tensor" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "reference_radius" ),
           py::arg( "output_normalized_coefficients" ) = true,
           R"doc(

 Function to compute degree-two spherical harmonic coefficients from an inertia tensor

 Function to compute degree-two spherical harmonic coefficients :math:`C_{20}`, :math:`C_{21}`, :math:`C_{22}`, :math:`S_{21}`, :math:`S_{22}` and from an inertia tensor :math:`\mathbf{I}`, according to the following relation"

 .. math::
     \mathbf{I}=M R^2\left(\left(\begin{array}{ccc} \frac{C_{20}}{3}-2 C_{22} & -2 S_{22} & -C_{21} \\ -2 S_{22} & \frac{C_{20}}{3}+2 C_{22} & -S_{21} \\ -C_{21} & -S_{21} & -\frac{2 C_{20}}{3} \end{array}\right)+\bar{I} \mathbf{1}_{3 \times 3}\right)

 with :math:`M` the mass of the body, and :math:`R` the reference radius of the spherical harmonic coefficients. The term :math:`\bar{I}` is the mean moment of inertia. The spherical harmonic
 coefficients in the above matrix are unnormalized.


 Parameters
 ----------
 inertia tensor : numpy.ndarray[numpy.float64[3, 3]]
     Full inertia tensor :math:`\mathbf{I}` of the body for which spherical harmonic coefficients are to be computed.
 gravitational_parameter : float
     Gravitational parameter :math:`\mu` of the body for which spherical harmonic coefficients are to be computed.
 reference_radius : float
     Reference radius w.r.t. which spherical harmonic coefficients are to be computed.
 output_normalized_coefficients : bool, default=True
     Boolean denoting whether the coefficients computed are normalized (if true) or unnormalized (if false)
 Returns
 -------
 tuple[numpy.ndarray, numpy.ndarray]
     Tuple of two matrices, containing the spherical harmonic coefficients :math:`{C}_{lm}` (first) and :math:`{S}_{lm}` (second) up to degree and order 2.
     The degree-two coefficients are computed from the inertia tensor, the degree-one coefficients are set to zero (and :math:`C_{00}=0`)







     )doc" );

    m.def( "transform_spherical_harmonic_coefficients_with_wigner_d",
           &tg::transformSphericalHarmonicCoefficientsWithWignerDPy,
           py::arg( "cosine_coefficients" ),
           py::arg( "sine_coefficients" ),
           py::arg( "wigner_matrices" ),
           py::arg( "are_coefficients_normalized" ) = true,
           R"doc(

 Function to rotate spherical harmonic coefficients with precomputed Wigner D-matrices

 Function to rotate real-valued spherical harmonic coefficients using precomputed Wigner D-matrices.
 The transformation follows the Wigner-D coefficient rotation formalism used by :cite:t:`boue2017`; see
 :class:`~tudatpy.astro.gravitation.WignerDMatricesCache` for the algorithm used to compute the Wigner D-matrices.
 For each degree :math:`l`, input coefficients :math:`C_{lk},S_{lk}` are mapped to
 :math:`C'_{lm},S'_{lm}` using the Wigner-D entries :math:`D^l_{m,k}` as

 .. math::

    \tilde{C}_{lm} =
    \eta_m\Re(D^l_{m,0})C_{l0}
    +\frac{1}{2}\sum_{k=1}^{l}\eta_{mk}
    \left[
    \left((-1)^k\Re(D^l_{m,k})+\Re(D^l_{m,-k})\right)C_{lk}
    +\left((-1)^k\Im(D^l_{m,k})-\Im(D^l_{m,-k})\right)S_{lk}
    \right]

 .. math::

    \tilde{S}_{lm} =
    \eta_m\Im(D^l_{m,0})C_{l0}
    +\frac{1}{2}\sum_{k=1}^{l}\eta_{mk}
    \left[
    \left((-1)^k\Im(D^l_{m,k})+\Im(D^l_{m,-k})\right)C_{lk}
    +\left(-(-1)^k\Re(D^l_{m,k})+\Re(D^l_{m,-k})\right)S_{lk}
    \right]

 .. math::

    C'_{lm}=(-1)^m(2-\delta_{0m})\tilde{C}_{lm},\qquad
    S'_{lm}=(-1)^{m+1}(2-\delta_{0m})\tilde{S}_{lm}.

 For normalized coefficients, :math:`\eta_m=1` when :math:`m=0`, :math:`\eta_m=1/\sqrt{2}` when
 :math:`m>0`, and :math:`\eta_{mk}=\sqrt{2}\eta_m`. For unnormalized coefficients,
 :math:`\eta_m=\sqrt{(l-m)!/(l+m)!}` and
 :math:`\eta_{mk}=\sqrt{(l+k)!/(l-k)!}\eta_m`.


 Parameters
 ----------
 cosine_coefficients : numpy.ndarray
     Matrix of cosine spherical harmonic coefficients.
 sine_coefficients : numpy.ndarray
     Matrix of sine spherical harmonic coefficients.
 wigner_matrices : list[numpy.ndarray]
     Complex Wigner D-matrices, one matrix per spherical harmonic degree.
 are_coefficients_normalized : bool, default=True
     Boolean denoting whether the input coefficients are fully normalized.
 Returns
 -------
 tuple[numpy.ndarray, numpy.ndarray]
     Tuple containing the rotated cosine and sine coefficient matrices.





     )doc" );
}

}  // namespace gravitation
}  // namespace astro
}  // namespace tudatpy
