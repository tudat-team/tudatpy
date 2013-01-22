/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      101116    K. Kumar          File created.
 *      101117    K. Kumar          Added getPotential( ) and getGradientOfPotential functions.
 *      101214    K. Kumar          Updated getGradientOfPotential( ) and getLaplacianOfPotential( ).
 *      110106    K. Kumar          Added order and degree of expansion variables and set/get
 *                                  functions.
 *      110202    K. Kumar          Updated code to make use of the CartesianPositionElements
 *                                  class.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *      110805    K. Kumar          Added predefined functionality with WGS-72 and WGS-84
 *                                  predefined Earth gravity fields.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References
 *      Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. Revisiting Spacetrack Report #3:
 *          Rev 1, Proceedings of the AIAA/AAS Astrodynamics Specialist Conference. Keystone, CO,
 *          2006.
 *
 *    Notes
 *      The coefficients J2, J3, and J4 have been hardcoded for now, but in future these variables
 *      should be removed and the data should be stored in data files.
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H

#include <iostream>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"

namespace tudat
{
namespace gravitation
{

//! SphericalHarmonicsGravityField class.
/*!
 * Spherical harmonics gravity field model class included in Tudat.
 */
class SphericalHarmonicsGravityField : public GravityFieldModel
{
public:

    //! Bodies with predefined spherical harmonics gravity fields.
    /*!
     * Bodies with predefined spherical harmonics gravity fields.
     */
    enum BodiesWithPredefinedSphericalHarmonicsGravityFields
    {
        earthWorldGeodeticSystem72,
        earthWorldGeodeticSystem84
    };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SphericalHarmonicsGravityField( const unsigned int degreeOfExpansion = 0,
                                    const unsigned int orderOfExpansion = 0,
                                    const double referenceRadius = 0.0 )
        : degreeOfExpansion_( degreeOfExpansion ),
          orderOfExpansion_( orderOfExpansion ),
          referenceRadius_( referenceRadius ),
          j2SphericalHarmonicsGravityFieldCoefficient_( TUDAT_NAN ),
          j3SphericalHarmonicsGravityFieldCoefficient_( TUDAT_NAN ),
          j4SphericalHarmonicsGravityFieldCoefficient_( TUDAT_NAN )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~SphericalHarmonicsGravityField( ) { }

    //! Set predefined spherical harmonics gravity field settings.
    /*!
     * Sets predefined spherical harmonics gravity field settings.
     * \param bodyWithPredefinedSphericalHarmonicsGravityField Body with
     *          predefined spherical harmonics gravity field.
     */
    void setPredefinedSphericalHarmonicsGravityFieldSettings(
        BodiesWithPredefinedSphericalHarmonicsGravityFields
        bodyWithPredefinedSphericalHarmonicsGravityField );

    //! Get the reference radius.
    /*!
     * Returns the reference radius used for the spherical harmonics expansion in meters.
     * \return Reference radius.
     */
    double getReferenceRadius( ) { return referenceRadius_; }

    //! Get degree of spherical harmonics gravity field expansion.
    /*!
     * Returns the degree of the spherical harmonics gravity field expansion.
     * \return Degree of spherical harmonics expansion.
     */
    double getDegreeOfExpansion( ) { return degreeOfExpansion_; }

    //! Get order of spherical harmonics gravity field expansion.
    /*!
     * Returns the order of the spherical harmonics gravity field expansion.
     * \return Order of spherical harmonics expansion.
     */
    double getOrderOfExpansion( ) { return orderOfExpansion_; }

    //! Get the gravitational potential.
    /*!
     * Returns the value of the gravitational potential, expressed in spherical harmonics for the
     * given position.
     * \param position Position at which potential is to be determined
     * \return Gravitational potential.
     */
    double getPotential( const Eigen::Vector3d& position )
    {
        relativePosition_ = position - positionOfOrigin_;
        return gravitationalParameter_ / relativePosition_.norm( );
    }

    //! Get the gradient of the gravitational potential.
    /*!
     * Returns the value of the gradient of the gravitational potential, expressed in spherical
     * harmonics for the given position.
     * \param position Position at which gradient of potential is to be determined
     * \return Gradient of gravitational potential.
     */
    Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& position )
    {
        relativePosition_ = position - positionOfOrigin_;
        return -gravitationalParameter_ * relativePosition_
                / pow( relativePosition_.norm( ), 3.0 );
    }

    //! Get gradient tensor of the gravitational potential.
    /*!
     * Returns the value of the gradient tensor of the gravitational potential expressed in
     * spherical harmonics for the given position.
     * \param position Position at which gradient tensor of potential is to be determined
     * \return Gradient tensor of gravitational potential.
     */
    Eigen::Matrix3d getGradientTensorOfPotential( const Eigen::Vector3d& position );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param sphericalHarmonicsGravityField Spherical harmonics gravity field.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     SphericalHarmonicsGravityField&
                                     sphericalHarmonicsGravityField );

protected:

    //! Degree of spherical harmonics expansion.
    /*!
     * Degree of spherical harmonics expansion.
     */
    unsigned int degreeOfExpansion_;

    //! Order of spherical harmonics expansion.
    /*!
     * Order of spherical harmonics expansion.
     */
    unsigned int orderOfExpansion_;

    //! Reference radius.
    /*!
     * The reference radius used for the spherical harmonics expansion in
     * meters.
     */
    double referenceRadius_;

private:

    //! J2 spherical harmonics gravity field coefficient.
    /*!
     * J2 spherical harmonics gravity field coefficient.
     */
    double j2SphericalHarmonicsGravityFieldCoefficient_;

    //! J3 spherical harmonics gravity field coefficient.
    /*!
     * J3 spherical harmonics gravity field coefficient.
     */
    double j3SphericalHarmonicsGravityFieldCoefficient_;

    //! J4 spherical harmonics gravity field coefficient.
    /*!
     * J4 spherical harmonics gravity field coefficient.
     */
    double j4SphericalHarmonicsGravityFieldCoefficient_;
};

//! Typedef for shared-pointer to SphericalHarmonicsGravityField object.
typedef boost::shared_ptr< SphericalHarmonicsGravityField > SphericalHarmonicsGravityFieldPointer;

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITY_FIELD_H
