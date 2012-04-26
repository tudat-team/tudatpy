/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      101110    K. Kumar          File created.
 *      101116    K. Kumar          Changed filename and class name.
 *      101215    K. Kumar          Added virtual functions and missing Doxygen comments.
 *      101216    K. Kumar          Added set/get functions for origin.
 *      110107    K. Kumar          Removed reference radius get/set virtual functions.
 *      110113    K. Kumar          Updated arguments of get-functions with const.
 *      110202    K. Kumar          Updated code to use the CartesianPositionElements class.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *
 *    References
 *
 */

#ifndef TUDAT_GRAVITY_FIELD_MODEL_H
#define TUDAT_GRAVITY_FIELD_MODEL_H

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{

//! GravityFieldModel class.
/*!
 * Gravity field model class included in Tudat.
 */
class GravityFieldModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GravityFieldModel( ) : gravitationalParameter_( TUDAT_NAN ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~GravityFieldModel( ) { }

    //! Set the gravitational parameter.
    /*!
     * Define the gravitational parameter in meter^3 per second^2.
     * \param gravitationalParameter
     */
    void setGravitationalParameter( const double gravitationalParameter )
    {
        gravitationalParameter_ = gravitationalParameter;
    }

    //! Set origin of gravity field.
    /*!
     * Set origin of gravity field.
     * \param positionOfOrigin Position of origin
     */
    void setOrigin( const Eigen::VectorXd& positionOfOrigin )
    {
        positionOfOrigin_ = positionOfOrigin;
    }

    //! Get the gravitational parameter.
    /*!
     * Return the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    double getGravitationalParameter( ) { return gravitationalParameter_; }

    //! Get origin of gravity field.
    /*!
     * Get origin of gravity field.
     * \return Position of origin
     */
    Eigen::VectorXd getOrigin( ) { return positionOfOrigin_; }

    //! Get the potential.
    /*!
     * Returns the potential for the gravity field selected.
     * \param position Position at which potential is to be determined
     * \return Potential.
     */
    virtual double getPotential( const Eigen::Vector3d& position ) = 0;

    //! Get the gradient of the potential.
    /*!
     * Returns the gradient of the potential for the gravity field selected.
     * \param position Position at which gradient of potential is to be determined
     * \return Gradient of potential.
     */
    virtual Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& position ) = 0;

    //! Get gradient tensor of the potential.
    /*!
     * Returns the gradient tensor of the potential for the gravity field
     * selected.
     * \param position Position at which gradient tensor of potential is to be determined
     * \return Gradient tensor of potential.
     */
    virtual Eigen::Matrix3d getGradientTensorOfPotential( const Eigen::Vector3d& position ) = 0;

protected:

    //! Gravitational parameter.
    /*!
     * The gravitational parameter in meter^3 per second^2.
     */
    double gravitationalParameter_;

    //! Origin of gravity field.
    /*!
     * Origin of gravity field given in Cartesian Elements as an Eigen::VectorXd
     */
    Eigen::VectorXd positionOfOrigin_;

    //! Relative position of
    /*!
     * Relative position given in Cartesian Elements as an Eigen::VectorXd
     */
    Eigen::VectorXd relativePosition_;

private:
};

} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITY_FIELD_MODEL_H
