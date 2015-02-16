/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      101110    K. Kumar          File created.
 *      101116    K. Kumar          Changed filename and class name.
 *      101215    K. Kumar          Added virtual functions and missing Doxygen comments.
 *      101216    K. Kumar          Added set/get functions for origin.
 *      110107    K. Kumar          Removed reference radius get/set virtual functions.
 *      110113    K. Kumar          Updated arguments of get-functions with const.
 *      110202    K. Kumar          Updated code to use the CartesianPositionElements class.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *      120502    K. Kumar          Added missing constructor initialization of position vectors.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      15YOLO    D. Dirkx          Changed some stuff.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_GRAVITY_FIELD_MODEL_H
#define TUDAT_GRAVITY_FIELD_MODEL_H

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

namespace tudat
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
     * \param gravitationalParameter Gravitational parameter associated with gravity field
     */
    GravityFieldModel( const double gravitationalParameter ):
        gravitationalParameter_( gravitationalParameter )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~GravityFieldModel( ) { }

    //! Set the gravitational parameter.
    /*!
     * Define the gravitational parameter in meter^3 per second^2.
     * \param gravitationalParameter New gravitational parameter associated with gravity field.
     */
    void resetGravitationalParameter( const double gravitationalParameter )
    {
        gravitationalParameter_ = gravitationalParameter;
    }

    //! Get the gravitational parameter.
    /*!
     * Return the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    double getGravitationalParameter( )
    {
        return gravitationalParameter_;
    }


    //! Get the gravitational potential at given body-fixed position.
    /*!
     * Return the gravitational potential at given body-fixed position.
     * \param bodyFixedPosition Position at which the gravitational potential is to be evaluated.
     * \return Gravitational potential.
     */
    double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        return gravitationalParameter_ / bodyFixedPosition.norm( );
    }

    //! Get the gradient of the potential.
    /*!
     * Returns the gradient of the potential for the gravity field selected.
     * \param bodyFixedPosition Position at which gradient of potential is to be determined.
     * \return Gradient of potential.
     */
    virtual Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        return computeGravitationalAcceleration( bodyFixedPosition, gravitationalParameter_ );
    }

protected:

    //! Gravitational parameter.
    /*!
     * The gravitational parameter in meter^3 per second^2.
     */
    double gravitationalParameter_;


private:
};

//! List of bodies for which predefined central gravity fields may be created through the
//! getPredefinedCentralGravityField function.
enum BodiesWithPredefinedCentralGravityFields
{
    sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune
};

//! Typedef for shared-pointer to GravityFieldModel object.
typedef boost::shared_ptr< GravityFieldModel > GravityFieldModelPointer;


//! Function to create a central gravity field model of one of the planets, moon or sun.
/*!
 *  Function to create a central gravity field model of one of the planets, moon or sun.
 *  \param bodyWithPredefinedCentralGravityField Identified determining for which body a
 *  gravity field is to be created.
 *  \return Central gravity field model of requested body.
 */
boost::shared_ptr< GravityFieldModel > getPredefinedCentralGravityField(
    BodiesWithPredefinedCentralGravityFields bodyWithPredefinedCentralGravityField );

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_GRAVITY_FIELD_MODEL_H
