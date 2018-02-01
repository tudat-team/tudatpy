/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_GRAVITY_FIELD_MODEL_H
#define TUDAT_GRAVITY_FIELD_MODEL_H

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

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
