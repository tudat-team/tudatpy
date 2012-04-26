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
 *      110623    K. Kumar          File created.
 *      110701    K. Kumar          Added Mercury, Saturn, Neptune.
 *      120327    D. Dirkx          Moved setting of predefined field to constructor.
 *
 *    References
 *
 */

#ifndef TUDAT_CENTRAL_GRAVITY_FIELD_H
#define TUDAT_CENTRAL_GRAVITY_FIELD_H

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{

//! Central gravity field class.
/*!
 * Class that defines a central gravity field.
 */
class CentralGravityField : public SphericalHarmonicsGravityField
{
public:

    //! Bodies with predefined central gravity fields.
    /*!
     * Bodies with predefined central gravity fields.
     */
    enum BodiesWithPredefinedCentralGravityFields
    {
        sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune
    };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CentralGravityField( BodiesWithPredefinedCentralGravityFields
                         bodyWithPredefinedCentralGravityField )
    {
        setPredefinedCentralGravityFieldSettings( bodyWithPredefinedCentralGravityField );
    }

protected:

private:

    //! Set predefined central gravity field settings.
    /*!
     * Sets predefined central gravity field settings.
     * \param bodyWithPredefinedCentralGravityField Body with predefined
     *          central gravity field.
     */
    void setPredefinedCentralGravityFieldSettings(
        BodiesWithPredefinedCentralGravityFields bodyWithPredefinedCentralGravityField );
};

} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_CENTRAL_GRAVITY_FIELD_H
