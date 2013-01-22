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
 *      110623    K. Kumar          File created.
 *      110701    K. Kumar          Added Mercury, Saturn, Neptune.
 *      120327    D. Dirkx          Moved setting of predefined field to constructor.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_CENTRAL_GRAVITY_FIELD_H
#define TUDAT_CENTRAL_GRAVITY_FIELD_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
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

//! Typedef for shared-pointer to CentralGravityField object.
typedef boost::shared_ptr< CentralGravityField > CentralGravityFieldPointer;

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_CENTRAL_GRAVITY_FIELD_H
