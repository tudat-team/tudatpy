/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      120827    K. Kumar          File created.
 *      121105    K. Kumar          Simplified base class definition.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITATIONAL_ACCELERATION_MODEL_BASE_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITATIONAL_ACCELERATION_MODEL_BASE_H

#include <boost/function.hpp>

namespace tudat
{
namespace gravitation
{

//! Template base class for spherical harmonics gravitational acceleration models.
/*!
 * This template class serves as the base class for the
 * SphericalHarmonicsGravitationalAccelerationModel, CentralGravitationalAccelerationModel,
 * CentralJ2GravitationalAccelerationModel, CentralJ2J3GravitationalAccelerationModel, and
 * CentralJ2J3J4GravitationalAccelerationModel classes.
 * \tparam DataType Data type.
 * \tparam PositionType Data type for position vector.
 */
template< typename DataType, typename PositionType >
class SphericalHarmonicsGravitationalAccelerationModelBase
{
protected:

    //! Typedef for a (DataType)-returning function.
    typedef boost::function< DataType( ) > DataTypeReturningFunction;

    //! Typedef for a position-returning function.
    typedef boost::function< PositionType( ) > PositionReturningFunction;

    //! Typedef for acceleration type.
    typedef PositionType AccelerationType;

public:

    //! Default constructor taking position of body subject to acceleration, variable
    //! gravitational parameter, and position of body exerting accelearation.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, a pointer to a function returning the gravitational parameter of
     * the body exerting the acceleration, and a pointer to a function returning the position of
     * the body exerting the  gravitational acceleration (typically the central body). The
     * constructor also updates all the internal members.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param gravitationalParameterFunction Pointer to function returning gravitational parameter
     *          of body exerting gravitational acceleration.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration.
     */
    SphericalHarmonicsGravitationalAccelerationModelBase(
            const PositionReturningFunction positionOfBodySubjectToAccelerationFunction,
            const DataTypeReturningFunction gravitationalParameterFunction,
            const PositionReturningFunction positionOfBodyExertingAccelerationFunction )
        : getPositionOfBodySubjectToAcceleration( positionOfBodySubjectToAccelerationFunction ),
          getGravitationalParameter( gravitationalParameterFunction ),
          getPositionOfBodyExertingAcceleration( positionOfBodyExertingAccelerationFunction )
    { }

    //! Virtual destructor.
    /*!
     * Base class virtual destructor.
     */
    virtual ~SphericalHarmonicsGravitationalAccelerationModelBase( ){ }

    //! Update base class members.
    /*!
     * Updates members of base class. In this case, the current position of the body subject to and
     * exerting the gravitational acceleration, as well as the gravitational parameter are updated
     * by calling the functions provided through the constructor.
     * \return True; this should be modified to return a flag indicating if the update was
     *          successful.
     */
    bool updateBaseMembers( );

protected:

    //! Position of body subject to acceleration.
    /*!
     * Current position vector of body subject to gravitational acceleration in inertial frame.
     */
    PositionType positionOfBodySubjectToAcceleration;

    //! Gravitational parameter.
    /*!
     * Current gravitational parameter of body exerting acceleration.
     */
    DataType gravitationalParameter;

    //! Position of body exerting acceleration.
    /*!
     * Current position vector of body exerting gravitational acceleration in inertial frame.
     */
    PositionType positionOfBodyExertingAcceleration;

    //! Pointer to function returning position of body subject to acceleration.
    /*!
     * Pointer to function that returns the current position of the body subject to the
     * gravitational acceleration.
     */
    const PositionReturningFunction getPositionOfBodySubjectToAcceleration;

    //! Pointer to function returning gravitatational parameter.
    /*!
     * Pointer to function that returns the current value of the gravitational parameter of the
     * body exerting the gravitational acceleration.
     */
    const DataTypeReturningFunction getGravitationalParameter;

    //! Pointer to function returning position of body exerting acceleration.
    /*!
     * Pointer to function that returns the current position of the body exerting the
     * gravitational acceleration.
     */
    const PositionReturningFunction getPositionOfBodyExertingAcceleration;

private:
};

// Template class source.
// The code given below is effectively the ".cpp file" for the template class definition, so you
// only need to look at the code below if you are interested in the source implementation.

//! Update base class members.
template< typename DataType, typename PositionType >
bool SphericalHarmonicsGravitationalAccelerationModelBase<
DataType, PositionType >::updateBaseMembers( )
{
    this->positionOfBodySubjectToAcceleration
            = this->getPositionOfBodySubjectToAcceleration( );
    this->gravitationalParameter = this->getGravitationalParameter( );
    this->positionOfBodyExertingAcceleration
            = this->getPositionOfBodyExertingAcceleration( );
    return true;
}

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITATIONAL_ACCELERATION_MODEL_BASE_H
