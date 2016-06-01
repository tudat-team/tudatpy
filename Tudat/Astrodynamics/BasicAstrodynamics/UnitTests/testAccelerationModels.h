/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120820    K. Kumar          Created based off of test acceleration models originally used
 *                                  in unitTestAccelerationModel.cpp.
 *      120821    K. Kumar          Updated AnotherDerivedAccelerationModel class; completed
 *                                  Doxygen documentation.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_TEST_ACCELERATION_MODELS_H
#define TUDAT_TEST_ACCELERATION_MODELS_H

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

namespace tudat
{
namespace unit_tests
{

//! Acceleration model derived class.
/*!
 * This class serves as an example of a derived class from the AccelerationModel base class. This
 * class contains an implementation that computes acceleration based on a given position and time.
 * The context of this class is somewhat specific to make it more understandable, though derived
 * classes of the AccelerationModel base class can be completely general. It should be noted that
 * this class should NOT be used as is in the rest of the library, and is only used in conjunction
 * with unit tests.
 * \tparam AccelerationDataType Data type used to represent accelerations
 *          (default=Eigen::Vector3d).
 * \tparam PositionDataType Data type used to represent positions (default=Eigen::Vector3d).
 * \tparam TimeDataType Data type used to represent time (default=double).
*/
template< typename AccelerationDataType = Eigen::Vector3d,
          typename PositionDataType = Eigen::Vector3d, typename TimeDataType = double >
class DerivedAccelerationModel
        : public basic_astrodynamics::AccelerationModel< AccelerationDataType >
{
private:

    //! Typedef for a pointer to a function that returns a position.
    typedef boost::function< PositionDataType( ) > PositionReturningFunction;

    //! Typedef for a pointer to a function that returns a time.
    typedef boost::function< TimeDataType( ) > TimeReturningFunction;

public:

    //! Default constructor.
    /*!
     * Default constructor that takes function-pointers as input. These function-pointers point to
     * functions that return a position and a time. Internally, the updateMembers() function also
     * gets called to ensure that all members are up-to-date after construction.
     * \param aPositionFunction A function-pointer that points to a function returning a position.
     * \param aTimeFunction A function-pointer that points to a function returning a time.
     */
    DerivedAccelerationModel( const PositionReturningFunction aPositionFunction,
                              const TimeReturningFunction aTimeFunction )
        : getPosition( aPositionFunction ),
          getTime( aTimeFunction )
    {
        updateMembers( );
    }

    //! Get acceleration.
    /*!
     * Returns acceleration. In this case, this functions returns an acceleration that is dependent
     * on the internally stored position and time members. This function merely serves as an
     * example, rather than representing real dynamics.
     * \return Computed acceleration.
     */
    AccelerationDataType getAcceleration( ) { return position / ( time * time ); }

    //! Update member variables used by the acceleration model.
    /*!
     * Updates member variables used by the acceleration model. In this case, the internally stored
     * position and time are updated by calling the function-pointers passed to the constructor.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN ) { position = getPosition( ); time = getTime( ); }

protected:

private:

    //! Pointer to position-returning function.
    /*!
     * Function-pointer that points to a function passed to the constructor that returns a
     * position. This function-pointer is used by the updateMembers() function to update the
     * internally stored position.
     */
    PositionReturningFunction getPosition;

    //! Pointer to time-returning function.
    /*!
     * Function-pointer that points to a function passed to the constructor that returns a time.
     * This function-pointer is used by the updateMembers() function to update the internally
     * stored time.
     */
    TimeReturningFunction getTime;

    //! Position.
    /*!
     * Internally stored current position, used by the getAcceleration() function.
     */
    PositionDataType position;

    //! Time.
    /*!
     * Internally stored current time, used by the getAcceleration() function.
     */
    TimeDataType time;
};

//! Another acceleration model derived class.
/*!
 * This class serves as an example of another derived class from the AccelerationModel base class.
 * This class contains an implementation that computes acceleration based on a given position,
 * velocity, and time. The context of this class is somewhat specific to make it more
 * understandable, though derived classes of the AccelerationModel base class can be completely
 * general. It should be noted that this class should NOT be used as is in the rest of the library,
 * and is only used in conjunction with unit tests.
 * \tparam AccelerationDataType Data type used to represent accelerations
 *          (default=Eigen::Vector3d).
 * \tparam PositionDataType Data type used to represent positions (default=Eigen::Vector3d).
 * \tparam VelocityDataType Data type used to represent velocities (default=Eigen::Vector3d).
 * \tparam TimeDataType Data type used to represent time (default=double).
*/
template< typename AccelerationDataType = Eigen::Vector3d,
          typename PositionDataType = Eigen::Vector3d, typename VelocityDataType = Eigen::Vector3d,
          typename TimeDataType = double >
class AnotherDerivedAccelerationModel
        : public basic_astrodynamics::AccelerationModel< AccelerationDataType >
{
public:

    //! Typedef for a pointer to a function that returns a position.
    typedef boost::function< PositionDataType( ) > PositionReturningFunction;

    //! Typedef for a pointer to a function that returns a position.
    typedef boost::function< VelocityDataType( ) > VelocityReturningFunction;

    //! Typedef for a pointer to a function that returns a time.
    typedef boost::function< TimeDataType( ) > TimeReturningFunction;

    //! Default constructor.
    /*!
     * Default constructor that takes function-pointers as input. These function-pointers point to
     * functions that return a position, a velocity, and a time. Internally, the updateMembers()
     * function also gets called to ensure that all members are up-to-date after construction.
     * \param aPositionFunction A function-pointer that points to a function returning a position.
     * \param aVelocityFunction A function-pointer that points to a function returning a velocity.
     * \param aTimeFunction A function-pointer that points to a function returning a time.
     */
    AnotherDerivedAccelerationModel( const PositionReturningFunction aPositionFunction,
                                     const VelocityReturningFunction aVelocityFunction,
                                     const TimeReturningFunction aTimeFunction )
        : getPosition( aPositionFunction ),
          getVelocity( aVelocityFunction ),
          getTime( aTimeFunction )
    {
        updateMembers( );
    }

    //! Get acceleration.
    /*!
     * Returns acceleration. In this case, this functions returns an acceleration that is dependent
     * on the internally stored position, velocity, and time members. This function merely serves
     * as an example, rather than representing real dynamics.
     * \return Computed acceleration.
     */
    AccelerationDataType getAcceleration( )
    {
        return 0.5 * position / ( 3.2 * ( time + 3.4 ) * time ) + velocity / time;
    }

    //! Update member variables used by the acceleration model.
    /*!
     * Updates member variables used by the acceleration model. In this case, the internally stored
     * position, velocity, and time are updated by calling the function-pointers passed to the
     * constructor.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        position = getPosition( );
        velocity = getVelocity( );
        time = getTime( );
    }

protected:

private:

    //! Pointer to position-returning function.
    /*!
     * Function-pointer that points to a function passed to the constructor that returns a
     * position. This function-pointer is used by the updateMembers() function to update the
     * internally stored position.
     */
    PositionReturningFunction getPosition;

    //! Pointer to velocity-returning function.
    /*!
     * Function-pointer that points to a function passed to the constructor that returns a
     * velocity. This function-pointer is used by the updateMembers() function to update the
     * internally stored velocity.
     */
    VelocityReturningFunction getVelocity;

    //! Pointer to time-returning function.
    /*!
     * Function-pointer that points to a function passed to the constructor that returns a time.
     * This function-pointer is used by the updateMembers() function to update the internally
     * stored time.
     */
    TimeReturningFunction getTime;

    //! Position.
    /*!
     * Internally stored current position, used by the getAcceleration() function.
     */
    PositionDataType position;

    //! Velocity.
    /*!
     * Internally stored current velocity, used by the getAcceleration() function.
     */
    VelocityDataType velocity;

    //! Time.
    /*!
     * Internally stored current time, used by the getAcceleration() function.
     */
    TimeDataType time;
};

} // namespace tudat
} // namespace unit_tests

#endif // TUDAT_TEST_ACCELERATION_MODELS_H
