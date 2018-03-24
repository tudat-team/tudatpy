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

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_ACCELERATION_MODELS_H
