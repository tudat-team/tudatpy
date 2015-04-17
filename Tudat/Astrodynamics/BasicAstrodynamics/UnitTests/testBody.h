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
 *      120821    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_TEST_BODY_H
#define TUDAT_TEST_BODY_H

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{

//! Test body class.
/*!
 * This class serves as an example of how a container can be constructed that stored state and
 * time information, which can be used in conjunction with acceleration models, the Cartesian state
 * derivative model, and the composite state derivative model. It should be noted that this class
 * should NOT be used as is in the rest of the library, and is only used in conjunction with unit
 * tests. Classes of this nature should go in user applications, since they are typically
 * application-specific, hence not general enough to be added to the Tudat libraries.
 * \tparam SpatialDimensions Spatial dimensions of position and velocity vectors. The state vector
 *          is taken to be twice the number of spatial dimensions. All vectors are represented by
 *          Eigen::Matrix-types.
 * \tparam DataType The data type of all of the internally stored state-related vectors and time.
 */
template< int SpatialDimensions = 3, typename DataType = double >
class TestBody
{
private:

    //! Typedef for time.
    typedef DataType TimeDataType;

    //! Typedef for position vectors.
    typedef Eigen::Matrix< DataType, SpatialDimensions, 1 > PositionVectorType;

    //! Typedef for velocity vectors.
    typedef Eigen::Matrix< DataType, SpatialDimensions, 1 > VelocityVectorType;

    //! Typedef for state vectors.
    typedef Eigen::Matrix< DataType, 2 * SpatialDimensions, 1 > StateVectorType;

public:

    //! Constructor taking a state and a time.
    /*!
     * Constructor taking an input state and time. The input state is used internally to
     * set the current position (taken as a segment of the input state given by the indices
     * (0, SpatialDimensions)) and the current velocity (taken as a segment of the input state
     * given by the indices (SpatialDimensions, SpatialDimensions).
     * \param aState An input state vector.
     * \param aTime An input time.
     */
    TestBody( const StateVectorType& aState, const TimeDataType aTime )
        : currentPosition( aState.segment( 0, SpatialDimensions ) ),
          currentVelocity( aState.segment( SpatialDimensions, SpatialDimensions ) ),
          currentTime( aTime )
    { }

    //! Set current time and state.
    /*!
     * Sets the current time, position and current velocity internally based on the input
     * arguments. The current position is taken as a segment of the input state given by the
     * indices (0, SpatialDimensions)), and the current velocity is taken as a segment of the input
     * state given by the indices (SpatialDimensions, SpatialDimensions).
     * \param aTime An input time.
     * \param aState An input state vector.
     */
    void setCurrentTimeAndState( const TimeDataType aTime, const StateVectorType& aState )
    {
        currentTime = aTime;
        currentPosition = aState.segment( 0, SpatialDimensions );
        currentVelocity = aState.segment( SpatialDimensions, SpatialDimensions );
    }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    PositionVectorType getCurrentPosition( ) { return currentPosition; }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    VelocityVectorType getCurrentVelocity( ) { return currentVelocity; }

    //! Get current time.
    /*!
     * Returns the internally stored current tome.
     * \return Current time.
     */
    TimeDataType getCurrentTime( ) { return currentTime; }

protected:

private:

    //! Current position.
    PositionVectorType currentPosition;

    //! Current position.
    VelocityVectorType currentVelocity;

    //! Current time.
    TimeDataType currentTime;
};

//! Test body class.
/*!
 * This class is a template specialization of the general TestBody class for scalar
 * positions/velocities (i.e. doubles, floats etc. instead if Eigen-types).
 * \tparam DataType The data type of all of the internally stored state-related vectors and time.
 */
template< typename DataType >
class TestBody< 1, DataType >
{
private:

    //! Typedef for time.
    typedef DataType TimeDataType;

    //! Typedef for position vectors.
    typedef DataType PositionVectorType;

    //! Typedef for velocity vectors.
    typedef DataType VelocityVectorType;

    //! Typedef for state vectors.
    typedef Eigen::Matrix< DataType, 2, 1 > StateVectorType;

public:

    //! Constructor taking a state and a time.
    /*!
     * Constructor taking an input state and time. The input state is used internally to set the
     * current position (taken as first element of state) and the current velocity (taken as second
     * element of state).
     * \param aState An input state vector.
     * \param aTime An input time.
     */
    TestBody( const StateVectorType& aState, const TimeDataType aTime )
        : currentPosition( aState( 0 ) ),
          currentVelocity( aState( 1 ) ),
          currentTime( aTime )
    { }

    //! Set current time and state.
    /*!
     * Sets the current time, position and current velocity internally based on the input
     * arguments. The current position is taken as the first element of the input state, and the
     * current velocity is taken as the second element of the input state.
     * \param aTime An input time.
     * \param aState An input state vector.
     */
    void setCurrentTimeAndState( const TimeDataType aTime, const StateVectorType& aState )
    {
        currentTime = aTime;
        currentPosition = aState( 0 );
        currentVelocity = aState( 1 );
    }

    //! Get current position.
    PositionVectorType getCurrentPosition( ) { return currentPosition; }

    //! Get current velocity.
    VelocityVectorType getCurrentVelocity( ) { return currentVelocity; }

    //! Get current time.
    TimeDataType getCurrentTime( ) { return currentTime; }

protected:

private:

    //! Current position.
    PositionVectorType currentPosition;

    //! Current position.
    VelocityVectorType currentVelocity;

    //! Current time.
    TimeDataType currentTime;
};

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_BODY_H
