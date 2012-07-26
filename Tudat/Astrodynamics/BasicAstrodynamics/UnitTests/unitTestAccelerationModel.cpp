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
 *      120720    A. Ronse          First creation of the unit test.
 *      120724    K. Kumar          Addition of extensive comments and tests for
 *                                  updateAndGetAcceleration functions.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

namespace tudat
{
namespace unit_tests
{

using tudat::basic_astrodynamics::acceleration_models::AccelerationModel;
using tudat::basic_astrodynamics::acceleration_models::updateAndGetAcceleration;

//! A test value of double-type used to demonstrate the updateMembers() function.
static double testDouble = 3.2;

//! Get test double.
/*!
 * Returns the value of the testDouble variable.
 * \sa testDouble.
 * \return Value of testDouble.
 */
double getTestDouble( ) { return testDouble; }

//! A test value of float-type used to demonstrate the updateMembers() function.
static float testFloat = 1.4f;

//! Get test float.
/*!
 * Returns the value of the testFloat variable.
 * \sa testFloat.
 * \return Value of testFloat.
 */
float getTestFloat( ) { return testFloat; }

//! A test value of int-type used to demonstrate the updateMembers() function.
static int testInt = 2;

//! Get test integer.
/*!
 * Returns the value of the testInt variable.
 * \sa testInt.
 * \return Value of testInt.
 */
int getTestInt( ) { return testInt; }

BOOST_AUTO_TEST_SUITE( test_accelerationModel )

//! Typedef to the function-type used for the verification (double).
/*!
 * Typedef to the function-type used for the verification. Consists of a function which returns a
 * double.
 */
typedef boost::function< double( ) > DoubleReturningFunction;

//! Typedef to the function-type used for the verification (float).
/*!
 * Typedef to the function-type used for the verification. Consists of a function which returns a
 * float.
 */
typedef boost::function< float( ) > FloatReturningFunction;

//! Typedef to the function-type used for the verification (int).
/*!
 * Typedef to the function type used for the verification. Consists of a function which returns an
 * integer.
 */
typedef boost::function< int( ) > IntReturningFunction;

//! (Dummy) 3D acceleration model class.
/*!
 * This class outputs a 3-dimensional acceleration, expressed using a vector of 3 doubles. The
 * updateMembers() function changes one of the vector's members.
 */
class DummyAccelerationModel3d : public AccelerationModel< 3, double >
{
public:

    //! Default constructor.
    /*!
     * Default constructor that takes a boost::function as input, which points to a function that
     * returns a test value. Internally, the updateMembers() function also gets called to ensure
     * that all members are up-to-date after construction.
     * \param dummyFunction Boost-function that points to dummy function that does nothing.
     */
    DummyAccelerationModel3d( DoubleReturningFunction dummyFunction )
        : dummyFunction_( dummyFunction )
    {
        updateMembers( );
    }

    //! Get acceleration.
    /*!
     * Returns acceleration. In this case, this functions simply returns a dummy acceleration
     * vector that is constructed internally.
     * \return Dummy acceleration vector.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return Eigen::Vector3d( currentTestVariable_, 0.0, 0.0 );
    }

    bool updateMembers( )
    {
        currentTestVariable_ = dummyFunction_( );
        return true;
    }

protected:

private:

    //! Pointer to dummy function.
    /*!
     * Boost-function that points to dummy function passed through constructor.
     * \sa DummyAccelerationModel3d().
     */
    DoubleReturningFunction dummyFunction_;

    //! Current test variable.
    /*!
     * Current value of test variable. This value is updated using the updateMembers() function.
     * \sa updateMembers().
     */
    double currentTestVariable_;
};

//! (Dummy) 3D acceleration model class.
/*!
 * This class outputs a 2-dimensional acceleration, expressed using a vector of 2 floats. The
 * updateMembers() function changes one of the vector's members.
 */
class DummyAccelerationModel2f: public AccelerationModel< 2, float >
{
public:

    //! Default constructor.
    /*!
     * Default constructor that takes a boost::function as input, which points to a function that
     * returns a test value. Internally, the updateMembers() function also gets called to ensure
     * that all members are up-to-date after construction.
     * \param dummyFunction Boost-function that points to dummy function that does nothing.
     */
    DummyAccelerationModel2f( FloatReturningFunction dummyFunction )
        : dummyFunction_( dummyFunction )
    {
        updateMembers( );
    }

    //! Get acceleration.
    /*!
     * Returns acceleration. In this case, this functions simply returns a dummy acceleration
     * vector that is constructed internally.
     * \return Dummy acceleration vector.
     */
    Eigen::Vector2f getAcceleration( ) { return Eigen::Vector2f( currentTestVariable_, 0.0 ); }

    //! Update members.
    /*!
     * Updates class members. This function ensures that the getAcceleration() function uses the
     * latest available values in computing an acceleration vector. In this case, the
     * dummyFunction_() is called, and the value is internally stored. The function returns true to
     * indicate that the update was successful.
     * \return Flag indicating if update was successful.
     */
    bool updateMembers( )
    {
        currentTestVariable_ = dummyFunction_( );
        return true;
    }

protected:

private:

    //! Pointer to dummy function.
    /*!
     * Boost-function that points to dummy function passed through constructor.
     * \sa DummyAccelerationModel3d().
     */
    FloatReturningFunction dummyFunction_;

    //! Current test variable.
    /*!
     * Current value of test variable. This value is updated using the updateMembers() function.
     * \sa updateMembers().
     */
    float currentTestVariable_;
};

//! (Dummy) 3D acceleration model class.
/*!
 * This class outputs a 1-dimensional acceleration, expressed using 1 integer. The
 * updateMembers() function changes it's value.
 */
class DummyAccelerationModel1i: public AccelerationModel< 1, int >
{
public:

    //! Default constructor.
    /*!
     * Default constructor that takes a boost::function as input, which points to a function that
     * returns a test value. Internally, the updateMembers() function also gets called to ensure
     * that all members are up-to-date after construction.
     * \param dummyFunction Boost-function that points to dummy function that does nothing.
     */
    DummyAccelerationModel1i(  IntReturningFunction dummyFunction )
        : dummyFunction_( dummyFunction )
    {
        updateMembers( );
    }

    //! Get acceleration.
    /*!
     * Returns acceleration. In this case, this functions simply returns a scalar acceleration.
     * \return Dummy acceleration vector.
     */
    int getAcceleration( ) { return currentTestVariable_; }

    //! Update members.
    /*!
     * Updates class members. This function ensures that the getAcceleration() function uses the
     * latest available values in computing an acceleration. In this case, the dummyFunction_()
     * is called, and the value is internally stored. The function returns true to
     * indicate that the update was successful.
     * \return Flag indicating if update was successful.
     */
    bool updateMembers( )
    {
        currentTestVariable_ = dummyFunction_( );
        return true;
    }

protected:

private:

    //! Pointer to dummy function.
    /*!
     * Boost-function that points to dummy function passed through constructor.
     * \sa DummyAccelerationModel3d().
     */
    IntReturningFunction dummyFunction_;

    //! Current test variable.
    /*!
     * Current value of test variable. This value is updated using the updateMembers() function.
     * \sa updateMembers().
     */
    int currentTestVariable_;
};

//! Test whether AccelerationModel derived classes function correctly.
BOOST_AUTO_TEST_CASE( test_derived3dAccelerationModel )
{
    // Create dummy acceleration model and pass getTestDouble() function as input.
    boost::shared_ptr< AccelerationModel< 3, double > > dummyAccelerationModel
            = boost::make_shared< DummyAccelerationModel3d >( getTestDouble );

    // Get acceleration vector before members are updated.
    Eigen::Vector3d accelerationBeforeUpdate = dummyAccelerationModel->getAcceleration( );

    // Update test double value.
    testDouble = 2.1;

    // Update acceleration model members. This should now result in an update of the internally
    // stored test double value.
    dummyAccelerationModel->updateMembers( );

    // Get acceleration vector, now after members have been updated.
    Eigen::Vector3d accelerationAfterUpdate = dummyAccelerationModel->getAcceleration( );

    // Update test double value.
    testDouble = -87.685;

    // Update and get acceleration with single function.
    Eigen::Vector3d accelerationAfterSecondUpdate = updateAndGetAcceleration(
                dummyAccelerationModel );

    // Check that the acceleration vectors before and after the update match expected values.
    BOOST_CHECK_EQUAL( accelerationBeforeUpdate( 0 ), 3.2 );
    BOOST_CHECK_EQUAL( accelerationAfterUpdate( 0 ), 2.1 );
    BOOST_CHECK_EQUAL( accelerationAfterSecondUpdate( 0 ), -87.685 );
}

BOOST_AUTO_TEST_CASE( test_derived2fAccelerationModel )
{
    // Create dummy acceleration model and pass getTestFloat() function as input.
    boost::shared_ptr< AccelerationModel< 2, float > > dummyAccelerationModel
            = boost::make_shared< DummyAccelerationModel2f >( getTestFloat );

    // Get acceleration vector before members are updated.
    Eigen::Vector2f accelerationBeforeUpdate = dummyAccelerationModel->getAcceleration( );

    // Update test float value.
    testFloat = -7.6f;

    // Update acceleration model members. This should now result in an update of the internally
    // stored test float value.
    dummyAccelerationModel->updateMembers( );

    // Get acceleration vector, now after members have been updated.
    Eigen::Vector2f accelerationAfterUpdate = dummyAccelerationModel->getAcceleration( );

    // Update test float value.
    testFloat = 10.65f;

    // Update and get acceleration with single function.
    Eigen::Vector2f accelerationAfterSecondUpdate = updateAndGetAcceleration(
                dummyAccelerationModel );

    // Check that the acceleration vectors before and after the update match expected values.
    // NOTE: It is imperative that the expected values are given with the "f" suffix, to ensure
    // that their precision is indeed of float-type.
    BOOST_CHECK_EQUAL( accelerationBeforeUpdate( 0 ), 1.4f );
    BOOST_CHECK_EQUAL( accelerationAfterUpdate( 0 ), -7.6f );
    BOOST_CHECK_EQUAL( accelerationAfterSecondUpdate( 0 ), 10.65f );
}

BOOST_AUTO_TEST_CASE( test_derived1iAccelerationModel )
{
    // Create dummy acceleration model and pass getTestInt() function as input.
    boost::shared_ptr< AccelerationModel< 1, int > > dummyAccelerationModel
            = boost::make_shared< DummyAccelerationModel1i >( getTestInt );

    // Get acceleration (scalar) before members are updated.
    int accelerationBeforeUpdate = dummyAccelerationModel->getAcceleration( );

    // Update test int value.
    testInt = 4;

    // Update acceleration model members. This should now result in an update of the internally
    // stored test int value.
    dummyAccelerationModel->updateMembers( );

    // Get acceleration vector, now after members have been updated.
    int accelerationAfterUpdate = dummyAccelerationModel->getAcceleration( );

    // Update test int value.
    testInt = -99;

    // Update and get acceleration with single function.
    int accelerationAfterSecondUpdate = updateAndGetAcceleration( dummyAccelerationModel );

    // Check that the acceleration vectors before and after the update match expected values.
    BOOST_CHECK_EQUAL( accelerationBeforeUpdate, 2 );
    BOOST_CHECK_EQUAL( accelerationAfterUpdate, 4 );
    BOOST_CHECK_EQUAL( accelerationAfterSecondUpdate, -99 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
