#ifndef ROTATIONMATRIXPARTIAL_H
#define ROTATIONMATRIXPARTIAL_H

#include <vector>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

namespace tudat
{

namespace observation_partials
{


Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
        const Eigen::Quaterniond intertialBodyFixedToIntegrationFrame,
        const double rotationRate, const double timeSinceEpoch );

std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,//right ascension, declination, longitude of prime meridian.
        const double rotationRate, const double timeSinceEpoch );

//! Base class for partial derivatives of rotation matrix from body fixed to inertial frame w.r.t. a parameter.
/*!
 *  Base class for partial derivatives of rotation matrix from body fixed to inertial frame w.r.t. a parameter. A derived class is implemented
 *  for each estimatable parameter which represents a property of a rotation matrix from a body-fixed frame (as defined by the
 *  isParameterRotationMatrixProperty function). Pointers to this class are used for both PositionPartials and the SphericalHarmonicsGravityPartial
 *  partial. In this manner, only a single partial object needs to be implemented for both the position and partial and the sh acceleration
 *  partial wrt a rotational parameter when such a parameter is added.
 */
class RotationMatrixPartial
{
public:

    //! Virtual destructor
    /*!
     *  Virtual destructor.
     */
    virtual ~RotationMatrixPartial( ){ }

    //! Function to calculate partials of rotation matrix wrt a parameter (for either double of Vector parameters).
    /*!
     *  Function to calculate partials of rotation matrix from the body-fixed to the inertial base frame wrt a parameter.
     *  For a double parameter, the size of the return vector is 1, for a VectorXd parameter, the size is equal to the size of the parameter
     *  and the entries of the vector returned from here correspond to partials of the same index in the parameter.
     *  \param time Time at which the partial(s) is(are) to be evaluated.
     *  \return Vector of partials of rotation matrix from body-fixed to inertial frame wrt parameter(s)
     */
    virtual std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time ) = 0;

    //! Function to calculate the partial of the position of a vector, which is given in a body-fixed frame, in the inertial frame wrt a parameter.
    /*!
     *  Function to calculate the partial of the position of a vector, which is given in a body-fixed frame, in the inertial frame wrt a parameter
     *  denoting a property of the rotation between a body-fixed and an inertial frame. The type of parameter is defined by the derived class and a
     *  a separate derived class is impleneted for each parameter. An example is the change in position of a ground station,
     *  as expressed in the inertial frame, wrt a rotational parameter. Each column of the return Matrix denotes a single entry of the parameter
     *  (so it is a Vector3d for a double parameter).
     *  \param time Time at which the partial is to be evaluated.
     *  \param vectorInLocalFrame Vector, expressed in the body-fixed frame, of which the partial in an inertial frame wrt parameter is to be determined.
     *  \return Partial of the value of the vector in an inertial frame wrt the parameter(s).
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfRotatedVector(
            const double time,
            const Eigen::Vector3d vectorInLocalFrame );

    virtual std::string getSecondaryIdentifier( )
    {
        return "";
    }

};

class RotationMatrixPartialWrtConstantRotationRate: public RotationMatrixPartial
{
public:

    RotationMatrixPartialWrtConstantRotationRate(
            const boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel ):
        bodyRotationModel_( bodyRotationModel ){ }

    ~RotationMatrixPartialWrtConstantRotationRate( ){ }

    //R^(I/B)
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        return boost::assign::list_of( calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
                                           bodyRotationModel_->getInitialRotationToTargetFrame( ).inverse( ),
                                           bodyRotationModel_->getRotationRate( ), time - bodyRotationModel_->getInitialSecondsSinceEpoch( ) ) );
    }
private:

    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel_;
};

class RotationMatrixPartialWrtPoleOrientation: public RotationMatrixPartial
{
public:

    RotationMatrixPartialWrtPoleOrientation( const boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel ):
        bodyRotationModel_( bodyRotationModel ){ }

    ~RotationMatrixPartialWrtPoleOrientation( ){ }

    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        return calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
                    bodyRotationModel_->getInitialEulerAngles( ),
                    bodyRotationModel_->getRotationRate( ), time - bodyRotationModel_->getInitialSecondsSinceEpoch( ) );
    }

private:

    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel_;
};


}

}

#endif
