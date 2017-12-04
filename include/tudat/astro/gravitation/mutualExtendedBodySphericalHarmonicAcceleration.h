#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
#define MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Basics/basicTypedefs.h>
#include <Tudat/Mathematics/BasicMathematics/legendrePolynomials.h>

#include "Tudat/Astrodynamics/Gravitation/mutualForcePotential.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"

namespace tudat
{

namespace gravitation
{

class MutualExtendedBodySphericalHarmonicAcceleration: public basic_astrodynamics::AccelerationModel3d
{

public:

    MutualExtendedBodySphericalHarmonicAcceleration(
            const boost::function< Eigen::Vector3d( ) > positionOfBody1Function,
            const boost::function< Eigen::Vector3d( ) > positionOfBody2Function,
            const boost::function< double( ) > gravitationalParameterFunction,
            const double equatorialRadiusOfBody1,
            const double equatorialRadiusOfBody2,
            const boost::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody1Function,
            const boost::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody1Function,
            const boost::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody2Function,
            const boost::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody2Function,
            const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
            const boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation,
            const boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation,
            const bool useCentraBodyFrame,
            const bool areCoefficientsNormalized = 1 );

    void updateMembers( const double currentTime = TUDAT_NAN );

    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }


    bool getUseCentraBodyFrame( )
    {
        return useCentraBodyFrame_;
    }

    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    Eigen::Vector3d getMutualPotentialGradient( )
    {
        return mutualPotentialGradient_;
    }

    Eigen::Quaterniond getCurrentRotationFromBody2ToBody1( )
    {
        return currentRotationFromBody2ToBody1_;
    }

    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > getEffectiveMutualPotentialField( )
    {
        return effectiveMutualPotentialField_;
    }

    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > getSphericalHarmonicsCache( )
    {
        return sphericalHarmonicsCache_;
    }

    double getEquatorialRadiusOfBody1( )
    {
        return equatorialRadiusOfBody1_;
    }

    double getEquatorialRadiusOfBody2( )
    {
        return equatorialRadiusOfBody2_;
    }

    double getRadius1Power( const int index )
    {
        return radius1Powers_.at( index );
    }

    double getRadius2Power( const int index )
    {
        return radius2Powers_.at( index );
    }

    std::vector< double > getRadius1Powers( )
    {
        return radius1Powers_;
    }

    std::vector< double > getRadius2Powers( )
    {
        return radius2Powers_;
    }

private:

    //! Function returning the position of body 1
    boost::function< Eigen::Vector3d( ) > positionOfBody1Function_;

    //! Function returning the position of body 2
    boost::function< Eigen::Vector3d( ) > positionOfBody2Function_;

    //! Function returning the effective gravitational parameter
    boost::function< double( ) > gravitationalParameterFunction_;

    //! Equatorial radius of body 1.
    double equatorialRadiusOfBody1_;

    //! Equatorial radius of body 2.
    double equatorialRadiusOfBody2_;

    //! List of degrees/orders that are to be used for the series expansion.
    /*!
     * List of degrees/orders that are to be used for the series expansion. Each tuple contains: (degree of body 1, order of
     * body 1, degree of body 2, order of body 2)
     */
    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    //! Function that returns the effective one-body spherical harmonic cosine coefficient
    /*!
     *  Function that returns the effective one-body spherical harmonic cosine coefficient as a function of (degree of body 1,
     *  order of body 1, degree of body 2, order of body 2)
     */
    boost::function< double( int, int, int, int ) > effectiveCosineCoefficientFunction_;

    //! Function that returns the effective one-body spherical harmonic sine coefficient
    /*!
     *  Function that returns the effective one-body spherical harmonic sine coefficient as a function of (degree of body 1,
     *  order of body 1, degree of body 2, order of body 2)
     */
    boost::function< double( int, int, int, int ) > effectiveSineCoefficientFunction_;

    //! Function that returns the rotation from the inertial to body-fixed frame of body 1
    boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation_;

    //! Function that returns the rotation from the inertial to body-fixed frame of body 2
    boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation_;

    //! Maximum degree of effective one-body spherical harmonic expansion
    unsigned int maximumDegree_;

    //! Maximum order of effective one-body spherical harmonic expansion
    unsigned int maximumOrder_;

    //! Object that computes the effect one-body coefficients from the two one-body gravoity fields and their relative orientation
    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;

    //! Function that computes the coefficients of body, in the body-fixed frame of body 1
    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    //! Acceleration, as computed by last call to updateMembers
    Eigen::Vector3d currentAcceleration_;

    //! Position of body 1 w.r.t. body 2, expressed in inertial frame, as computed by last call to updateMembers
    Eigen::Vector3d currentRelativePosition_;

    //! Position of body 1 w.r.t. body 2, expressed in body-fixed frame of body 1, as computed by last call to updateMembers
    Eigen::Vector3d currentBodyFixedRelativePosition_;

    //! Current gradient of mutual potential, in body-fixed frame of body 1 in Cartesian coordinates
    Eigen::Vector3d mutualPotentialGradient_;

    //! Current rotation from body-fixed frame of body 2 to that of body 1, as computed by last call to updateMembers
    Eigen::Quaterniond currentRotationFromBody2ToBody1_;

    //! Variable denoting whether acceleration is expressed in body-fixed frame of body 1, or in inertial frame.
    bool useCentraBodyFrame_;

    //! Boolean denoting whether spherical harmonic coefficients are normalized
    bool areCoefficientsNormalized_;

    //! List of integer powers of (distance / equatorial radius of body 1)
    std::vector< double > radius1Powers_;

    //! List of integer powers of (distance / equatorial radius of body 2)
    std::vector< double > radius2Powers_;


};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
