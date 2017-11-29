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

    Eigen::Vector3d getCurrentEulerAngles( )
    {
        return currentEulerAngles_;
    }

    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > getEffectiveMutualPotentialField( )
    {
        return effectiveMutualPotentialField_;
    }

    basic_mathematics::LegendreCache* getLegendreCache( )
    {
        return legendreCache_;
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

    boost::function< Eigen::Vector3d( ) > positionOfBody1Function_;
    boost::function< Eigen::Vector3d( ) > positionOfBody2Function_;

    boost::function< double( ) > gravitationalParameterFunction_;

    double equatorialRadiusOfBody1_;
    double equatorialRadiusOfBody2_;

    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    boost::function< double( int, int, int, int ) > effectiveCosineCoefficientFunction_;
    boost::function< double( int, int, int, int ) > effectiveSineCoefficientFunction_;

    boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation_;
    boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation_;

    int maximumDegree_;
    int maximumOrder_;

    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;
    basic_mathematics::LegendreCache* legendreCache_;

    Eigen::Vector3d currentAcceleration_;
    Eigen::Vector3d currentRelativePosition_;
    Eigen::Vector3d currentBodyFixedRelativePosition_;
    Eigen::Vector3d mutualPotentialGradient_;

    Eigen::Quaterniond currentRotationFromBody2ToBody1_;

    Eigen::Vector3d currentEulerAngles_;

    bool useCentraBodyFrame_;
    bool areCoefficientsNormalized_;

    std::vector< double > radius1Powers_;
    std::vector< double > radius2Powers_;


};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
