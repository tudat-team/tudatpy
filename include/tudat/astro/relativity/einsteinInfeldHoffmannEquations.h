#ifndef EINSTEININFELDHOFFMANNEQUATIONS_H
#define EINSTEININFELDHOFFMANNEQUATIONS_H

#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace gravitation
{

class EinsteinInfeldHoffmannEquations
{
public:
    EinsteinInfeldHoffmannEquations( const std::vector< std::string > acceleratedBodies,
                                     const std::vector< std::function< double( ) > > gravitationalParameterFunction,
                                     const std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > bodyStateFunctions,
                                     const std::function< double( ) > ppnGammaFunction,
                                     const std::function< double( ) > ppnBetaFunction );

    void update( const double currentTime );

    Eigen::Vector3d getAccelerationOfBody( const int bodyIndex )
    {
        return currentAccelerations_[ bodyIndex ];
    }

    Eigen::Vector3d getAccelerationOfBody( const std::string bodyName )
    {
        return getAccelerationOfBody( acceleratedBodyMap_.at( bodyName ) );
    }

private:

    void calculateAccelerations( );

    std::vector< std::string > acceleratedBodies_;

    std::vector< std::function< double( ) > > gravitationalParameterFunction_;

    std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > bodyStateFunctions_;

    std::function< double( ) > ppnGammaFunction_;

    std::function< double( ) > ppnBetaFunction_;





    std::map< std::string, int > acceleratedBodyMap_;

    std::vector< double > currentGravitationalParameters_;

    std::vector< Eigen::Vector3d > currentVelocities_;

    std::vector< Eigen::Vector3d > currentPositions_;

    std::vector< double > currentSquareSpeeds_;

    std::vector< std::vector< Eigen::Vector3d > > currentRelativePositions_;

    std::vector< Eigen::Vector3d > totalPointMassAccelerations_;

    std::vector< std::vector< Eigen::Vector3d > > singlePointMassAccelerations_;

    std::vector< std::vector< double > > currentRelativeDistances_;

    std::vector< std::vector< double > > secondaryCentralBodyTerms_; //sum mu_{l}/r_{il} (i=1st, l=2nd)

    std::vector< std::vector< double > > relativePositionVelocityProduct_; //(r_{i}-r_{j})*\dot r_{j} (i=1st, j=2nd)

    std::vector< std::vector< double > > velocityInnerProducts_;

    std::vector< std::vector< Eigen::Vector3d > > currentSingleAccelerations_;

    std::vector< Eigen::Vector3d > currentAccelerations_;


    double currentPpnGamma_;

    double currentPpnBeta_;


    double currentTime_;
};

}

}

#endif // EINSTEININFELDHOFFMANNEQUATIONS_H
