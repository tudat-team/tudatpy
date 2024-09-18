#ifndef EINSTEININFELDHOFFMANNEQUATIONS_H
#define EINSTEININFELDHOFFMANNEQUATIONS_H

#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace relativity
{

class EinsteinInfeldHoffmannEquations
{
public:
    EinsteinInfeldHoffmannEquations( const std::vector< std::string > acceleratedBodies,
                                     const std::vector< std::string > acceleratingBodies,
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

    Eigen::Vector3d& getRelativePositions( const int bodyUndergoing, const int bodyExerting )
    {
        return currentRelativePositions_.at( bodyUndergoing ).at( bodyExerting );
    }

    double getRelativeDistance( const int bodyUndergoing, const int bodyExerting )
    {
        return currentRelativeDistances_.at( bodyUndergoing ).at( bodyExerting );
    }

    double getInverseSquareDistance( const int bodyUndergoing, const int bodyExerting )
    {
        return currentInverseSquareDistances_.at( bodyUndergoing ).at( bodyExerting );
    }


    Eigen::Vector3d& getRelativeVelocity( const int bodyUndergoing, const int bodyExerting )
    {
        return currentRelativeVelocities_.at( bodyUndergoing ).at( bodyExerting );
    }



    Eigen::Vector3d& getVelocity( const int bodyIndex )
    {
        return currentVelocities_.at( bodyIndex );
    }

    double getGravitationalParameter( const int bodyIndex )
    {
        return currentGravitationalParameters_.at( bodyIndex );
    }



    double getTotalScalarTermCorrection( const int bodyUndergoing, const int bodyExerting )
    {
        return totalScalarTermCorrection.at( bodyUndergoing ).at( bodyExerting );
    }

    Eigen::Vector3d& getTotalVectorTermCorrection( const int bodyUndergoing, const int bodyExerting )
    {
        return totalVectorTermCorrection_.at( bodyUndergoing ).at( bodyExerting );
    }




    double getSingleSourceLocalPotential( const int bodyUndergoing, const int bodyExerting )
    {
        return currentSingleSourceLocalPotential_.at( bodyUndergoing ).at( bodyExerting );
    }

    Eigen::Vector3d& getSinglePointMassAccelerations( const int bodyUndergoing, const int bodyExerting )
    {
        return singlePointMassAccelerations_.at( bodyUndergoing ).at( bodyExerting );
    }

    double getLineOfSighSpeed( const int bodyUndergoing, const int bodyExerting )
    {
        return lineOfSightSpeed_.at( bodyUndergoing ).at( bodyExerting );
    }

    double getLocalPotential( const int bodyIndex )
    {
        return currentLocalPotentials_.at( bodyIndex );
    }

    Eigen::Vector3d& getTotalPointMassAcceleration( const int bodyIndex )
    {
        return totalPointMassAccelerations_.at( bodyIndex );
    }

    std::vector< std::vector< std::vector< double > > > getScalarEihCorrections( )
    {
        return scalarEihCorrections_;
    }

    double getScalarEihCorrection( const int k, const int bodyUndergoing, const int bodyExerting )
    {
        return scalarEihCorrections_.at( k ).at( bodyUndergoing ).at( bodyExerting );
    }


    std::vector< std::vector< std::vector< Eigen::Vector3d > > > getVectorEihCorrections( )
    {
        return vectorEihCorrections_;
    }

    Eigen::Vector3d getVectorEihCorrection( const int k, const int bodyUndergoing, const int bodyExerting )
    {
        return vectorEihCorrections_.at( k ).at( bodyUndergoing ).at( bodyExerting );
    }

    std::vector< std::string > getBodiesUndergoingAcceleration( )
    {
        return acceleratedBodies_;
    }

    std::vector< std::string > getBodiesExertingAcceleration( )
    {
        return acceleratingBodies_;
    }

    double getScalarTermMultiplier( const int k )
    {
        return scalarTermMultipliers_.at( k );
    }

    double getVectorTermMultiplier( const int k )
    {
        return vectorTermMultipliers_.at( k );
    }

    bool getOmitMainTerm( )
    {
        return omitMainTerm_;
    }

    void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }

    std::map< std::string, int > getAcceleratedBodyMap( )
    {
        return acceleratedBodyMap_;
    }

    std::map< std::string, int > getAcceleratingBodyMap( )
    {
        return acceleratingBodyMap_;
    }


    void recomputeExpansionMultipliers( );


private:


    void calculateAccelerations( );

    std::vector< std::string > acceleratedBodies_;

    std::vector< std::string > acceleratingBodies_;

    std::vector< std::function< double( ) > > gravitationalParameterFunction_;

    std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > bodyStateFunctions_;

    std::function< double( ) > ppnGammaFunction_;

    std::function< double( ) > ppnBetaFunction_;

    bool omitMainTerm_;




    std::map< std::string, int > acceleratedBodyMap_;

    std::map< std::string, int > acceleratingBodyMap_;


    // mu_{i}
    std::vector< double > currentGravitationalParameters_;

    // v_{i}
    std::vector< Eigen::Vector3d > currentVelocities_;

    // r_{i}
    std::vector< Eigen::Vector3d > currentPositions_;

    // v_{i} * v_{i}
    std::vector< double > currentSquareSpeeds_;



    // sum_(j not i) ( mu_j / ||r_{ij}|| ) = sum( currentSingleSourceLocalPotential_ )
    std::vector< double > currentLocalPotentials_;

    // sum_(j not i) ( mu_{j} * r_{ij} / ||r_{ij}||^3 ) = sum( singlePointMassAccelerations_ )
    std::vector< Eigen::Vector3d > totalPointMassAccelerations_;



    // r_{ij} = r_{j} - r_{i}
    std::vector< std::vector< Eigen::Vector3d > > currentRelativePositions_;

    // v_{ij} = v_{j} - v_{i}
    std::vector< std::vector< Eigen::Vector3d > > currentRelativeVelocities_;

    // || r_{ij} ||
    std::vector< std::vector< double > > currentRelativeDistances_;

    // 1  / || r_{ij} ||^2
    std::vector< std::vector< double > > currentInverseSquareDistances_;

    // r_{ij} * v_{j}
    std::vector< std::vector< double > > lineOfSightSpeed_;

    // v_{i} * v_{j}
    std::vector< std::vector< double > > velocityInnerProducts_;

    // mu_j / ||r_{ij}||
    std::vector< std::vector< double > > currentSingleSourceLocalPotential_;

    // mu_{j} * r_{ij} / ||r_{ij}||^3
    std::vector< std::vector< Eigen::Vector3d > > singlePointMassAccelerations_;



    std::vector< std::vector< Eigen::Vector3d > > currentSingleAccelerations_;

    std::vector< std::vector< double > > totalScalarTermCorrection;

    std::vector< std::vector< Eigen::Vector3d > > totalVectorTermCorrection_;

    std::vector< Eigen::Vector3d > currentAccelerations_;



    std::vector< std::vector< std::vector< double > > > scalarEihCorrections_;

    std::vector< std::vector< std::vector< Eigen::Vector3d > > > vectorEihCorrections_;


    double currentPpnGamma_;

    double currentPpnBeta_;


    std::vector< double > scalarTermMultipliers_;

    std::vector< double > vectorTermMultipliers_;

    double currentTime_;
};

}

}

#endif // EINSTEININFELDHOFFMANNEQUATIONS_H
