#ifndef TUDAT_SPHERICALHARMONICCOSINECOEFFICIENTS_H
#define TUDAT_SPHERICALHARMONICCOSINECOEFFICIENTS_H

#include <map>

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{


std::map< int, std::pair< int, int > > convertShDegreeAndOrderLimitsToBlockIndices(
        const int minimumDegree, const int minimumOrder, const int maximumDegree, const int maximumOrder,
        int& totalSetSize );

std::map< int, std::pair< int, int > > convertShDegreeAndOrderLimitsToBlockIndices(
        const int minimumDegree, const int minimumOrder, const int maximumDegree, const int maximumOrder );

class SphericalHarmonicsCosineCoefficients: public EstimatableParameter< Eigen::VectorXd >
{

public:
    SphericalHarmonicsCosineCoefficients(
            boost::function< Eigen::MatrixXd( ) > getCosineCoefficients,
            boost::function< void( Eigen::MatrixXd ) > setCosineCoefficients,
            const int minimumDegree, const int minimumOrder,
            const int maximumDegree, const int maximumOrder,
            std::string& associatedBody );

    ~SphericalHarmonicsCosineCoefficients( ) { }

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( const Eigen::VectorXd parameterValue );

    int getParameterSize( ) { return parameterSize_; }

    std::map< int, std::pair< int, int > > getBlockIndices( ){ return blockIndices_; }

    int getMnimumDegree( )
    {
        return minimumDegree_;
    }

    int getMinimumOrder( )
    {
        return minimumOrder_;
    }
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }
    int getMaximumOrder( )
    {
        return maximumOrder_;
    }

protected:

private:
    boost::function< Eigen::MatrixXd( ) > getCosineCoefficients_;
    boost::function< void( Eigen::MatrixXd ) > setCosineCoefficients_;

    int minimumDegree_;
    int minimumOrder_;
    int maximumDegree_;
    int maximumOrder_;

    std::map< int, std::pair< int, int > > blockIndices_;// degree , (start order, number of entries)

    int parameterSize_;
};

}

}

#endif // SPHERICALHARMONICCOSINECOEFFICIENTS_H
