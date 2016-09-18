#ifndef TUDAT_SPHERICALHARMONICSINECOEFFICIENTS_H
#define TUDAT_SPHERICALHARMONICSINECOEFFICIENTS_H

#include <map>

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Class for spherical harmonic sine coefficient(s)
/*!
 *  Class for spherical harmonic sine coefficient(s). Specific parameter is a block from the sine coefficient matrix of SphericalHarmonicGravityField,
 *  with size set by constructor.
 */
class SphericalHarmonicsSineCoefficients: public EstimatableParameter< Eigen::VectorXd >
{

public:

    SphericalHarmonicsSineCoefficients(
            const boost::function< Eigen::MatrixXd( ) > getSineCoefficients,
            const boost::function< void( Eigen::MatrixXd ) > setSineCoefficients,
            const int minimumDegree, const int minimumOrder,
            const int maximumDegree, const int maximumOrder,
            const std::string& associatedBody );

    ~SphericalHarmonicsSineCoefficients( ) { }

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( const Eigen::VectorXd parameterValue );

    int getParameterSize( ) { return parameterSize_; }

    std::map< int, std::pair< int, int > > getBlockIndices( )
    {
        return blockIndices_;
    }


protected:

private:

    Eigen::MatrixXd convertSineCoefficientParameterToCoefficientBlock( Eigen::VectorXd parameterVector );


    boost::function< Eigen::MatrixXd( ) > getSineCoefficients_;
    boost::function< void( Eigen::MatrixXd ) > setSineCoefficients_;

    int minimumDegree_;

    int minimumOrder_;

    int maximumDegree_;

    int maximumOrder_;

    std::map< int, std::pair< int, int > > blockIndices_;// degree , (start order, number of entries)

    int parameterSize_;
};

}

}

#endif // SPHERICALHARMONICSINECOEFFICIENTS_H
