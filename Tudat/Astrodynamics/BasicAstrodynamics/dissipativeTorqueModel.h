#ifndef DISSIPATIVETORQUEMODEL_H
#define DISSIPATIVETORQUEMODEL_H

#include <iomanip>
#include <functional>

#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{

namespace basic_astrodynamics
{

class DissipativeTorqueModel: public TorqueModel
{
public:

    DissipativeTorqueModel(
            const std::function< Eigen::Vector3d( ) >& bodyFixedRotationVectorFunction,
            const std::function< Eigen::Matrix3d( ) >& dampingMatrixFunction,
            const double bodyMeanRotationRate ):
        bodyFixedRotationVectorFunction_( bodyFixedRotationVectorFunction ),
        dampingMatrixFunction_( dampingMatrixFunction ),
    bodyMeanRotationRate_( bodyMeanRotationRate ){ }

    ~DissipativeTorqueModel( ) { }

    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    void updateMembers( const double currentTime )
    {
        currentBodyRotationPerturbationVector_ = bodyFixedRotationVectorFunction_( );

        currentBodyRotationPerturbationVector_( 2 ) -= bodyMeanRotationRate_;
        currentTorque_ =  -dampingMatrixFunction_( ) * currentBodyRotationPerturbationVector_;
    }

    void setDampingMatrixFunction( const Eigen::Matrix3d& dampingMatrix )
    {
        dampingMatrixFunction_ = boost::lambda::constant( dampingMatrix );
    }

protected:

private:

    std::function< Eigen::Vector3d( ) > bodyFixedRotationVectorFunction_;

    std::function< Eigen::Matrix3d( ) > dampingMatrixFunction_;

    double bodyMeanRotationRate_;

    Eigen::Vector3d currentBodyRotationPerturbationVector_;

    Eigen::Vector3d currentTorque_;

};

}

}

#endif // DISSIPATIVETORQUEMODEL_H
