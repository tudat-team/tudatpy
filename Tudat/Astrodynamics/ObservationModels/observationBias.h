#ifndef TUDAT_OBSERVATIONBIAS_H
#define TUDAT_OBSERVATIONBIAS_H

#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_models
{

template< int ObservationSize = 1 >
class ObservationBias
{
public:
    ObservationBias( ){ }

    virtual ~ObservationBias( ){ }

    virtual Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& observationTimes ) = 0;

    int getObservationSize( )
    {
        return ObservationSize;
    }
};

template< int ObservationSize = 1 >
class ConstantObservationBias: public ObservationBias< ObservationSize >
{
public:

    ConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > observationBias ):
        observationBias_( observationBias ){ }

    ~ConstantObservationBias( ){ }

    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& observationTimes =  std::vector< double >( ) )
    {
        return observationBias_;
    }

    Eigen::Matrix< double, ObservationSize, 1 > getConstantObservationBias( )
    {
        return observationBias_;
    }

    void setConstantObservationBias( const Eigen::VectorXd& observationBias )
    {
        if( observationBias.rows( ) != ObservationSize )
        {
            std::cerr<<"Error when resetting observation bias, sizes are incompatible"<<std::endl;
        }
        observationBias_ = observationBias;
    }


private:

    Eigen::Matrix< double, ObservationSize, 1 > observationBias_;
};

}

}
#endif // TUDAT_OBSERVATIONMODEL_H
