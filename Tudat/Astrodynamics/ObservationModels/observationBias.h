#ifndef OBSERVATIONBIAS_H
#define OBSERVATIONBIAS_H

#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_models
{


struct ObservationBiasInterface
{
public:
    ObservationBiasInterface( const int observationSize ): observationSize_( observationSize )
    {
        if( observationSize_ > 0 )
        {
            observationBias_ = Eigen::VectorXd::Constant( observationSize_, 0.0 );
        }
    }

    ObservationBiasInterface( const Eigen::VectorXd observationBias ):
        observationBias_( observationBias ), observationSize_( observationBias.rows( ) ){ }

    virtual ~ObservationBiasInterface( ){ }

    virtual Eigen::VectorXd getObservationBias( const std::vector< double >& observationTimes )
    {
        return observationBias_;
    }

    Eigen::VectorXd getConstantObservationBias( )
    {
        return observationBias_;
    }

    void setConstantObservationBias( const Eigen::VectorXd& observationBias )
    {
        if( observationBias.rows( ) != observationSize_ )
        {
            std::cerr<<"Error when resetting observation bias, sizes are incompatible"<<std::endl;
        }
        observationBias_ = observationBias;
    }

    int getObservationSize( )
    {
        return observationSize_;
    }

protected:

    Eigen::VectorXd observationBias_;

    int observationSize_;
};

struct RangeObservationBiasInterface: public ObservationBiasInterface
{
public:
    RangeObservationBiasInterface(
            const std::map< int, boost::shared_ptr< hardware_models::TimingSystem > > timingSystems, // Boolean denotes whether to add or subtract from range
            const std::vector< bool > addContributions,
            const bool useStochasticError,
            const Eigen::VectorXd observationBias ):ObservationBiasInterface( observationBias ),
        timingSystems_( timingSystems ), addContributions_( addContributions ), useStochasticError_( useStochasticError ){ }

    ~RangeObservationBiasInterface( ){ }

    Eigen::VectorXd getObservationBias( const std::vector< double >& observationTimes );

    void resetUseStochasticError( const bool useStochasticError )
    {
        useStochasticError_ = useStochasticError;
    }

private:
    std::map< int, boost::shared_ptr< hardware_models::TimingSystem > > timingSystems_;

    std::map< int, boost::shared_ptr< hardware_models::TimingSystem > >::iterator timingSystemIterator_;

    std::vector< bool > addContributions_;

    bool useStochasticError_;
};

extern std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationBiasInterface > > > observationBiasInterfaceList;

void setObservationBiases( const std::map< ObservableType, std::map< LinkEnds, double > >& biasesForDoubleObservables );

void setObservationBiases( const std::map< ObservableType, std::vector< LinkEnds > >& linkEndList,
                           const std::map< ObservableType, double >& biasPerObservables );

void setObservationBiases( const std::map< ObservableType, std::vector< LinkEnds > >& linkEndList,
                           const double biasForAllObservations );

std::vector< bool > getClockErrorAdditions( const ObservableType observableType, const LinkEnds& linkEnds );

void setObservationBiasesWithTimingErrors( const std::map< ObservableType, std::vector< LinkEnds > >& linkEndList,
                                           const NamedBodyMap& bodyMap,
                                           const bool includeStochasticNoise = 0,
                                           const std::map< ObservableType, std::map< LinkEnds, double > >& biasesForDoubleObservables =
        std::map< ObservableType, std::map< LinkEnds, double > >( ) );

}

}
#endif // OBSERVATIONMODEL_H
