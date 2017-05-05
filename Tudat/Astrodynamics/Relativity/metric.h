#ifndef TUDAT_METRIC_H
#define TUDAT_METRIC_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace relativity
{

const static Eigen::Matrix4d minkowskiMetric = ( Eigen::Matrix4d( ) <<
                                                 -1.0, 0.0, 0.0, 0.0,
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 0.0, 1.0 ).finished( );

struct PPNParameterSet
{
public:
    PPNParameterSet( const double parameterGamma, const double parameterBeta ):
        parameterGamma_( parameterGamma ), parameterBeta_( parameterBeta )
    { }

    ~PPNParameterSet( ){ }

    double getParameterGamma( )
    {
        return parameterGamma_;
    }

    double getParameterBeta( )
    {
        return parameterBeta_;
    }


    void setParameterGamma( const double parameterGamma )
    {
        parameterGamma_ = parameterGamma;
    }

    void setParameterBeta( const double parameterBeta )
    {
        parameterBeta_ = parameterBeta;
    }

protected:
    double parameterGamma_;

    double parameterBeta_;

};

extern boost::shared_ptr< PPNParameterSet > ppnParameterSet;


}

}

#endif // TUDAT_METRIC_H
